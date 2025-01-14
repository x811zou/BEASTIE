#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

import logging
import os
import sys
import pandas as pd
from pkg_resources import resource_filename
import gzip
import csv
from .helpers import chrRange, tabix_regions
from .misc_tools.GffTranscriptReader import GffTranscriptReader
from .misc_tools.Pipe import Pipe
import subprocess

def normalize_chromosome_name(chrom):
    """
    Normalize chromosome names by removing the 'chr' prefix (if present)
    and converting to lowercase for consistent comparison.
    """
    return chrom.strip().lstrip("chr").lower()


def chrRange_species_aware(chr_start, chr_end, available_chromosomes):
    """
    Generate a list of target chromosomes based on the range and available chromosomes.
    """
    chroms = [str(i) for i in range(chr_start, chr_end + 1)]
    special_chroms = [chrom for chrom in available_chromosomes if chrom not in chroms]
    chroms.extend(special_chroms)
    return chroms


def get_available_chromosomes_from_vcf(vcf_filename):
    """
    Parse the VCF file to detect available chromosomes.
    Return a dictionary mapping normalized names to original names.
    """
    chromosomes = {}
    with gzip.open(vcf_filename, "rt") as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue  # Skip header lines
            fields = line.split("\t")
            if len(fields) > 0:
                chrom = fields[0].strip()
                normalized_chrom = normalize_chromosome_name(chrom)
                chromosomes[normalized_chrom] = chrom  # Map normalized to original case
    return chromosomes


def get_available_chromosomes(gtf_filename):
    """
    Parse the GTF file to detect available chromosomes.
    Return a dictionary mapping normalized names to original names.
    """
    chromosomes = {}
    with gzip.open(gtf_filename, "rt") as gtf_file:
        for line in gtf_file:
            if line.startswith("#"):
                continue  # Skip header lines
            fields = line.split("\t")
            if len(fields) > 0:
                chrom = fields[0].strip()
                normalized_chrom = normalize_chromosome_name(chrom)
                chromosomes[normalized_chrom] = chrom  # Map normalized to original case
    return chromosomes


def chunk_iter(iter, n):
    """Yield successive n-sized chunks from iter."""
    res = []
    try:
        while True:
            while len(res) < n:
                v = next(iter)
                res.append(v)
            yield res
            res = []
    except StopIteration:
        if len(res) > 0:
            yield res

""" Check if a genotype is heterozygous by testing whether them match with the 6 types of homozygous options
"""
HOMOZYGOUS_GENOTYPES = [
    "0|0", "1|1", "2|2", "3|3", "4|4", "5|5", "6|6", "7|7",
    "0/0", "1/1", "2/2"
]

def isHeterozygous(genotype):
    return not genotype in HOMOZYGOUS_GENOTYPES

def make_vcfline_processor(require_pass):
    def vcfline_processor(line):
        nonlocal require_pass
        fields = line.split("\t")
        if len(fields) < 10:
            logging.info(f"bad field line: '{line}'")
        assert len(fields) >= 10
        if require_pass and fields[6] != "PASS":
            return None
        genotype = fields[9].split(":")[0]

        if not isHeterozygous(genotype):
            return None

        pos = int(fields[1])
        rs = fields[2]

        return (pos, rs, genotype)

    return vcfline_processor

def check(hetSNP, VCF):
    if ((os.path.exists(hetSNP) and os.path.getsize(hetSNP) > 0) and
        (os.path.getmtime(hetSNP) > os.path.getmtime(VCF))):
        logging.info(f"..... Check: {os.path.basename(hetSNP)} exists and non-empty, and is newer than {os.path.basename(VCF)} -- skipping hetSNP file generation")
        return True
    else:
        logging.info(f"..... Check: {os.path.basename(hetSNP)} does not exist or is not newer than {os.path.basename(VCF)} -- generating hetSNP file")
    return False

def count_all_het_sites(
    sample,
    vcfFilename,
    outputFilename,
    chr_start,
    chr_end,
    genecode_gz,
    require_pass,
    DEBUG_GENES=None,
    **kwargs,
):
    logging.basicConfig(
        format="%(asctime)-15s [%(levelname)s] %(message)s",
        level=logging.INFO,
    )

    logging.info(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> extractHets: Start sample {sample} - require PASS :{require_pass}")
    if check(outputFilename, vcfFilename):
        return

    include_x_chromosome = kwargs.pop("include_x_chromosome", False)

    # Step 1: Get available chromosomes from the GTF file
    available_chromosomes = get_available_chromosomes(genecode_gz)
    logging.info(f"Chromosomes parsed from GTF: {available_chromosomes}")

    # Step 2: Get available chromosomes from the VCF
    available_chromosomes_vcf = get_available_chromosomes_from_vcf(vcfFilename)
    # After loading chromosomes from VCF
    logging.info(f"Chromosomes parsed from VCF: {available_chromosomes_vcf}")   

    # Step 3: Get the intersection of chromosomes in GTF and VCF
    common_chromosomes = set(available_chromosomes.keys()).intersection(set(available_chromosomes_vcf.keys()))
    logging.info(f"Common chromosomes for processing: {common_chromosomes}")

    # Step 4: Generate the list of target chromosomes
    # Only include chromosomes that exist in both VCF and GTF
    target_chromosomes = [
        chrom for chrom in chrRange_species_aware(chr_start, chr_end, common_chromosomes)
        if chrom in available_chromosomes_vcf
    ]

    if include_x_chromosome and 'x' in available_chromosomes_vcf:
        target_chromosomes.append('x')

    logging.info(f"Target chromosomes for extraction: {target_chromosomes}")

    out_stream = open(outputFilename, "w")
    out_stream.write(
        "chr\tchrN\tpos\tgeneID\ttranscriptID\ttranscript_pos\tSNP_id\tgenotype\n"
    )

    vcfline_processor = make_vcfline_processor(require_pass)

    reader = GffTranscriptReader()
    logging.info(f"..... Loading GTF file: {genecode_gz}")
    full_geneList = reader.loadGenes(genecode_gz)

    for chrId in target_chromosomes:
        chrom_name = available_chromosomes_vcf[chrId]  # Get the original chromosome name from VCF
        chr_normalized = normalize_chromosome_name(chrom_name)

        # Initialize geneList as an empty list
        geneList = []

        if chr_normalized not in available_chromosomes:
            logging.info(f"Skipping {chrId} as it is not present in the GTF file")
            continue

        # Filter genes for the current chromosome
        geneList = list(
            filter(lambda gene: normalize_chromosome_name(gene.getSubstrate()) == chr_normalized, full_geneList)
        )

        logging.info(f"..... Start loading gencode annotation {chrId} with {len(geneList)} genes")

        exon_region_to_transcripts = {}
        exon_region_to_snp_infos = {}  # Ensure this variable is always initialized

        if geneList:  # Proceed only if geneList is not empty
            for gene in geneList:
                Num_transcript = gene.getNumTranscripts()
                for n in range(Num_transcript):
                    transcript = gene.getIthTranscript(n)
                    rawExons = transcript.getRawExons()
                    for exon in rawExons:
                        begin = exon.getBegin() + 1
                        end = exon.getEnd()
                        exon_region = f"{chrId}:{begin}-{end}"
                        if exon_region not in exon_region_to_transcripts:
                            exon_region_to_transcripts[exon_region] = []
                        exon_region_to_transcripts[exon_region].append(transcript)

            exon_region_to_snp_infos = tabix_regions(
                exon_region_to_transcripts.keys(), vcfline_processor, vcfFilename
            )

            logging.info(f"Extracted SNPs for {chrId}: {len(exon_region_to_snp_infos)} regions")
        else:
            logging.info(f"No genes found for chromosome {chrId}, skipping SNP extraction")


        data = []
        variant_to_transcript_info = {}
        for exon_region in exon_region_to_snp_infos:
            snp_infos = exon_region_to_snp_infos[exon_region]
            transcripts = exon_region_to_transcripts[exon_region]
            for snp_info in snp_infos:
                pos, rsid, genotype = snp_info
                for transcript in transcripts:
                    if normalize_chromosome_name(chrom_name) == normalize_chromosome_name(transcript.getSubstrate()):
                        chr_pos = f"{chrom_name}_{pos}"
                        if chr_pos not in variant_to_transcript_info:
                            variant_to_transcript_info[chr_pos] = []
                        variant_to_transcript_info[chr_pos].append(
                            (transcript, pos, rsid, genotype)
                        )

        for chr_pos in variant_to_transcript_info:
            gene_ids = set([x[0].getGeneId() for x in variant_to_transcript_info[chr_pos]])
            if len(gene_ids) == 1:
                transcript, pos, rsid, genotype = variant_to_transcript_info[chr_pos][0]
                transcriptCoord = transcript.mapToTranscript(int(pos))

                # Use the original chromosome name with 'chr' for the 'chr' column
                chr_original = available_chromosomes[chrId]
                
                # Strip 'chr' from the chromosome name for the 'chrN' column
                chrN = chr_original.lstrip("chr")

                data.append(
                    [chr_original, chrN, pos, transcript.getGeneId(), transcript.getId(), transcriptCoord, rsid, genotype]
                )

        data.sort(key=lambda r: (normalize_chromosome_name(r[1]), r[2]))

        for r in data:
            out_stream.write("\t".join(map(str, r)))
            out_stream.write("\n")
    out_stream.close()
    logging.info(f"..... Finish writing to {outputFilename}")
    logging.info(f"..... Done sample {sample}")

def main():
    if len(sys.argv) != 8:
        print(
            "Usage: extractHets.py <sample> <in_vcfgz> <annotation_gz> <chr_start> <chr_end> <out_hetSNP> <skip_require_pass> <debug_gene>"
        )
        sys.exit(1)
    
    sample = sys.argv[1]
    vcfgz_path_filename = sys.argv[2]
    annotation_gz = sys.argv[3]
    chr_start = int(sys.argv[4])
    chr_end = int(sys.argv[5])
    hetSNP_filename = sys.argv[6]
    skip_require_pass = sys.argv[7]
    debug_gene = None
    count_all_het_sites(sample, vcfgz_path_filename, hetSNP_filename, chr_start, chr_end, annotation_gz, not skip_require_pass, debug_gene)

if __name__ == "__main__":
    main()
