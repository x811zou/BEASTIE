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
    Remove 'chr' prefix and any leading/trailing whitespace from chromosome names.
    """
    return chrom.strip().lstrip("chr")

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
    "0|0",
    "1|1",
    "2|2",
    "3|3",
    "4|4",
    "5|5",
    "6|6",
    "7|7",
    "0/0",
    "1/1",
    "2/2",
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

    out_stream = open(outputFilename, "w")
    out_stream.write(
        "chr\tchrN\tpos\tgeneID\ttranscriptID\ttranscript_pos\tSNP_id\tgenotype\n"
    )
    vcfline_processor = make_vcfline_processor(require_pass)

    reader = GffTranscriptReader()
    logging.info(
        f"..... Loading GTF file: {genecode_gz}")
    full_geneList = reader.loadGenes(genecode_gz)
    
    for chrId in chrRange(chr_start, chr_end, include_x_chromosome):
        
        chr = f"chr{chrId}"
        # normalized_chr = normalize_chromosome_name(chr)  # Strip 'chr'
        geneList = list(
        filter(lambda gene: gene.getSubstrate() == chr, full_geneList)
    )
        logging.info(
            f"..... Start loading gencode annotation chr {chrId} with {len(geneList)} genes")

        """
        construct a dict with unique exon region with transcript {exon_region:transcript list} 
        potential issue: 
        within the same gene, multiple transcript may cover the same exon region
        between genes, different transcripts may cover the same exon region
        """
        exon_region_to_transcripts = {}
        if DEBUG_GENES is not None:
            print(f"Debugging specific genes: {DEBUG_GENES}")
            geneList = list(filter(lambda x: x.getId() in DEBUG_GENES, geneList))

        for gene in geneList:
            # gene.getSubstrate() chr21
            # gene.getSubstrate().strip("chr")) 21
            if DEBUG_GENES is not None:
                logging.info(f">>>>>>>> DEBUGGING GENE {gene.getID()}")
                logging.info(f">>>>>>>>")
                logging.info(
                    f">>>>>>>> DEBUG1: Check every exon start-end region on each transcript"
                )
            if str(gene.getSubstrate().strip("chr")) == chrId:
                Num_transcript = gene.getNumTranscripts()
                chrom = gene.getSubstrate()  # column 1 --> chr21
                for n in range(Num_transcript):
                    transcript = gene.getIthTranscript(n)
                    if (
                        DEBUG_GENES is not None
                        and transcript.getId() == "ENST00000202917.5"
                    ):
                        print(f"{gene.getID()}-transcript {n} {transcript.getId()}")
                    # rawExons = transcript.UTR
                    # rawExons = transcript.exons
                    rawExons = transcript.getRawExons()
                    for exon in rawExons:

                        begin = exon.getBegin() + 1  # column 7
                        # print(begin)
                        end = exon.getEnd()  # column 8
                        exon_region = f"{chrId}:{begin}-{end}"
                        if not exon_region in exon_region_to_transcripts:
                            exon_region_to_transcripts[exon_region] = []
                        exon_region_to_transcripts[exon_region].append(transcript)
                        if DEBUG_GENES is not None:
                            print(f"{exon_region}")
                            # uncomment to debug duplicate transcript IDs as a possible optimization
                            # for existing in exon_region_to_transcripts[exon_region]:
                            # print(
                            #     f"      existing transcript @ {exon_region} {existing.getTranscriptId()} != {transcript.getTranscriptId()}"
                            # )

        """
        construct a dict with unique exon region with VCF records {exon_region:VCF records}
        """
        if DEBUG_GENES is not None:
            logging.info(f">>>>>>>>")
            logging.info(f">>>>>>>> DEBUG2: list all unique exon regions")
            logging.info(f">>>>>>>>")
            print(exon_region_to_transcripts.keys())
        # print(vcfFilename)
        exon_region_to_snp_infos = tabix_regions(
            exon_region_to_transcripts.keys(), vcfline_processor, vcfFilename
        )
        if DEBUG_GENES is not None:
            logging.info(f">>>>>>>>")
            logging.info(f">>>>>>>> DEBUG3: list all snp info corresponding to exon region")
            logging.info(f">>>>>>>>")
            print(exon_region_to_snp_infos)
        data = []
        variant_to_transcript_info = {}
        for exon_region in exon_region_to_snp_infos:
            snp_infos = exon_region_to_snp_infos[exon_region]
            transcripts = exon_region_to_transcripts[exon_region]
            # transcript list contains transcripts from the same gene or different genes
            for snp_info in snp_infos:  # loop through each variant record
                pos, rsid, genotype = snp_info
                if DEBUG_GENES is not None:
                    print(pos)
                for transcript in transcripts:
                    if chr == transcript.getSubstrate():
                        chr_pos = f"{chr}_{pos}"
                        if chr_pos not in variant_to_transcript_info:
                            variant_to_transcript_info[chr_pos] = []
                        variant_to_transcript_info[chr_pos].append(
                            (transcript, pos, rsid, genotype)
                        )
                    # else:
                    #     logging.info(f"Exception: chr {chr} - chr_pos {chr_pos} - rsid {rsid} - transcript {transcript.getSubstrate()} - transcriptID {transcript.getId()}")

        #print(">> dict variant_to_transcript_info")
        #print(f"len of dic: {len(variant_to_transcript_info)}")
        # print(">> print")
        for chr_pos in variant_to_transcript_info:
            gene_ids = set(
                [x[0].getGeneId() for x in variant_to_transcript_info[chr_pos]]
            )
            # print(f"{chr_pos} ---- {gene_ids}")
            if len(gene_ids) == 1:
                transcript, pos, rsid, genotype = variant_to_transcript_info[chr_pos][0]
                # chrom = transcript.getSubstrate()
                # chromN = chrom.strip("chr")
                transcriptCoord = transcript.mapToTranscript(int(pos))
                # print(
                #     f"......... write {transcript.getId()},{transcript.getGeneId()},{rsid}, {genotype}"
                # )
                data.append(
                    [
                        chrom,
                        chrId,
                        pos,
                        transcript.getGeneId(),
                        transcript.getId(),
                        transcriptCoord,
                        rsid,
                        genotype,
                    ]
                )
            # else:
            #     print("......... SKIPPING")
        data.sort(key=lambda r: (r[1], r[2]))
        for r in data:
            out_stream.write("\t".join(map(str, r)))
            out_stream.write("\n")
    out_stream.close()
    logging.info(f"..... Finish writing to {outputFilename}")
    logging.info(f"..... Done sample {sample}")

def count_all_het_sites_forpeaks(vcfFilename, outputFilename, annotation_file):
    out_stream = open(outputFilename, "w")
    out_stream.write(
        "chrN\tchr\tpos\tSNP_id\tgenotype\tpeak_start\tpeak_end\tpeakID\tgeneID\n"
    )

    """
    1. construct a dict with unique peak region with IDs {peak_region: chr,peak_start_pos,peak_end_pos,geneID, peakID}
    """
    peak_region_to_IDs = {}
    with gzip.open(annotation_file, "rt", newline="\n") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            key = f"{row['chrN']}:{row['peak_start']}-{row['peak_end']}"
            peak_region_to_IDs[key] = row

    """
    2. construct a dict with unique peak region with VCF records {peak_region:VCF records}
    """
    vcfline_processor = make_vcfline_processor(False)
    peak_region_to_snp_infos = tabix_regions(
        peak_region_to_IDs.keys(),
        vcfline_processor,
        vcfFilename,
    )
    """
    3. write output data
    """
    data = []
    for peak_region in peak_region_to_snp_infos:
        snp_infos = peak_region_to_snp_infos[peak_region]
        IDs = peak_region_to_IDs[peak_region]
        for snp_info in snp_infos:  # loop through each variant record
            pos, rsid, genotype = snp_info
            chrId = IDs["chrN"]
            chrom = "chr" + chrId
            peak_start = IDs["peak_start"]
            peak_end = IDs["peak_end"]
            peakID = IDs["peakID"]
            geneID = IDs["geneID"]
            data.append(
                [
                    chrId,
                    chrom,
                    pos,
                    rsid,
                    genotype,
                    peak_start,
                    peak_end,
                    peakID,
                    geneID,
                ]
            )
    data.sort(key=lambda r: (r[1], r[2]))
    for r in data:
        out_stream.write("\t".join(map(str, r)))
        out_stream.write("\n")
    out_stream.close()



#python extractHets.py /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21_hetSNP.tsv /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21.bihets.vcf.gz 21 21 /data2/reference/reference/gencode_chr

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