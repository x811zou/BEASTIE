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
            print(f"bad field line: '{line}'")
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


def count_all_het_sites(
    vcfFilename,
    outputFilename,
    chr_start,
    chr_end,
    gencode_path,
    require_pass,
    DEBUG_GENES=None,
    **kwargs,
):
    include_x_chromosome = kwargs.pop("include_x_chromosome", False)

    out_stream = open(outputFilename, "w")
    out_stream.write(
        "chr\tchrN\tpos\tgeneID\ttranscriptID\ttranscript_pos\tSNP_id\tgenotype\n"
    )

    vcfline_processor = make_vcfline_processor(require_pass)

    for chrId in chrRange(chr_start, chr_end, include_x_chromosome):
        reader = GffTranscriptReader()
        geneFile = gencode_path + f"/gencode.chr{chrId}.gtf.gz"
        geneList = reader.loadGenes(geneFile)
        chr = f"chr{chrId}"
        geneList = list(filter(lambda gene: gene.getSubstrate() == chr, geneList))
        logging.info(
            "..... start loading gencode annotation chr {0} with {1} genes".format(
                chrId, len(geneList)
            )
        )

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
                print(f">>>>>>>> DEBUGGING GENE {gene.getID()}")
                print(f">>>>>>>>")
                print(
                    f">>>>>>>> DEBUG1: Check every exon start-end region on each transcript"
                )
                print(f">>>>>>>>")
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
            print(f">>>>>>>>")
            print(f">>>>>>>> DEBUG2: list all unique exon regions")
            print(f">>>>>>>>")
            print(exon_region_to_transcripts.keys())
        # print(vcfFilename)
        exon_region_to_snp_infos = tabix_regions(
            exon_region_to_transcripts.keys(), vcfline_processor, vcfFilename
        )
        if DEBUG_GENES is not None:
            print(f">>>>>>>>")
            print(f">>>>>>>> DEBUG3: list all snp info corresponding to exon region")
            print(f">>>>>>>>")
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
                    # geneID = transcript.getGeneId()
                    assert chr == transcript.getSubstrate()
                    chr_pos = f"{chr}_{pos}"
                    if chr_pos not in variant_to_transcript_info:
                        variant_to_transcript_info[chr_pos] = []
                    variant_to_transcript_info[chr_pos].append(
                        (transcript, pos, rsid, genotype)
                    )

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
