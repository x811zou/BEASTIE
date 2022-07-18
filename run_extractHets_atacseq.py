#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

import os
import sys
import argparse
import gzip
import csv
from BEASTIE.helpers import tabix_regions
from BEASTIE.extractHets import vcfline_processor_atacseq


def count_all_het_sites_forpeaks(
    vcfFilename,
    outputFilename,
    annotation_file,
):
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
    peak_region_to_snp_infos = tabix_regions(
        peak_region_to_IDs.keys(), vcfline_processor_atacseq, vcfFilename
    )
    """
    3. write output data
    """
    variant_to_peak_info = {}
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


# =========================================================================
# module load htslib/1.11
# python3 /hpc/group/allenlab/scarlett/script/BEASTIE/run_extractHets_atacseq.py /hpc/group/allenlab/scarlett/ATACseq/HG00096.no_chr.content.SNPs.hets.filtered.vcf.gz /hpc/group/allenlab/scarlett/ATACseq/peak_annotation.gz /hpc/group/allenlab/scarlett/ATACseq/hetSNP.tsv
# =========================================================================
# parameter arguments
parser = argparse.ArgumentParser()
parser.add_argument("vcfFilename", help="vcf file containing only heterozygous sites")
parser.add_argument("annotation_file", help="gzipped annotation files")
parser.add_argument("outputFilename", help="output file name")
# parser.add_argument("--chr_start", help="chr_start", default=1, type=int)
# parser.add_argument("--chr_end", help="chr end", default=22, type=int)
args = parser.parse_args()
vcfFilename = args.vcfFilename
annotation_file = args.annotation_file
outputFilename = args.outputFilename
#
print("derp")
count_all_het_sites_forpeaks(vcfFilename, outputFilename, annotation_file)
