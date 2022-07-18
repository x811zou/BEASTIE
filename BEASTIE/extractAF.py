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
from BEASTIE.helpers import tabix_regions
from .misc_tools.GffTranscriptReader import GffTranscriptReader
from .misc_tools.Pipe import Pipe
from glob import glob

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

def vcfline_processor(line):
    fields = line.split("\t")
    if len(fields) < 10:
        print(f"bad field line: '{line}'")
    assert len(fields) >= 10

    if fields[6] != "PASS":
        return None

    AF = fields[7]

    pos = int(fields[1])
    rs = fields[2]

    return (pos, rs, AF)


def count_all_sites(
    vcfFile_path,
    outputFilename,
    chr_start,
    chr_end,
    gencode_path,
    DEBUG_GENES=None
):
    filename = os.path.splitext(str(outputFilename))[0]
    count = 0
    if chr_end == "X":
        chromosome_list=list(range(1,23,1))
        chromosome_list.append("X")
        chromosome_list.append("Y")
    if chr_start == "X":
        chromosome_list=[]
        chromosome_list.append("X")
        chromosome_list.append("Y")
    for Num in chromosome_list:
        vcfFilename = glob(vcfFile_path + f"/ALL.chr{Num}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes*vcf.gz")[0]
        reader = GffTranscriptReader()
        geneFile = gencode_path + f"gencode.chr{Num}.gtf.gz" 
        geneList = reader.loadGenes(geneFile)
        outputFile = filename + ".chr" + str(Num) + ".TEMP.tsv"
        if not os.path.isfile(outputFile):
            chr = f"chr{Num}"
            geneList = list(filter(lambda gene: gene.getSubstrate() == chr, geneList))
            logging.info(
                "..... start loading gencode annotation chr {0} with {1} genes".format(
                    Num, len(geneList)
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
                print(f">>>> {gene.getID()}")
                if str(gene.getSubstrate().strip("chr")) == str(Num):
                    Num_transcript = gene.getNumTranscripts()
                    chrom = gene.getSubstrate()  # column 1 --> chr21
                    chromN = Num
                    for n in range(Num_transcript):
                        transcript = gene.getIthTranscript(n)
                        # rawExons = transcript.UTR
                        # rawExons = transcript.exons
                        rawExons = transcript.getRawExons()
                        for exon in rawExons:
                            begin = exon.getBegin()  # column 7
                            end = exon.getEnd()  # column 8
                            exon_region = f"{chromN}:{begin}-{end}"
                            if not exon_region in exon_region_to_transcripts:
                                exon_region_to_transcripts[exon_region] = []
                            exon_region_to_transcripts[exon_region].append(transcript)
                            # print(
                            #     f"{gene.getID()}-transcript {n} {transcript.getId()} {chromN}:{begin}-{end}"
                            # )
                            # uncomment to debug duplicate transcript IDs as a possible optimization
                            #     for existing in region_str_to_transcripts[region_str]:
                            #         print(f"existing transcript @ {region_str} {existing.getTranscriptId()} != {transcript.getTranscriptId()}")

            """
            construct a dict with unique exon region with VCF records {exon_region:VCF records}
            """
            #print(exon_region_to_transcripts.keys())
            #(['1:249200570-249200665', '1:249211477-249211600'])
            #print(vcfFilename)
            #['/datacommons/allenlab/scarlett/data/VCF/1000_genome/20130502/bgzip/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz']
            
            exon_region_to_snp_infos = tabix_regions(
                exon_region_to_transcripts.keys(), vcfline_processor, vcfFilename
            )

            data = []
            variant_to_transcript_info = {}
            for exon_region in exon_region_to_snp_infos:
                snp_infos = exon_region_to_snp_infos[exon_region]
                transcripts = exon_region_to_transcripts[exon_region]
                # transcript list contains transcripts from the same gene or different genes
                for snp_info in snp_infos:  # loop through each variant record
                    pos, rsid, AF = snp_info
                    #AC=25;AF=0.00499201;AN=5008;NS=2504;DP=12821;AMR_AF=0.0043;AFR_AF=0.0038;EUR_AF=0.0099;SAS_AF=0.0061;EAS_AF=0.001;AA=.|||;VT=SNP
                    for transcript in transcripts:
                        # geneID = transcript.getGeneId()
                        assert chr == transcript.getSubstrate()
                        chr_pos = f"{chr}_{pos}"
                        if chr_pos not in variant_to_transcript_info:
                            variant_to_transcript_info[chr_pos] = []
                        variant_to_transcript_info[chr_pos].append(
                            (transcript, pos, rsid, AF)
                        )
            print(">> dict variant_to_transcript_info")
            print(f"len of dic: {len(variant_to_transcript_info)}")
            print(">> print")
            for chr_pos in variant_to_transcript_info:
                gene_ids = set(
                    [x[0].getGeneId() for x in variant_to_transcript_info[chr_pos]]
                )
                print(f"{chr_pos} ---- {gene_ids}")
                if len(gene_ids) == 1:
                    pos, rsid, AF = variant_to_transcript_info[
                        chr_pos
                    ][0]
                    # chrom = transcript.getSubstrate()
                    # chromN = chrom.strip("chr")
                    data.append(
                        [
                            chrom,
                            rsid,
                            pos,
                            AF,
                        ]
                    )
                # else:
                #     print("......... SKIPPING")
            data.sort(key=lambda r: (r[1], r[2]))
            out_stream = open(outputFile, "w")
            #EAS_AF,AMR_AF,AFR_AF,EUR_AF,SAS_AF
            out_stream.write(
                "chr,rsid,pos,AF\n"
            )
            for r in data:
                out_stream.write("\t".join(map(str, r)))
                out_stream.write("\n")
            out_stream.close()
        else:
            logging.info("..... chr{0} existed".format(Num))
        count += 1
        if count == 1:
            data = pd.read_csv(outputFile, sep="\t", header=0, index_col=False)
            data0 = data
        else:
            data = pd.read_csv(outputFile, sep="\t", header=0, index_col=False)
            data0 = pd.concat([data0, data])
    # data0 = data0[
    #     ["chr", "chrN", "pos", "transcriptID"，"transcript_pos", "geneID", "SNP_id", "genotype"]
    # ]
    data0.drop_duplicates()
    data0.to_csv(outputFilename, sep="\t", header=True, index=False)
    for files in os.listdir(os.path.dirname(outputFilename)):
        if "TEMP" in files:
            logging.info("..... remove created TEMP files: {0}".format(files))
            os.remove(os.path.dirname(outputFilename) + "/" + files)
