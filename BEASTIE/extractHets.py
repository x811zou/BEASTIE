#!/usr/bin/env python
# =========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
# =========================================================================

import logging
import os
import sys
import pandas as pd
from pkg_resources import resource_filename
import gzip
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


def count_all_het_sites(
    sample,
    vcfFilename,
    outputFilename,
    chr_start,
    chr_end,
    gencode_path,
    DEBUG_GENES=None,
):

    logging.info("..... We are looking at individual: {0}".format(sample))
    filename = os.path.splitext(str(outputFilename))[0]
    count = 0

    for Num in range(chr_start, chr_end + 1):
        reader = GffTranscriptReader()
        geneFile = gencode_path + f"/gencode.chr{Num}.gtf.gz"
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

            byGene = {}
            total_biSNP = 0
            region_str_to_transcripts = {}
            if DEBUG_GENES is not None:
                print(f"Debugging specific genes: {DEBUG_GENES}")
                geneList = list(filter(lambda x: x.getId() in DEBUG_GENES, geneList))
            # else:
            #     print("Not debugging mode")
            for gene in geneList:
                # print(gene.getSubstrate().strip("chr"))
                if str(gene.getSubstrate().strip("chr")) == str(Num):
                    Num_transcript = gene.getNumTranscripts()
                    for n in range(Num_transcript):
                        transcript = gene.getIthTranscript(n)
                        geneID = transcript.getGeneId()  # column 5
                        byGene[geneID] = byGene.get(geneID, set())
                        chrom = transcript.getSubstrate()  # column 1
                        chromN = chrom.strip("chr")
                        # rawExons = transcript.UTR
                        rawExons = transcript.getRawExons()
                        # rawExons = transcript.exons
                        for exon in rawExons:
                            begin = exon.getBegin()  # column 7
                            end = exon.getEnd()  # column 8
                            region_str = f"{chromN}:{begin}-{end}"
                            print(region_str)
                            #
                            print(
                                f"{geneID}-transcript {n} {transcript.getId()} {chromN}:{begin}-{end}"
                            )
                            if not region_str in region_str_to_transcripts:
                                region_str_to_transcripts[region_str] = []
                            # uncomment to debug duplicate transcript IDs as a possible optimization
                            # else:
                            #     for existing in region_str_to_transcripts[region_str]:
                            #         print(f"existing transcript @ {region_str} {existing.getTranscriptId()} != {transcript.getTranscriptId()}")
                            region_str_to_transcripts[region_str].append(transcript)
            # print(region_str_to_transcripts)
            CHUNK_SIZE = 1000
            outputs = []
            for x in chunk_iter(iter(region_str_to_transcripts.keys()), CHUNK_SIZE):
                regions = " ".join(x)
                cmd = f"tabix --separate-regions {vcfFilename} {regions}"
                output = Pipe.run(cmd)
                if len(output) > 0:
                    outputs.append(output)

            data = []
            for output in outputs:
                lines = output.split("\n")
                transcripts = None
                for line in lines:
                    if line.startswith("#"):
                        region_str = line[1:]
                        assert region_str in region_str_to_transcripts
                        transcripts = region_str_to_transcripts[region_str]
                        continue

                    assert transcripts is not None

                    fields = line.split("\t")
                    assert len(fields) >= 10

                    if fields[6] != "PASS":
                        continue

                    genotype = fields[9].split(":")[0]

                    if not isHeterozygous(genotype):
                        continue

                    pos = int(fields[1])
                    # print(f"pos: {pos} region_str: {region_str}")
                    rs = fields[2]

                    for transcript in transcripts:
                        transcriptCoord = transcript.mapToTranscript(pos)
                        chrom = transcript.getSubstrate()
                        chromN = int(chrom.strip("chr"))
                        total_biSNP += 1
                        chr_pos = f"{chromN}_{pos}"
                        # currently we are not includeing the variants that are already covered by genes processed before to ensure that not two genes share same variants
                        if not chr_pos in byGene[geneID]:
                            print(f"{transcript.getGeneId()} - {chr_pos}")
                            byGene[geneID].add(chr_pos)
                            data.append(
                                [
                                    chrom,
                                    chromN,
                                    pos,
                                    transcript.getGeneId(),
                                    transcript.getTranscriptId(),
                                    transcriptCoord,
                                    rs,
                                    genotype,
                                ]
                            )
                            # break

            data.sort(key=lambda r: (r[1], r[3]))
            # print(data)

            out_stream = open(outputFile, "w")
            out_stream.write(
                "chr\tchrN\tpos\tgeneID\ttranscriptID\ttranscript_pos\tSNP_id\tgenotype\n"
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
