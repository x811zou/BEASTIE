#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================

from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe
import time
import logging
import pandas as pd

""" Check if a genotype is heterozygous by testing whether them match with the 6 types of homozygous options
"""
def isHeterozygous(genotype):
    Homo = ["0|0","1|1","2|2","3|3","4|4","5|5","6|6","7|7","0/0","1/1","2/2"]
    if genotype in Homo:return False
    else:return True

def count_all_het_sites(sample,vcfFilename,file_dir,outputFilename,chr_start,chr_end):
    logging.info('..... We are looking at individual: {0}'.format(sample))
    out_stream = open(outputFilename, "w")
    out_stream.write("chr\tchrN\tgeneID\tpos\ttranscriptID\ttranscript_pos\tSNP_id\tgenotype\n")
    for Num in range(chr_start,chr_end+1):
        reader=GffTranscriptReader()
        geneList=reader.loadGenes(str(file_dir)+"/chr"+str(Num))
        logging.info('..... Working on chr {0} with {1} genes'.format(Num,len(geneList)))
        byGene={}
        byGene_trans={}
        total_biSNP = 0
        for gene in geneList:
            Num_transcript = gene.getNumTranscripts()
            for n in range(Num_transcript):
                transcript=gene.getIthTranscript(n)
                transID = transcript.getTranscriptId()# column 6
                geneID=transcript.getGeneId()      # column 5
                byGene[geneID]=byGene.get(geneID, set())
                byGene_trans[geneID]=byGene_trans.get(geneID, set())
                chrom=transcript.getSubstrate()      # column 1
                chromN=chrom.strip("chr")
                begin=gene.getBegin()    # column 7
                end=gene.getEnd()      # column 8
                cmd = "tabix " + vcfFilename + " "+chromN+":"+str(begin)+"-"+str(end)#         
                #tabix /data/allenlab/scarlett/data/VCF/GSD/DNA_vcf/125249.vcf.recode.vcf.gz 1:10042358-10045556
                output=Pipe.run(cmd)
                if(not output):
                    continue
                else:
                    if(len(output)<=9):continue
                    lines=output.split("\n")
                    for line in lines:
                        fields=line.split("\t")
                        if(fields[6]!="PASS"): continue
                        pos=fields[1]                         # column 9
                        rs = fields[2]                        # column 10
                        genotype = fields[9].split(':')[0] # for HG00097
                        if(not isHeterozygous(str(genotype))):continue # go back to the begining of the loop
                        transcriptCoord=transcript.mapToTranscript(int(pos))#    
                        total_biSNP += 1 
                        chr_pos=chromN+"_"+pos
                        byGene[geneID].add(chr_pos)
                        out_stream.write("\t".join([str(chrom),str(chromN),str(geneID),str(pos),str(transID),str(transcriptCoord),str(rs),str(genotype),"\n"]))
        time.sleep(1/len(geneList))
        # write up the basic information 
        all_set_list = []
        non_empty_count = 0
        non_empty_count_len = 0
        for each_key in byGene:
            all_set_list.extend(list(byGene[each_key]))
            if len(byGene[each_key])!= 0:
                non_empty_count += 1
                non_empty_count_len += len(byGene[each_key])
    out_stream.close()



