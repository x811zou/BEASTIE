#=======================================================================================
# Run code as : python hetsMeta_for_mpileup_GSD.py 125249
# updated and saved in hardac on 07/05

##### make sure VCF file is zipped and has index: 
#module load vcftools/0.1.15-gcb01
#sample="383581"
#sample2="123667"
#outDir=/data/allenlab/scarlett/data/VCF/GSD/DNA_vcf
#vcftools --gzvcf /data/reddylab/GSD/kishnani_uwcmg_gsd_1.HF.final.vcf.gz --indv ${sample} --out $outDir/${sample2}.vcf --recode-INFO-all --recode
#bgzip -c ${sample2}.vcf.recode.vcf > ${sample2}.vcf.recode.vcf.gz
#tabix -p vcf ${sample2}.vcf.recode.vcf.gz


####### hetsMeta_GSD saves information for all hets by chr (Allchr_hets_all_transcript.tsv) for mpileup
####### hetsDict_GSD saves pickle file storing the gene:(pos set) for annotation (in altratio files)
#=======================================================================================

from __future__ import (absolute_import, division, print_function,
  unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
  chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import ProgramName
from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe
#import pandas as pd
#import numpy as np
import pickle
#from scipy import stats
from Gene import Gene
import sys
import time
#from HARDAC_hetsMeta_for_mpileup_GSD import isHeterozygous
import gzip
import os
""" Check if a genotype is heterozygous by testing whether them match with the 6 types of homozygous options
"""
def isHeterozygous(genotype):
    Homo = ["0|0","1|1","2|2","3|3","4|4","5|5","6|6","7|7","0/0","1/1","2/2"]
    if genotype in Homo:return False
    else:return True

def count_all_het_sites(sample,genom_dir,trans_dir,vcf_dir,file_dir,outDir):
    print(">>> sample %s"%(sample))
    vcfFilename = vcf_dir
    outputdir=outDir+"/hetsMeta/"
    ######
    #if str(sample) == "122687_2":
    #    vcfFilename="/data/reddylab/scarlett/1000G/data/VCF/GSD/"+str(sample)+".buffycoat.rnaseq.untrt.rep1.uniq.unmapped.variant_filtered.vcf.gz"
    #else:
    #    vcfFilename="/data/reddylab/scarlett/1000G/data/VCF/GSD/"+str(sample)+".buffycoat.rnaseq.untrt.rep1.unmapped.variant_filtered.vcf.gz"
    ######
    outputFilename=outputdir+"/Allchr_hets_all_transcript.tsv"
    out_stream = open(outputFilename, "w")
    out_stream.write("chr\tchrN\tgeneID\tgenomicCoord_pos\ttranscriptCoord\tSNP_id\tgenotype\n")
    #
    #output2Filename="/data/reddylab/scarlett/1000G/result/chrGeneHet_GSD/all_chr.tsv"
    #out_stream2 = open(output2Filename, "w")
    #out_stream2.write("Chr\tTotal_gene\tGene_w_hets\tTotal_biSNP\tBi_hets\n")
  # iterate chr1-22
    print("- all transcripts:")
    for Num in range(1,23):
        reader=GffTranscriptReader()
        print("we are working with chr" + str(Num))
        geneList=reader.loadGenes(file_dir+"chr"+str(Num))
        print(len(geneList),"Genes loaded")  #439  genes for chr22
        byGene={}
        byGene_trans={}
        total_biSNP = 0
        for gene in geneList:
            #Num_transcript = gene.getNumTranscripts()
            #for n in range(Num_transcript):
                #transcript=gene.getIthTranscript(n)
            #transID = transcript.getTranscriptId()# column 6
            transcript=gene.longestTranscript()
            geneID=transcript.getGeneId()      # column 5
            byGene[geneID]=byGene.get(geneID, set())
            byGene_trans[geneID]=byGene_trans.get(geneID, set())
            chrom=transcript.getSubstrate()      # column 1
            chromN=chrom.strip("chr")
                #rawExons=transcript.getRawExons()
                #for exon in rawExons:
                    # exon=rawExons[0]  ## for debugging only
            begin=gene.getBegin()    # column 7
            end=gene.getEnd()      # column 8
            cmd = "tabix " + vcfFilename + " "+chromN+":"+str(begin)+"-"+str(end)#         
                    #tabix /data/allenlab/scarlett/data/VCF/GSD/DNA_vcf/125249.vcf.recode.vcf.gz 1:10042358-10045556
            output=Pipe.run(cmd)
            if(len(output)<=9):continue
            lines=output.split("\n")
            for line in lines:
                fields=line.split("\t")
                #print(fields)
                if(fields[6]!="PASS"): continue
                pos=fields[1]                         # column 9
                #if(fields[2]=="."): continue
                rs = fields[2]                        # column 10
                genotype = fields[9].split(':')[0] # for HG00097
                if(not isHeterozygous(str(genotype))):continue # go back to the begining of the loop
                transcriptCoord=transcript.mapToTranscript(int(pos))#    
                total_biSNP += 1 ###########################
                chr_pos=chromN+"_"+pos
                byGene[geneID].add(chr_pos)
                trans_chr_pos = chromN+"_"+str(transcriptCoord)
                byGene_trans[geneID].add(trans_chr_pos)
                #if(not isHeterozygous(genotype)): continue # go back to the begining of the loop
                out_stream.write("\t".join([str(chrom),str(chromN),str(geneID),str(pos),str(transcriptCoord),str(rs),str(genotype),"\n"]))
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
    #    out_stream2.write("\t".join([str(Num),str(len(byGene)),str(non_empty_count),str(total_biSNP),str(len(all_set_list)),"\n"]))
    # loop finishes here
    #out_stream.close()
    #out_stream2.close()
    #  Done
        print("we finished with chr" + str(Num))
        with open(genom_dir+'chr'+str(Num)+'.pickle', 'wb') as handle:
            pickle.dump(byGene, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(trans_dir+'chr'+str(Num)+'.pickle', 'wb') as handle:
            pickle.dump(byGene_trans, handle, protocol=pickle.HIGHEST_PROTOCOL)
    out_stream.close()

if __name__ == "__main__":
    sample = sys.argv[1]
    vcfDir = sys.argv[2]
    outDir = sys.argv[3]
    # There is only one VCF file for GSD case
    genom_dir = outDir+'/hetsDict/genom_all/'
    trans_dir = outDir +'/hetsDict/trans_all/'
    file_dir = '/data/allenlab/scarlett/data/coding-noncoding/'
    count_all_het_sites(sample,genom_dir,trans_dir,vcfDir,file_dir,outDir)

