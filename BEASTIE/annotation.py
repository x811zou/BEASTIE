#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
import sys
import os
import calendar
import time
import pandas as pd
import logging

def annotateAF(ancestry,AF_file,hetSNP_intersect_unique,hetSNP_AF_file):
    logging.info('..... start reading 1000 Genome AF annotation file')
    AF=pd.read_csv(AF_file, header=0,sep='\t') 
    logging.info('..... finish reading 1000 Genome AF annotation file')
    data=pd.read_csv(hetSNP_intersect_unique,sep="\t",header=0,index_col=False)
    data=data.drop(['if_Indel','if_SV','if_SNP','if_biallelic','altRatio','variantID'],axis=1)
    if ancestry == "EUR":
        AF=AF[['chr','pos','rsid','EUR_AF']]
        AF=AF.rename(columns={'EUR_AF': 'AF'})
    elif ancestry == "AFR":
        AF=AF[['chr','pos','rsid','AFR_AF']]
        AF=AF.rename(columns={'AFR_AF': 'AF'})
    elif ancestry == "EAS":
        AF=AF[['chr','pos','rsid','EAS_AF']]  
        AF=AF.rename(columns={'EAS_AF': 'AF'})
    elif ancestry == "AMR":
        AF=AF[['chr','pos','rsid','AMR_AF']]
        AF=AF.rename(columns={'AMR_AF': 'AF'})  
    elif ancestry == "SAS":
        AF=AF[['chr','pos','rsid','SAS_AF']]  
        AF=AF.rename(columns={'SAS_AF': 'AF'})
    data_AF=pd.merge(data,AF,on=['chr','pos'],how='left')
    data_AF=data_AF.drop_duplicates()
    data_AF.to_csv(hetSNP_AF_file,sep='\t')


def run(prefix,tmp_dir,ancestry,ref_dir,LD_token,chr_start,chr_end,meta_file,hetSNP_intersect_unique):
    AF_file = os.path.join(ref_dir,"AF_1_22.tsv")
    out_AF = '{0}_hetSNP_intersect_unique_AF.tsv'.format(os.path.join(tmp_dir,prefix))
    if not os.path.isfile(out_AF):
        annotateAF(ancestry,AF_file,hetSNP_intersect_unique,out_AF)
        logging.info('..... finish annotating AF for SNPs, file save at {0}'.format(out_AF))
    else:
        logging.info('..... skip annotating AF for SNPs, file already saved at {0}'.format(out_AF)) 
    if not os.path.isfile(meta_file):
        print(prefix)
        print(out_AF)
        tmp_dir=tmp_dir+"/"
        print(tmp_dir)
        cmd="Rscript --vanilla annotate_LD_new.R %s %s %s %s %s %d %d"%(prefix,ancestry,out_AF,tmp_dir,LD_token,int(chr_start),int(chr_end))
        os.system(cmd)
        logging.info('..... finish annotating LD for SNP pairs, file save at {0}'.format(meta_file))
    else:
        logging.info('..... skip annotating LD for SNP pairs, file already saved at {0}'.format(meta_file))
