#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
import logging
import os

import pandas as pd

from pkg_resources import resource_filename

def annotateAF(ancestry, hetSNP, out_AF, ref_dir):
    AF_file = os.path.join(ref_dir, "AF_1_22.tsv")
    if not os.path.isfile(out_AF):
        logging.info('..... start reading 1000 Genome AF annotation file')
        AF=pd.read_csv(AF_file, header=0,sep='\t')
        logging.info('..... finish reading 1000 Genome AF annotation file')
        data=pd.read_csv(hetSNP,sep="\t",header=0,index_col=False)
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
        data_AF.to_csv(out_AF,sep='\t')
        logging.info('..... finish annotating AF for SNPs, file save at {0}'.format(out_AF))
    else:
        logging.info('..... skip annotating AF for SNPs, file already saved at {0}'.format(out_AF))

def annotateLD(prefix,ancestry,hetSNP_intersect_unique,out,LD_token,chr_start,chr_end,meta):
    annotate_ld_new = resource_filename('BEASTIE', 'annotate_LD_new.R')

    if not os.path.isfile(meta):
        cmd = f"Rscript --vanilla {annotate_ld_new} {prefix} {ancestry} {hetSNP_intersect_unique} {out} {LD_token} {chr_start} {chr_end} {meta}"
        os.system(cmd)
        logging.info(f"..... finish annotating LD for SNP pairs, file save at {meta}")
    else:
        logging.info(f"..... skip annotating LD for SNP pairs, file already saved at {meta}")
