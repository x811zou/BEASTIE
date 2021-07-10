#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
import os
import logging
import sys
import pandas as pd
from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe
from extractHets import count_all_het_sites
from parse_mpileup import Parse_mpileup_allChr
from intersect_hets import Intersect_exonicHetSnps
import annotation

def check_file_existence(prefix,in_path,out,model,vcf,ref_dir,pileup,hetSNP,hetSNP_intersect_unique,parsed_pileup,meta,alpha):
    out,tmp = create_output_directory(in_path,out)
    split=os.path.split(model)
    STAN=split[0]
    modelName=split[1]
    ##### STAN model
    stan_model = False
    if os.path.exists(STAN):
        for filename in os.listdir(STAN):
            if '.stan' in filename:
                stan_model = True
    else:
        logging.error('Oops! STAN path {0} doesn\'t exist. Please try again ...'.format(STAN))
        sys.exit(1)
    if stan_model is False:
        logging.error('Oops! STAN model {0} doesn\'t exist in {1}. Please try again ...'.format(modelName,STAN))
        sys.exit(1)
    ##### vcf
    if not os.path.isfile(vcf):
        logging.error('Oops! vcf file {0} doesn\'t exist. Please try again ...'.format(vcf))
        exit(1)
    ##### vcfgz
    vcfgz = '{0}.gz'.format(vcf)
    if not os.path.isfile(vcfgz):
        logging.warning('Oops! VCFgz file {0} not found. We will generate that for you ...'.format(vcfgz))
        cmd="bgzip -c %s > %s"%(vcf,vcfgz)
        os.system(cmd)
        cmd="tabix -p vcf %s"%(vcfgz)
        os.system(cmd)
    ##### reference directory : AF file and gencode directory
    AF_file = False
    if os.path.exists(ref_dir):
        for filename in os.listdir(ref_dir):
            if 'AF_1_22.tsv' in filename:
                AF_file = True
                AF_file_name=filename
            if 'gencode' in filename:
                gencode_dir = True
                gencode_dir_name=os.path.join(ref_dir,"gencode_chr")
    else:
        logging.error('Oops! REFERENCE path {0} doesn\'t exist. Please try again ...'.format(ref_dir))
        sys.exit(1)
    if AF_file is False:
        logging.error('Oops! AF file {0} doesn\'t exist in {1}. Please try again ...'.format(AF_file_name,ref_dir))
        sys.exit(1)
    if gencode_dir is False:
        logging.error('Oops! gencode directory {0} doesn\'t exist. Please try again ...'.format(gencode_dir_name))
        sys.exit(1)
    ##### pileup_file
    pileup_file = os.path.join(in_path,pileup)
    if not os.path.isfile(pileup_file):
        logging.error('Oops!  pileup file {0} doesn\'t exist in {1}. Please try again ...'.format(pileup,in_path))
    ##### hetSNP_file
    if hetSNP is not None:
        hetSNP_file = in_path+hetSNP
        if not os.path.isfile(hetSNP_file):
            logging.warning('Alright, hetSNP file {0} doesn\'t exist in {1}. We will generate that for you ...'.format(hetSNP,in_path))
        else:
            logging.warning('Found existed hetSNP file {0} doesn\'t exist in {1}.'.format(hetSNP,in_path))
    else:
        hetSNP_file = '{0}_hetSNP.tsv'.format(os.path.join(tmp,prefix))
        logging.info('We will generate {0} for you ...'.format(hetSNP_file))
    ##### hetSNP_filec
    if hetSNP_intersect_unique is not None:
        hetSNP_intersect_unique_file = in_path+hetSNP_intersect_unique
        if not os.path.isfile(hetSNP_intersect_unique_file):
            logging.warning('Alright, hetSNP file {0} doesn\'t exist in {1}. We will generate that for you ...'.format(hetSNP,in_path))
        else:
            logging.warning('Found existed hetSNP file {0} doesn\'t exist in {1}.'.format(hetSNP_intersect_unique,in_path))
    else:
        hetSNP_intersect_unique_file = '{0}_hetSNP_intersect_unique.tsv'.format(os.path.join(tmp,prefix))
        logging.info('We will generate {0} for you ...'.format(hetSNP_intersect_unique_file))
        hetSNP_intersect_unique_forlambda_file = '{0}_hetSNP_intersect_unique_forLambda.tsv'.format(os.path.join(tmp,prefix))
        hetSNP_intersect_unique_lambdaPredicted_file= '{0}_hetSNP_intersect_unique_alpha'.format(os.path.join(tmp,prefix))+str(alpha)+'_lambdaPredicted.tsv'
    ##### meta_file
    if meta is not None:
        meta_file = os.path.join(in_path,meta)
        if not os.path.isfile(meta_file):
            logging.warning('Alright, meta file {0} doesn\'t exist in {1}. We will generate that for you ...'.format(meta,in_path))
    else:
            meta_file = '{0}_logisticReg_input.tsv'.format(os.path.join(tmp,prefix))
            logging.info('We will generate {0} for you ...'.format(meta_file))
    ##### parsed_pileup_file
    if parsed_pileup is not None:
        parsed_pileup_file = os.path.join(in_path,parsed_pileup)
        if not os.path.isfile(parsed_pileup_file):
            logging.warning('Alright, parsed_pileup file {0} doesn\'t exist in {1}. We will generate that for you ...'.format(parsed_pileup,in_path))
    else:
            parsed_pileup_file = '{0}_parsed_pileup.tsv'.format(os.path.join(tmp,prefix))
            logging.info('We will generate {0} for you ...'.format(parsed_pileup_file))
    return out,tmp,vcfgz,pileup_file,hetSNP_file,meta_file,hetSNP_intersect_unique_file,hetSNP_intersect_unique_forlambda_file,hetSNP_intersect_unique_lambdaPredicted_file,parsed_pileup_file,gencode_dir_name


def create_output_directory(in_path,out):
   out = os.path.join(in_path,out)
   if not os.path.isdir(out):
      logging.info('Creating output directory : {0}'.format(out))
      os.mkdir(out)
   else:
      logging.info('Found existed output directory : {0}'.format(out))
   tmp = os.path.join(out,"TEMP")
   if not os.path.isdir(tmp):
      logging.info('Creating temporary output directory : {0}'.format(tmp))
      os.mkdir(tmp)
   else:
      logging.info('Found existed temporary output directory : {0}'.format(tmp))
   return out,tmp


def run(prefix,vcf_sample_name,in_path,out,model,vcf,ref_dir,ancestry,chr_start,chr_end,min_total_cov,min_single_cov,read_length,LD_token,pileup,hetSNP,hetSNP_intersect_unique,parsed_pileup,meta,alpha):
    out,tmp,vcfgz,pileup,hetSNP,meta,hetSNP_intersect_unique,hetSNP_intersect_unique_forlambda_file,hetSNP_intersect_unique_lambdaPredicted_file,parsed_pileup,gencode_dir = check_file_existence(prefix,in_path,out,model,vcf,ref_dir,pileup,hetSNP,hetSNP_intersect_unique,parsed_pileup,meta,alpha)
    ##### Generate hetSNP file: extract heterozygous bi-allelic SNPs for specific chromosomes from all gencode transcripts
    if not os.path.exists(hetSNP):
        logging.info('>>>>>>>>>>')
        logging.info('>>>>>>>>>> Starting step 1.1 : extract heterozygous bi-allelic SNPs')
        count_all_het_sites(prefix,vcfgz,gencode_dir,hetSNP,int(chr_start),int(chr_end))
        logging.info('..... file save at {0}'.format(hetSNP))
    else:
        logging.info('>>>>>>>>>>')
        logging.info('>>>>>>>>>> Skipping step 1.1 : extract heterozygous bi-allelic SNPs')
        data=pd.read_csv(hetSNP,sep="\t",header=0,index_col=False)
        if data.shape[1]<2:
            #print(data.head())
            os.remove(hetSNP)
            logging.error('Existed hetSNP file is empty, please try again!')
            sys.exit(1)
        else:
            logging.info('..... hetSNP file exists at {0}'.format(hetSNP))            

    ##### Generate parsed pileup file: parse samtools mpile up output files
    if not os.path.exists(parsed_pileup):
        logging.info('>>>>>>>>>>')
        logging.info('>>>>>>>>>> Starting step 1.2 : parse samtools pileup result')
        Parse_mpileup_allChr(vcf_sample_name,vcfgz,pileup,min_total_cov,min_single_cov,parsed_pileup)
        logging.info('..... file save at {0}'.format(parsed_pileup))
    else:
        logging.info('>>>>>>>>>>')
        logging.info('>>>>>>>>>> Skipping step 1.2 : parse samtools pileup result')
        logging.info('..... parsed pileup file exists at {0}'.format(parsed_pileup))

    ##### Thinning reads: one reads only count once
    if (not os.path.exists(hetSNP_intersect_unique)) or (not os.path.exists(hetSNP_intersect_unique_forlambda_file)):
        logging.info('>>>>>>>>>>')
        logging.info('>>>>>>>>>> Starting step 1.3 : keeping unique reads')
        Intersect_exonicHetSnps(parsed_pileup,hetSNP,read_length,min_total_cov,min_single_cov,hetSNP_intersect_unique,hetSNP_intersect_unique_forlambda_file)
        logging.info('..... hetSNP file with filtered sites save at {0}'.format(hetSNP_intersect_unique))
        logging.info('..... hetSNP file with filtered sites prepared for lambda model save at {0}'.format(hetSNP_intersect_unique_forlambda_file))
    else:
        logging.info('>>>>>>>>>>')
        logging.info('>>>>>>>>>> Skipping step 1.3 : keeping unique reads')
        logging.info('..... hetSNP file with filtered sites exists at {0}'.format(hetSNP_intersect_unique))
        logging.info('..... hetSNP file with filtered sites prepared for lambda model exists at {0}'.format(hetSNP_intersect_unique_forlambda_file))

   ##### Annotation: AF and LD
    if not os.path.isfile(meta):
        logging.info('>>>>>>>>>>')
        logging.info('>>>>>>>>>> Starting step 1.4 : annotate AF and LD information')
        annotation.run(prefix,tmp,ancestry,ref_dir,LD_token,chr_start,chr_end,meta,hetSNP_intersect_unique)
        logging.info('>>>>>>>>>> Finishing step 1.4 : annotated file save at {0}'.format(meta))
    else:
        logging.info('>>>>>>>>>>')
        logging.info('>>>>>>>>>> Skipping step 1.4 : annotated file exists at {0}'.format(meta)) 

    return hetSNP_intersect_unique,meta,hetSNP_intersect_unique_forlambda_file,hetSNP_intersect_unique_lambdaPredicted_file
