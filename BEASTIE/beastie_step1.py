#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
import os
import logging
import sys
import pandas as pd
from misc_tools.GffTranscriptReader import GffTranscriptReader
from misc_tools.Pipe import Pipe
from extractHets import count_all_het_sites
from parse_mpileup import Parse_mpileup_allChr
from intersect_hets import Intersect_exonicHetSnps
import annotation
import shutil

def check_file_existence(prefix,in_path,out,model,vcf,ref_dir,pileup,hetSNP,parsed_pileup,sigma,alpha,WARMUP,KEEPER,min_single_cov,min_total_cov,chr_start,chr_end):
    out,common = create_output_directory(in_path,out,sigma,alpha,WARMUP,KEEPER,min_single_cov,min_total_cov)
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
    ##### vcf & vcfgz
    vcfgz = '{0}.gz'.format(vcf)
    vcfgztbi = '{0}.tbi'.format(vcfgz)
    if (not os.path.isfile(vcf)) :
        logging.error('Oops! vcf file {0} doesn\'t exist. Please try again ...'.format(vcf,vcfgz))
        exit(1)
    elif (os.path.isfile(vcf) and (not os.path.isfile(vcfgz)) and (not os.path.isfile(vcfgztbi))):
        logging.warning('We will generate vcfgz for you ...'.format(vcfgz))
        cmd="bgzip -c %s > %s"%(vcf,vcfgz)
        os.system(cmd)
        cmd="tabix -p vcf %s"%(vcfgz)
        os.system(cmd)
    elif (os.path.isfile(vcfgz)) and (not os.path.isfile(vcf)):
        logging.warning('Oops! VCF file {0} not found. We will generate that for you ...'.format(vcf))
        cmd="gunzip %s > %s"%(vcfgz,vcf)
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
        hetSNP_file = '{0}_hetSNP_chr{1}-{2}.tsv'.format(os.path.join(common,prefix),chr_start,chr_end)
        logging.info('We will generate {0} for you ...'.format(hetSNP_file))
    ##### TEMP output generation: hetSNP_file
    hetSNP_intersect_unique_file = '{0}_hetSNP_intersected_filtered.TEMP.tsv'.format(os.path.join(out,prefix))
    hetSNP_intersect_unique_forlambda_file = '{0}_hetSNP_intersected_filtered_forLambda.TEMP.tsv'.format(os.path.join(out,prefix))
    hetSNP_intersect_unique_lambdaPredicted_file= '{0}_hetSNP_intersected_filtered_alphaPredicted.TEMP.tsv'.format(os.path.join(out,prefix))
    logging.info('We will generate intermediate file {0} ...'.format(hetSNP_intersect_unique_file))
    logging.info('We will generate intermediate file {0} ...'.format(hetSNP_intersect_unique_forlambda_file))
    logging.info('We will generate intermediate file {0} ...'.format(hetSNP_intersect_unique_lambdaPredicted_file))
    ##### TEMP output generation: meta_file
    meta_file = '{0}_meta.TEMP.tsv'.format(os.path.join(out,prefix))
    logging.info('We will generate intermedate file {0} for you ...'.format(meta_file))
    ##### parsed_pileup_file
    if parsed_pileup is not None:
        parsed_pileup_file = os.path.join(in_path,parsed_pileup)
        if not os.path.isfile(parsed_pileup_file):
            logging.warning('Alright, parsed_pileup file {0} doesn\'t exist in {1}. We will generate that for you ...'.format(parsed_pileup,in_path))
    else:
        parsed_pileup_file = '{0}_parsed_pileup_chr{1}-{2}.tsv'.format(os.path.join(common,prefix),chr_start,chr_end)
        logging.info('We will generate {0} for you ...'.format(parsed_pileup_file))
    return out,common,vcfgz,pileup_file,hetSNP_file,meta_file,hetSNP_intersect_unique_file,hetSNP_intersect_unique_forlambda_file,hetSNP_intersect_unique_lambdaPredicted_file,parsed_pileup_file,gencode_dir_name

def create_output_directory(in_path,out,sigma,alpha,WARMUP,KEEPER,min_single_cov,min_total_cov):
   common = os.path.join(in_path,out)
   if not os.path.isdir(common):
       os.mkdir(common)
       logging.info('Creating output directory : {0}'.format(common))
   else:
       logging.info('Found existed output directory : {0}'.format(common))
   out=common+"/s-"+str(sigma)+"_a-"+str(alpha)+"_sinCov"+str(min_single_cov)+"_totCov"+str(min_total_cov)+"_W"+str(WARMUP)+"K"+str(KEEPER)
   if not os.path.isdir(out):
      logging.info('Creating specific output directory : {0}'.format(out))
      os.mkdir(out)
   else:
      #shutil.rmtree(out)
      #os.mkdir(out)
      logging.info('Found existed specific output directory : {0}, remove content ...'.format(out))
   return out,common


def run(sigma,alpha,WARMUP,KEEPER,prefix,vcf_sample_name,in_path,out,model,vcf,ref_dir,ancestry,chr_start,chr_end,min_total_cov,min_single_cov,read_length,LD_token,pileup,hetSNP,parsed_pileup):
    out,common,vcfgz,pileup,hetSNP,meta,hetSNP_intersect_unique,hetSNP_intersect_unique_forlambda_file,hetSNP_intersect_unique_lambdaPredicted_file,parsed_pileup,gencode_dir = check_file_existence(prefix,in_path,out,model,vcf,ref_dir,pileup,hetSNP,parsed_pileup,sigma,alpha,WARMUP,KEEPER,min_single_cov,min_total_cov,chr_start,chr_end)
    ##### 1.1 Generate hetSNP file: extract heterozygous bi-allelic SNPs for specific chromosomes from all gencode transcripts
    if not os.path.exists(hetSNP):
        logging.info('=================')
        logging.info('================= Starting common step 1.1')
        logging.info('..... start extracting heterozygous bi-allelic SNPs')
        count_all_het_sites(common,prefix,vcfgz,gencode_dir,hetSNP,int(chr_start),int(chr_end))
    else:
        logging.info('=================')
        logging.info('================= Skipping common step 1.1')
    data11=pd.read_csv(hetSNP,sep="\t",header=0,index_col=False)
    if data11.shape[0]<2:
        os.remove(hetSNP)
        logging.error('..... existed hetSNP file is empty, please try again!')
        sys.exit(1)
    else:
        logging.info('..... generated heterozygous bi-allelic SNPs data extracted from VCF can be found at {0}'.format(hetSNP))

    ##### 1.2 Generate parsed pileup file: parse samtools mpile up output files
    if not os.path.exists(parsed_pileup):
        logging.info('=================')
        logging.info('================= Starting common step 1.2')
        logging.info('..... start parsing samtools pileup result')
        Parse_mpileup_allChr(vcf_sample_name,vcfgz,pileup,min_total_cov,min_single_cov,parsed_pileup)
    else:
        logging.info('=================')
        logging.info('================= Skipping common step 1.2')
    data12=pd.read_csv(parsed_pileup,sep="\t",header=0,index_col=False)
    if data12.shape[1]<2:
        os.remove(parsed_pileup)
        logging.error('..... existed parsed pileup file is empty, please try again!')
        sys.exit(1)
    else:
        logging.info('..... generated parsed pileup file can be found at {0}'.format(parsed_pileup))

   ##### 1.3 Annotation: AF
    hetSNP_AF = '{0}_AF.tsv'.format(os.path.splitext(hetSNP)[0])
    print(hetSNP_AF)
    if not os.path.isfile(hetSNP_AF):
        logging.info('=================')
        logging.info('================= Starting common step 1.3')
        logging.info('..... start annotating AF and LD information')
        annotation.annotateAF(ancestry,hetSNP,hetSNP_AF,ref_dir)
    else:
        logging.info('=================')
        logging.info('================= Skipping common step 1.3')
        data13=pd.read_csv(hetSNP_AF,sep="\t",header=0,index_col=False)
    if data13.shape[0]<2:
        os.remove(hetSNP_AF)
        logging.error('..... existed annotated hetSNP with AF is empty, please try again!')
        sys.exit(1)
    else:
        logging.info('..... generated annotated hetSNP with AF file can be found at {0}'.format(hetSNP_AF))

    ##### 1.4 Thinning reads: one reads only count once
    if (not os.path.exists(hetSNP_intersect_unique)) or (not os.path.exists(hetSNP_intersect_unique_forlambda_file)):
        logging.info('=================')
        logging.info('================= Starting specific step 1.4')
        logging.info('..... start filtering variants, we only keep unique reads')
        Intersect_exonicHetSnps(parsed_pileup,hetSNP_AF,read_length,min_total_cov,min_single_cov,hetSNP_intersect_unique,hetSNP_intersect_unique_forlambda_file)
    else:
        logging.info('=================')
        logging.info('================= Skipping specific step 1.4')
    data14_1=pd.read_csv(hetSNP_intersect_unique,sep="\t",header=0,index_col=False)
    data14_2=pd.read_csv(hetSNP_intersect_unique_forlambda_file,sep="\t",header=0,index_col=False)
    if data14_1.shape[1]<2:
        os.remove(hetSNP_intersect_unique)
        logging.error('..... existed hetSNP file with filtered sites file is empty, please try again!')
        sys.exit(1)
    else:
        logging.info('..... generated hetSNP file after intersection with pileup file size {0} with filtered sites can be found at {1}'.format(data14_1.shape[1],hetSNP_intersect_unique))
    if data14_2.shape[0]<2:
        os.remove(hetSNP_intersect_unique)
        logging.error('..... existed hetSNP file with filtered sites prepared for lambda model file is empty, please try again!')
        sys.exit(1)
    else:
        logging.info('..... generated hetSNP file with filtered sites prepared for lambda model can be found at {0}'.format(hetSNP_intersect_unique_forlambda_file))

   ##### 1.5 Annotation LD
    if not os.path.isfile(meta):
        logging.info('=================')
        logging.info('================= Starting specific step 1.5')
        logging.info('..... start annotating LD information')
        annotation.annotateLD(prefix,ancestry,hetSNP_intersect_unique,out,LD_token,chr_start,chr_end,meta)
    else:
        logging.info('=================')
        logging.info('================= Skipping specific step 1.5')
    data14=pd.read_csv(meta,sep="\t",header=0,index_col=False)
    if data14.shape[0]<2:
        os.remove(meta)
        logging.error('..... existed meta file with filtered sites file is empty, please try again!')
        sys.exit(1)
    else:
        logging.info('..... generated annotated meta file can be found at {0}'.format(meta))

    logging.info('================= finish step1! ')
    return hetSNP_intersect_unique,meta,hetSNP_intersect_unique_forlambda_file,hetSNP_intersect_unique_lambdaPredicted_file
