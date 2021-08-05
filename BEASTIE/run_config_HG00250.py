#!/usr/bin/env python

"""
usage: python run_config_HG00250.py
"""
import configparser
import argparse
import os
import os.path
import subprocess

# If you do not change names of parameters config file, you do not need to modify this run_config.py
#====================================================
config = configparser.ConfigParser()
config.read('parameters_HG00250.cfg')               # read in parameters defined in the config you want: parameters_HG_chr20.cfg

#==================================================== No thing has to be changed below this line
inputs = config['inputs']
outputs = config['outputs']

[inputs]
prefix=inputs['prefix']
vcf_file_name=inputs['vcf_file_name']
vcf_sample_name=inputs['vcf_sample_name']
pileup_file_name=inputs['pileup_file_name']
ancestry=inputs['ancestry']
min_total_cov = inputs['min_total_cov']
min_single_cov = inputs['min_single_cov']
sigma = inputs['sigma']
cutoff=inputs['cutoff']
alpha=inputs['alpha']
chr_start=inputs['chr_start']
chr_end=inputs['chr_end']
read_length=inputs['read_length']
LD_token=inputs['LD_token']
modelName=inputs['modelName']
STAN=inputs['STAN']
work_dir=inputs['work_dir']
ref_dir=inputs['ref_dir']
in_path=inputs['input_dir']
SAVE_INT=inputs['SAVE_INT']

[outputs]
out_path=outputs['out_path']

#==================================================== Pre-requisite: pre-defined directories and input files and reference files/directories
in_path=in_path+str(prefix)+"/"
out_dir=out_path+str(prefix)+"/"
tmp_dir=out_dir+"TEMP/"
vcf_file=in_path+vcf_file_name
#CHROM   | POS  |    ID     |  REF   |  ALT  |  QUAL |  FILTER |     INFO      |           FORMAT             |                HG001
#1       837214  rs72631888      G       C       50      PASS    platforms=3;   GT:DP:ADALL:AD:GQ:IGT:IPS:PS    0|1:558:122,138:135,166:581:0/1:.:PATMAT

pileup_file=in_path+pileup_file_name
#chr10	323283	A	1	g	E	~
model = STAN+str(modelName)+"/"+str(modelName)

from datetime import date
today = date.today()
stdout_file=in_path+"output/"+str(prefix)+"-"+str(today.strftime("%b-%d-%Y"))+".stdout"
# #==================================================== Intermediate files we will generate (with default names), or users can generate using their own modifications (please modify file names here if needed)
# hetSNP_file=str(prefix)+"_hetSnps.tsv"
# #chr   | chrN  |     geneID     | genomicCoord_pos | transcriptCoord |    SNP_id |  genotype
# #chr1    1       ENSG00000227232.4       14930             -1         rs75454623      1|0  
# meta_file=str(prefix)+"_logisticReg_input.tsv"
# parsed_pileup_file=str(prefix)+"_parsed_pileup.tsv"

#==================================================== check the existence of required input files


# if not os.path.isfile(pileup_file):
#     print(f"Oops! samtools pileup file ({pileup_file}) not found. Please try again ...")
#     exit(1)

# if not os.path.isfile(AF_file):
#     print(f"Oops! reference AF file ({AF_file}) not found. Please try again ...")
#     exit(1)

# if not os.path.exists(gencode_ref):
#     print(f"Oops! reference gencode directory ({gencode_ref}) not found. Please try again ...")
#     exit(1)

# # output from step1.1
# if not os.path.isfile(hetSNP_file):
#     print(f"..... Alright, het SNP file ({hetSNP_file}) not found. We will generate that for you ...")

# # output from step1.3
# if not os.path.isfile(meta_file):
#     print(f"..... Alright, logistic regression input file ({meta_file}) not found. We will generate that for you ...")

# # output from step1.2
# if not os.path.isfile(parsed_pileup_file):
#     print(f"..... Alright, parsed pileup file ({parsed_pileup_file}) not found. We will generate that for you ...")

###############################################
# BEASTIE
###############################################
cmd = f"python BEASTIE.py build --prefix {prefix} --vcf_sample_name {vcf_sample_name} --ref_dir {ref_dir} --vcf {vcf_file} --pileup {pileup_file} --in_path {in_path} --ancestry {ancestry} --chr_start {chr_start} --chr_end {chr_end} --read_length {read_length} --LD_token {LD_token} --model {model}"
#print(cmd)
#print(stdout_file)
subprocess.call(cmd,shell=True,stdout=open(stdout_file, 'w'), stderr=open(stdout_file, 'w'))