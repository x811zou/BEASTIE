#!/usr/bin/env python

"""
usage: python run_config_HG_chr20.py
"""
import configparser
import argparse
import os
import os.path
from subprocess import check_call, CalledProcessError

# If you do not change names of parameters config file, you do not need to modify this run_config.py
#====================================================
config = configparser.ConfigParser()
config.read('parameters_HG_chr20.cfg')               # read in parameters defined in the config you want: parameters_HG_chr20.cfg

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
pileup_file=in_path+pileup_file_name
model = STAN+str(modelName)+"/"+str(modelName)

###############################################
# BEASTIE
###############################################

cmd = f"python BEASTIE.py build --prefix {prefix} --vcf_sample_name {vcf_sample_name} --ref_dir {ref_dir} --vcf {vcf_file} --pileup {pileup_file} --in_path {in_path} --ancestry {ancestry} --chr_start {chr_start} --chr_end {chr_end} --read_length {read_length} --LD_token {LD_token} --model {model} "
os.system(cmd)