#!/usr/bin/env python

###############################################
# usage: python run_ConfigFile.py
###############################################

import configparser
import os
import os.path
from subprocess import check_call, CalledProcessError

###############################################
# read in parameters defined in parameters.cfg
###############################################
config = configparser.ConfigParser()
config.read('parameters.cfg')

inputs = config['inputs']
outputs = config['outputs']
programs = config['programs']

min_total_cov = inputs['min_total_cov']
min_single_cov = inputs['min_single_cov']
sigma = inputs['sigma']
cutoff=inputs['cutoff']
alpha=inputs['alpha']
prefix=inputs['prefix']
modelName=inputs['modelName']
STAN=inputs['STAN']
local_dir=inputs['local_dir']
vcfgz_file=inputs['vcfgz_file']
pileup_file=inputs['pileup_file']
meta_file=inputs['meta_file']

[outputs]
out=outputs['out']

[programs]
model = programs['model']

###############################################
# data processing step by step
###############################################
# step0. check data requirement
### check the existence of RNAseq fastq file and VCF file path
### fastq path requirement: fastq path should have files containing '*1.fastq.gz' and '*2.fastq.gz'
stan_model = False

if os.path.exists(STAN):
    for filename in os.listdir(STAN):
        if '.stan' in filename:
            stan_model = True
else:
    print(f"Oops! STAN path ({STAN}) doesn't exist. Please try again ...")
    exit(1)

if not os.path.isfile(vcfgz_file):
    print(f"Oops! VCF file ({vcfgz_file}) not found. Please try again ...")
    exit(1)

# if not os.path.isfile(pileup_file):
#     print(f"Oops! samtools pileup file ({pileup_file}) not found. Please try again ...")
#     exit(1)

# if not os.path.isfile(meta_file):
#     print(f"Oops! meta data file ({meta_file})2878 not found. Please try again ...")
#     exit(1)
###############################################
# BEASTIE
###############################################
try:
    cmd = f"python BEASTIE.py {local_dir} {vcfgz_file} {pileup_file} {meta_file} {min_total_cov} {min_single_cov} {prefix} {out} {modelName} {STAN} {alpha} {sigma} {cutoff}"
    # cmd = f"./Phasing/step1_prepareVCF.sh \"{sample_name}\" \"{phasing_vcf_file}\" \"{phasing_vcf_dir}\""  
    check_call(cmd, shell=True)
except CalledProcessError as cpe:
    print(cpe.stderr)
    exit(cpe.returncode)