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
annotationDir = config['annoDir']
bam = config['bam']
fastq_path = inputs['fastqPath']
hets_meta_positions = inputs['hetsMetaPositions']
illuminaAdapters = inputs['illuminaAdapterFasta']
mismatch_n = inputs['mismatchN']
model_input_path = inputs['modelInputPath']
model_input = inputs['modelInput']
ref_genome = inputs['refGenome']
sample_name = inputs['sample']
sigma = inputs['sigma']
star_ind = inputs['startIndex']
vcf_dir = inputs['vcfDir']
vcf_file = inputs['vcfFile']

outputs = config['outputs']
rna_pipeline_dir = outputs['rnaPipelineDirectory']
model_output_folder = outputs['modelOutputFolder']
vcf_pipeline_dir = outputs['vcfPipelineDirectory']

programs = config['programs']
picard = programs['picard']
trimmomatic = programs['trimmomatic']

###############################################
# data processing step by step
###############################################
# step0. check data requirement
### check the existence of RNAseq fastq file and VCF file path
### fastq path requirement: fastq path should have files containing '*1.fastq.gz' and '*2.fastq.gz'
rnaSeq1_fastq = False
rnaSeq2_fastq = False

if os.path.exists(fastq_path):
    for filename in os.listdir(fastq_path):
        if '1.fastq' in filename:
            rnaSeq1_fastq = True
        if '2.fastq' in filename:
            rnaSeq2_fastq = True
else:
    print(f"Oops! fastq path ({fastq_path}) doesn't exist. Please try again ...")
    exit(1)

if not (rnaSeq1_fastq or rnaSeq2_fastq):
    print("""Oops! fastq files could not be located. The fastq directory should contain files with names matching
        *1.fastq.gz
        *2.fastq.gz
    Please try again ...""")
    exit(1)

if not os.path.isfile(vcf_file):
    print(f"Oops! VCF file ({vcf_file}) not found. Please try again ...")
    exit(1)

###############################################
# Process RNAseq & VCF
###############################################
#### step1.1 trimming fastq
try:
    cmd = f"./Process_RNAseq/Process_RNAseq_pipeline_I_trim.sh \"{trimmomatic}\" \"{illuminaAdapters}\" \"{sample_name}\" \"{fastq_path}\" \"{rna_pipeline_dir}\""
    check_call(cmd, shell=True)
except CalledProcessError as cpe:
    print(cpe.stderr)
    exit(cpe.returncode)

#### step1.2 RNAseq fastq file alignment
try:
    cmd = f"./Process_RNAseq/Process_RNAseq_pipeline_II_align.sh \"{sample_name}\" \"{ref_genome}\" \"{star_ind}\" \"{mismatch_n}\" \"{picard}\" \"{rna_pipeline_dir}\""
    check_call(cmd, shell=True)
except CalledProcessError as cpe:
    print(cpe.stderr)
    exit(cpe.returncode)

#### step1.3 clean VCF files
try:
    cmd = f"./Process_VCF/Process_VCF_pipeline_I.extractVCF.sh \"{sample_name}\" \"{vcf_dir}\" \"{vcf_pipeline_dir}\""
    check_call(cmd, shell=True)
except CalledProcessError as cpe:
    print(cpe.stderr)
    exit(cpe.returncode)

#### step1.4 extract het sites information and store it in a file for mpileup
cmd = "./Process_VCF/Process_VCF_pipeline_II_hetsMeta.sh %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
check_call(cmd, shell=True)

#### step1.5 mpileup
try:
    cmd = f"./Process_RNAseq/Process_RNAseq_pipeline_III_mpileup.sh \"{sample_name}\" \"{ref_genome}\" \"{mismatch_n}\" \"{hets_meta_positions}\" \"{bam}\" \"{rna_pipeline_dir}\""
    check_call(cmd, shell=True)
except CalledProcessError as cpe:
    print(cpe.stderr)
    exit(cpe.returncode)

###############################################
# Phasing
###############################################
#### step2.1 prepare VCF
cmd = "./Phasing/step1_prepareVCF.sh %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
check_call(cmd, shell=True)

#### step2.2 Phasing
cmd = "./Phasing/step2_phasing.sh %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
check_call(cmd, shell=True)

###############################################
# BEASTIE
###############################################
### (for now) requires user to make the input file in a format that can be used for BEASTIE
#### step3.1 run BEASTIE
cmd = f"python ./BEASTIE/wrapper.py {model_input} {sigma} 0 BEASTIE {model_input_path} {model_output_folder}"
check_call(cmd, shell=True)
