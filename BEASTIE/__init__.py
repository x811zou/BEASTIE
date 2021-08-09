#!/usr/bin/env python

"""
usage: python run_config.py config_file.cfg
"""
import argparse
import configparser
import os.path
import subprocess
import sys

from collections import namedtuple
from datetime import date

ConfigurationData = namedtuple("ConfigurationData", [
    "prefix", "vcf_file_name", "vcf_sample_name", "pileup_file_name", "ancestry", "min_total_cov", "min_single_cov", "sigma", "cutoff", "alpha", "chr_start", "chr_end", "read_length", "LD_token", "modelName", "STAN", "work_dir", "ref_dir", "input_dir", "SAVE_INT", "WARMUP", "KEEPER", "output_dir"
])

def check_arguments():
    parser = argparse.ArgumentParser("beastie", description="Bayesian Estimation of Allele Specific Transcription Integrating across Exons")
    parser.add_argument("-c", "--config", help="File containing parameters for running BEASTIE", required=True)

    return parser.parse_args()


def load_configuration(config_file):
    # If you do not change names of parameters config file, you do not need to modify this run_config.py
    #====================================================
    config = configparser.ConfigParser()
    config.read(config_file)               # read in parameters defined in the config you want: parameters_HG_chr20.cfg

    #==================================================== Nothing has to be changed below this line
    inputs = config['inputs']
    outputs = config['outputs']


    config = ConfigurationData(
        prefix=inputs['prefix'],
        vcf_file_name=inputs['vcf_file_name'],
        vcf_sample_name=inputs['vcf_sample_name'],
        pileup_file_name=inputs['pileup_file_name'],
        ancestry=inputs['ancestry'],
        min_total_cov=inputs['min_total_cov'],
        min_single_cov=inputs['min_single_cov'],
        sigma=inputs['sigma'],
        cutoff=inputs['cutoff'],
        alpha=inputs['alpha'],
        chr_start=inputs['chr_start'],
        chr_end=inputs['chr_end'],
        read_length=inputs['read_length'],
        LD_token=inputs['LD_token'],
        modelName=inputs['modelName'],
        STAN=inputs['STAN'],
        work_dir=inputs['work_dir'],
        ref_dir=inputs['ref_dir'],
        input_dir=inputs['input_dir'],
        SAVE_INT=inputs['SAVE_INT'],
        WARMUP=inputs['WARMUP'],
        KEEPER=inputs['KEEPER'],
        output_dir=outputs['out_path'],
    )

    #======== Pre-requisite: pre-defined directories and input files and reference files/directories

    return config

###############################################
# BEASTIE
###############################################

def run(config):
    in_path = os.path.join(config.input_dir, config.prefix)
    vcf_file = os.path.join(in_path, config.vcf_file_name)
    pileup_file = os.path.join(in_path, config.pileup_file_name)
    model = os.path.join(f"{config.STAN}{config.modelName}", config.modelName)
    stdout_file = os.path.join(in_path, "output", f"{config.prefix}-{date.today().strftime('%b-%d-%Y')}.stdout")

    cmd = f"python BEASTIE.py build --prefix {config.prefix} --vcf_sample_name {config.vcf_sample_name} --ref_dir {config.ref_dir} --vcf {vcf_file} --pileup {pileup_file} --in_path {in_path} --ancestry {config.ancestry} --chr_start {config.chr_start} --chr_end {config.chr_end} --read_length {config.read_length} --LD_token {config.LD_token} --model {model} --cutoff {config.cutoff} --alpha {config.alpha} --WARMUP {config.WARMUP} --KEEPER {config.KEEPER}"
    subprocess.call(cmd, shell=True, stdout=open(stdout_file, 'w'), stderr=open(stdout_file, 'w'))

if __name__ == "__main__":
    args = check_arguments()
    config = load_configuration(args.config)
    run(config)
