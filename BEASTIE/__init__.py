#!/usr/bin/env python

"""
usage: python run_config.py config_file.cfg
"""
import argparse
import configparser
import logging
import os.path
from collections import namedtuple
from datetime import date

from . import beastie_step1, beastie_step2

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
        cutoff=float(inputs['cutoff']),
        alpha=float(inputs['alpha']),
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
    today = date.today()

    logname = os.path.join(in_path, "output", f"{config.prefix}-{today.strftime('%b-%d-%Y')}.log")
    if os.path.isfile(logname):
        os.remove(logname)
    logging.basicConfig(filename=logname,
                            filemode='a',
                            format='%(asctime)-15s [%(levelname)s] %(message)s',
                            level=logging.DEBUG)

    logging.info(">> Starting running BEASTIE")

    logging.info('========================================')
    logging.info('======================================== step1: Processing raw data & annotating LD and AF information')
    logging.info('======================================== ')
    hetSNP_intersect_unique, meta, hetSNP_intersect_unique_forlambda_file, hetSNP_intersect_unique_lambdaPredicted_file = beastie_step1.run(
        0.5,  # sigma
        config.alpha,
        config.WARMUP,
        config.KEEPER,
        config.prefix,
        config.vcf_sample_name,
        in_path,
        "output",  # out
        model,
        vcf_file,
        config.ref_dir,
        config.ancestry,
        config.chr_start,
        config.chr_end,
        1,  # min_total_cov
        0,  # min_single_cov
        config.read_length,
        config.LD_token,
        pileup_file,
        None,  # hetSNP
        None  # parsed_pileup
    )

    logging.info('======================================== ')
    logging.info('======================================== step2: Preparing input in a format required for BEASTIE model')
    logging.info('======================================== ')
    beastie_step2.run(
        hetSNP_intersect_unique,
        meta,
        hetSNP_intersect_unique_forlambda_file,
        hetSNP_intersect_unique_lambdaPredicted_file,
        config.prefix,
        config.alpha,
        model,
        0.5,  # sigma
        in_path,
        "output",  # out
        config.cutoff,
        "False",  # SAVE_INT
        config.WARMUP,
        config.KEEPER,
        1,  # min_total_cov
        0,  # min_single_cov
    )

if __name__ == "__main__":
    args = check_arguments()
    config = load_configuration(args.config)
    run(config)
