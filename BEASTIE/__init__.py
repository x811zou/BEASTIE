#!/usr/bin/env python

"""
usage: python run_config.py config_file.cfg
"""
import argparse
from cgitb import handler
import configparser
import logging
import os.path
import sys
from pkg_resources import resource_filename
from collections import namedtuple
from pathlib import Path
from datetime import date
import os
from . import runModel

ConfigurationData = namedtuple(
    "ConfigurationData",
    [
        "prefix",
        "vcfgz_file",
        "vcf_sample_name",
        "shapeit2",
        "collected_alignmentBias_file",
        "simulation_pileup_file",
        "filtered_het_snp_file",
        "ancestry",
        "min_total_cov",
        "min_single_cov",
        "read_length",
        "sigma",
        "alignBiasP_cutoff",
        "ase_cutoff",
        "alpha",
        "chr_start",
        "chr_end",
        "LD_token",
        "modelName",
        "nophasing",
        "atacseq",
        "STAN",
        "SAVE_INT",
        "WARMUP",
        "KEEPER",
        "output_dir",
        "ldlink_cache_dir",
        "ldlink_token_db",
        "gam_model_name",
    ],
)


def load_config_from_args(args):
    if not args.ld_token and not args.ldlink_token_db:
        print("ERROR: ld-token or ldlink-token-db required")
        sys.exit(1)

    return ConfigurationData(
        prefix=args.prefix if args.prefix is not None else args.vcf_sample_name,
        vcfgz_file=args.vcfgz_file,
        vcf_sample_name=args.vcf_sample_name,
        shapeit2=args.shapeit2_phasing_file,
        simulation_pileup_file=args.simulation_pileup_file,
        collected_alignmentBias_file=args.collected_alignmentBias_file,
        filtered_het_snp_file=args.filtered_het_snp_file,
        ancestry=args.ancestry,
        min_total_cov=args.min_total_cov,
        min_single_cov=args.min_single_cov,
        read_length=args.read_length,
        sigma=args.sigma,
        alignBiasP_cutoff=args.alignBiasP_cutoff,
        ase_cutoff=args.ase_cutoff,
        alpha=args.alpha,
        chr_start=args.chr_start,
        chr_end=args.chr_end,
        LD_token=args.ld_token,
        modelName=args.model,
        nophasing=args.nophasing,
        atacseq=args.atacseq,
        STAN=args.STAN,
        SAVE_INT=args.save_intermediate,
        WARMUP=args.warmup,
        KEEPER=args.keeper,
        output_dir=args.output_dir,
        ldlink_cache_dir=os.path.expanduser(args.ldlink_cache_dir),
        ldlink_token_db=os.path.expanduser(args.ldlink_token_db) if args.ldlink_token_db else None,
        gam_model_name=args.gam_model_name,
    )


###############################################
# BEASTIE
###############################################


def run(config):
    output_path = config.output_dir
    today = date.today()
    specification = f"chr{config.chr_start}-{config.chr_end}_alignBiasp{config.alignBiasP_cutoff}_s{config.sigma}_a{config.alpha}_sinCov{config.min_single_cov}_totCov{config.min_total_cov}_W{config.WARMUP}K{config.KEEPER}"
    specification_path = os.path.join(output_path, specification)
    log_path = os.path.join(specification_path, "log")
    tmp_path = os.path.join(specification_path, "tmp")
    
    Path(output_path).mkdir(parents=True, exist_ok=True)
    Path(specification_path).mkdir(parents=True, exist_ok=True)
    Path(log_path).mkdir(parents=True, exist_ok=True)
    Path(tmp_path).mkdir(parents=True, exist_ok=True)

    if config.atacseq is True:
        atacseq=True
        nophasing = True
        modelName = "BEASTIE3-fix-uniform"
        ancestry = None
    else:
        modelName = config.modelName
        nophasing = config.nophasing
        ancestry == config.ancestry
    result_path = os.path.join(specification_path, modelName)
    Path(result_path).mkdir(parents=True, exist_ok=True)
    log_filename = f"{config.prefix}-{today.strftime('%b-%d-%Y')}-{modelName}"
    # stdout_stderr_filepath = os.path.join(log_path, f"{log_filename}.output")
    # stdout_stderr_file = open(stdout_stderr_filepath, "w")
    # sys.stdout = stdout_stderr_file
    # sys.stderr = stdout_stderr_file

    logname = os.path.join(log_path, f"{log_filename}.log")
    if os.path.isfile(logname):
        os.remove(logname)

    file_handler = logging.FileHandler(logname, mode="w", delay=False)
    stream_handler = logging.StreamHandler(sys.stdout)
    logging.basicConfig(
        # filename=logname,
        # filemode="a",
        handlers=[file_handler, stream_handler],
        format="%(asctime)-15s [%(levelname)s] %(message)s",
        level=logging.DEBUG,
    )

    logging.info(">> Starting running BEASTIE")
    logging.info("======================================== ")
    logging.info(
        "======================================== Processing simulation raw data & Preparing input in a format required for BEASTIE model"
    )
    logging.info("======================================== ")
    runModel.run(
        config.prefix,
        config.vcfgz_file,
        config.vcf_sample_name,
        config.collected_alignmentBias_file,
        config.simulation_pileup_file,
        config.filtered_het_snp_file,
        output_path,
        tmp_path,
        result_path,
        config.shapeit2,
        config.alignBiasP_cutoff,
        config.ase_cutoff,
        os.path.join(config.STAN, modelName),
        config.ancestry,
        config.chr_start,
        config.chr_end,
        config.min_total_cov,
        config.min_single_cov,
        config.read_length,
        config.alpha,
        config.sigma,
        config.SAVE_INT,
        nophasing,
        atacseq,
        config.WARMUP,
        config.KEEPER,
        config.LD_token,
        config.ldlink_cache_dir,
        config.ldlink_token_db,
        config.gam_model_name,
    )
