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
from . import extractHets, beastie_step1, beastie_step2

ConfigurationData = namedtuple(
    "ConfigurationData",
    [
        "prefix",
        "vcfgz_file",
        "vcf_sample_name",
        "pileup_file",
        "simulation_pileup_file",
        "shapeit2",
        "het_snp_file",
        "gencode_path",
        "af_path",
        "ancestry",
        "min_total_cov",
        "min_single_cov",
        "sigma",
        "cutoff",
        "alpha",
        "chr_start",
        "chr_end",
        "include_x_chromosome",
        "read_length",
        "LD_token",
        "modelName",
        "STAN",
        "SAVE_INT",
        "WARMUP",
        "KEEPER",
        "output_dir",
        "ldlink_cache_dir",
        "ldlink_token_db",
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
        pileup_file=args.pileup_file,
        shapeit2=args.shapeit2_phasing_file,
        simulation_pileup_file=args.simulation_pileup_file,
        het_snp_file=args.het_snp_file,
        gencode_path=args.gencode_dir,
        af_path=args.af_dir,
        ancestry=args.ancestry,
        min_total_cov=args.min_total_cov,
        min_single_cov=args.min_single_cov,
        sigma=args.sigma,
        cutoff=args.cutoff,
        alpha=args.alpha,
        chr_start=args.chr_start,
        chr_end=args.chr_end,
        include_x_chromosome=args.include_x_chromosome,
        read_length=args.read_length,
        LD_token=args.ld_token,
        modelName=args.model_name,
        STAN=args.STAN,
        SAVE_INT=args.save_intermediate,
        WARMUP=args.warmup,
        KEEPER=args.keeper,
        output_dir=args.output_dir,
        ldlink_cache_dir=os.path.expanduser(args.ldlink_cache_dir),
        ldlink_token_db=os.path.expanduser(args.ldlink_token_db)
        if args.ldlink_token_db
        else None,
    )


###############################################
# BEASTIE
###############################################


def run(config):
    out_dir = config.output_dir
    vcfgz_file = config.vcfgz_file
    pileup_file = config.pileup_file
    shapeit2_file = config.shapeit2
    simulation_pileup_file = config.simulation_pileup_file
    gencode_path = config.gencode_path
    af_path = config.af_path
    model = os.path.join(config.STAN, config.modelName)
    today = date.today()

    specification = f"chr{config.chr_start}-{config.chr_end}_s{config.sigma}_a{config.alpha}_sinCov{config.min_single_cov}_totCov{config.min_total_cov}_W{config.WARMUP}K{config.KEEPER}"
    output_path = out_dir
    specification_path = os.path.join(output_path, specification)
    log_path = os.path.join(specification_path, "log")
    tmp_path = os.path.join(specification_path, "tmp")
    result_path = os.path.join(specification_path, "result")

    Path(output_path).mkdir(parents=True, exist_ok=True)
    Path(specification_path).mkdir(parents=True, exist_ok=True)
    Path(log_path).mkdir(parents=True, exist_ok=True)
    Path(tmp_path).mkdir(parents=True, exist_ok=True)
    Path(result_path).mkdir(parents=True, exist_ok=True)

    log_filename = f"{config.prefix}-{today.strftime('%b-%d-%Y')}"

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

    logging.info("========================================")
    logging.info(
        "======================================== step1: Processing raw data & annotating LD and AF information"
    )
    logging.info("======================================== ")
    hetSNP_intersect_unique, hetSNP_intersect_unique_sim = beastie_step1.run(
        config.prefix,
        config.vcf_sample_name,
        output_path,
        tmp_path,
        gencode_path,
        af_path,
        model,
        vcfgz_file,
        config.ancestry,
        config.chr_start,
        config.chr_end,
        config.min_total_cov,
        config.min_single_cov,
        config.read_length,
        pileup_file,
        simulation_pileup_file,
        shapeit2_file,
        config.het_snp_file,
    )
    logging.info("======================================== ")
    logging.info(
        "======================================== step2: Preparing input in a format required for BEASTIE model"
    )
    logging.info("======================================== ")
    beastie_step2.run(
        shapeit2_file,
        hetSNP_intersect_unique,
        hetSNP_intersect_unique_sim,
        config.prefix,
        config.alpha,
        model,
        config.sigma,
        tmp_path,
        result_path,
        config.cutoff,
        config.SAVE_INT,
        config.WARMUP,
        config.KEEPER,
        config.min_total_cov,
        config.min_single_cov,
        config.ancestry,
        config.LD_token,
        config.ldlink_cache_dir,
        config.ldlink_token_db,
    )
