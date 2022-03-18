#!/usr/bin/env python

"""
usage: python run_config.py config_file.cfg
"""
import argparse
import configparser
import logging
import os.path
import sys
from pkg_resources import resource_filename
from collections import namedtuple
from pathlib import Path
from datetime import date
from . import beastie_step1, beastie_step2

ConfigurationData = namedtuple(
    "ConfigurationData",
    [
        "prefix",
        "vcf_file",
        "vcf_sample_name",
        "pileup_file",
        "simulation_pileup_file",
        "shapeit2",
        "gencode_path",
        "ancestry",
        "min_total_cov",
        "min_single_cov",
        "sigma",
        "cutoff",
        "alpha",
        "chr_start",
        "chr_end",
        "read_length",
        "LD_token",
        "modelName",
        "STAN",
        "ref_dir",
        "SAVE_INT",
        "WARMUP",
        "KEEPER",
        "output_dir",
        "ldlink_cache_dir",
    ],
)


def check_arguments():
    parser = argparse.ArgumentParser(
        "beastie",
        description="Bayesian Estimation of Allele Specific Transcription Integrating across Exons",
    )
    parser.add_argument(
        "--prefix",
        help="Prefix for output files, default is --vcf-sample-name value.",
    )
    parser.add_argument(
        "--vcf-file",
        help="Path to VCF file",
        required=True,
    )
    parser.add_argument(
        "--vcf-sample-name",
        help="Name of sample in VCF file (ex HG00096).",
        required=True,
    )
    parser.add_argument(
        "--pileup-file",
        help="Path to pileup file.",
        required=True,
    )
    parser.add_argument(
        "--shapeit2-phasing-file",
        help="Path to shapeit2 generated phasing data.",
    )
    parser.add_argument(
        "--simulation-pileup-file",
        help="Path to simulated data pileup file.",
    )
    parser.add_argument(
        "--gencode-dir",
        help="Path to gencode reference directory.",
        default=None,
    )
    parser.add_argument(
        "--ancestry",
        help="Ancestry abbreviation (ex YRI, CEU, FIN ..).",
        required=True,
    )
    parser.add_argument(
        "--min-total-cov",
        help="Minimum coverage requirement for total read counts on one site.",
        default=1,
        type=int,
    )
    parser.add_argument(
        "--min-single-cov",
        help="Minimum coverage requirement for each REF/ALT allele.",
        default=0,
        type=int,
    )
    parser.add_argument(
        "--sigma",
        help="Dispersion parameter for BEASTIE STAN model.",
        default=0.5,
        type=float,
    )
    parser.add_argument(
        "--cutoff", help="Binomial test p-value cutoff.", default=0.05, type=float
    )
    parser.add_argument(
        "--alpha",
        help="Type-1 error rate for lambda prediction model.",
        default=0.05,
        type=float,
    )
    parser.add_argument(
        "--chr-start", help="Starting chromosome number", default=1, type=int
    )
    parser.add_argument(
        "--chr-end", help="Ending chromosome number", default=22, type=int
    )
    parser.add_argument(
        "--read-length",
        help="Average length of reads for input fastq data.",
        type=int,
        required=True,
    )
    parser.add_argument(
        "--ld-token",
        help="LDlink API Token.  Register at https://ldlink.nci.nih.gov/?tab=apiaccess",
        required=True,
    )
    parser.add_argument(
        "--model-name",
        help="Name of stan model to use.",
        default="iBEASTIE2",
    )
    parser.add_argument("--STAN", help="Path to STAN model.", required=True)
    parser.add_argument(
        "--ref-dir",
        help="Path to reference directory containing AF annotation and gencode",
        default="UNSET_REF_DIR",
    )
    parser.add_argument(
        "--save-intermediate", help="Keep intermediate files.", action="store_true"
    )
    parser.add_argument(
        "--warmup",
        help="Number of warmup estimates to burn in STAN model.",
        default=1000,
        type=int,
    ),
    parser.add_argument(
        "--keeper",
        help="Number of estimates to keep from STAN model.",
        default=1000,
        type=int,
    )
    parser.add_argument(
        "--output-dir",
        help="Path to directory where output will be placed.",
        required=True,
    )
    parser.add_argument(
        "--ldlink-cache-dir",
        help="Path to directory to save ldlink cache database.",
        default="~/.beastie",
    )

    return parser.parse_args()


def load_config_from_args(args):
    return ConfigurationData(
        prefix=args.prefix if args.prefix is not None else args.vcf_sample_name,
        vcf_file=args.vcf_file,
        vcf_sample_name=args.vcf_sample_name,
        pileup_file=args.pileup_file,
        shapeit2=args.shapeit2_phasing_file,
        simulation_pileup_file=args.simulation_pileup_file,
        gencode_path=args.gencode_dir,
        ancestry=args.ancestry,
        min_total_cov=args.min_total_cov,
        min_single_cov=args.min_single_cov,
        sigma=args.sigma,
        cutoff=args.cutoff,
        alpha=args.alpha,
        chr_start=args.chr_start,
        chr_end=args.chr_end,
        read_length=args.read_length,
        LD_token=args.ld_token,
        modelName=args.model_name,
        STAN=args.STAN,
        ref_dir=args.ref_dir,
        SAVE_INT=args.save_intermediate,
        WARMUP=args.warmup,
        KEEPER=args.keeper,
        output_dir=args.output_dir,
        ldlink_cache_dir=os.path.expanduser(args.ldlink_cache_dir),
    )

    # ======== Pre-requisite: pre-defined directories and input files and reference files/directories

    return config


###############################################
# BEASTIE
###############################################


def run(config):
    out_dir = config.output_dir
    vcf_file = config.vcf_file
    pileup_file = config.pileup_file
    shapeit2_file = config.shapeit2
    simulation_pileup_file = config.simulation_pileup_file
    gencode_path = (
        config.gencode_path
        if config.gencode_path is not None
        else resource_filename("BEASTIE", "reference/gencode_chr")
        # else config.ref_dir + "/gencode_chr"
    )
    model = os.path.join(config.STAN, config.modelName)
    today = date.today()

    specification = f"s{config.sigma}_a{config.alpha}_sinCov{config.min_single_cov}_totCov{config.min_total_cov}_W{config.WARMUP}K{config.KEEPER}"
    output_path = os.path.join(out_dir, "output")
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

    stdout_stderr_filepath = os.path.join(log_path, f"{log_filename}.output")
    stdout_stderr_file = open(stdout_stderr_filepath, "w")
    sys.stdout = stdout_stderr_file
    sys.stderr = stdout_stderr_file

    logname = os.path.join(log_path, f"{log_filename}.log")
    if os.path.isfile(logname):
        os.remove(logname)
    logging.basicConfig(
        filename=logname,
        filemode="a",
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
        model,
        vcf_file,
        config.ref_dir,
        config.ancestry,
        config.chr_start,
        config.chr_end,
        config.min_total_cov,
        config.min_single_cov,
        config.read_length,
        pileup_file,
        simulation_pileup_file,
        shapeit2_file,
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
        config.chr_start,
        config.chr_end,
        config.ldlink_cache_dir,
    )
