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
        "vcf_file_name",
        "vcf_sample_name",
        "pileup_file_name",
        "simulation_pileup_file_name",
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
        "input_dir",
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
        help="Input is read from <input_dir>/<prefix> and output written to <output_dir>/<prefix>.",
        required=True,
    )
    parser.add_argument(
        "--vcf-filename",
        help="Name of VCF file in <input_dir>/<prefix>.",
        required=True,
    )
    parser.add_argument(
        '--vcf-sample-name',
        help="Name of sample in VCF file (ex HG00096).",
        required=True,
    )
    parser.add_argument(
        '--pileup-filename',
        help="Name of pileup file in <input_dir>/<prefix>.",
        required=True,
    )
    parser.add_argument(
        '--shapeit2-phasing-filename',
        help="Name of the shapeit2 generated phasing data in <input_dir>/<prefix>.",
    )
    parser.add_argument(
        '--simulation-pileup-filename',
        help="Name of simulated data pileup file in <input_dir>/<prefix>.",
    )
    parser.add_argument(
        '--gencode-dir',
        help="Path to gencode reference directory.",
        default=None,
    )
    parser.add_argument(
        '--ancestry',
        help="Ancestry abbreviation (ex GBR, YRI, CEU, ...).",
        required=True
    )
    parser.add_argument(
        "--min-total-cov",
        help="Minimum coverage requirement for total read counts on one site.",
        default=1,
        type=int
    )
    parser.add_argument(
        "--min-single-cov",
        help="Minimum coverage requirement for each REF/ALT allele.",
        default=0,
        type=int
    )
    parser.add_argument(
        "--sigma",
        help="Dispersion parameter for BEASTIE STAN model.",
        default=0.5,
        type=float
    )
    parser.add_argument(
        "--cutoff",
        help="Binomial test p-value cutoff.",
        default=0.05,
        type=float
    )
    parser.add_argument(
        "--alpha",
        help="Type-1 error rate for lambda prediction model.",
        default=0.05,
        type=float
    )
    parser.add_argument(
        '--chr-start',
        help='Starting chromosome number',
        default=1,
        type=int
    )
    parser.add_argument(
        '--chr-end',
        help='Ending chromosome number',
        default=22,
        type=int
    )
    parser.add_argument(
        '--read-length',
        help="Average length of reads for input fastq data.",
        type=int,
        required=True
    )
    parser.add_argument(
        '--ld-token',
        help='LDlink API Token.  Register at https://ldlink.nci.nih.gov/?tab=apiaccess',
        required=True,
    )
    parser.add_argument(
        '--model-name',
        help='Name of stan model to use.',
        default='iBEASTIE2',
    )
    parser.add_argument(
        '--STAN',
        help='Path to STAN model.',
        required=True
    )
    parser.add_argument(
        '--ref-dir',
        help='Path to reference directory containing AF annotation and gencode',
        default="UNSET_REF_DIR",
    )
    parser.add_argument(
        '--input-dir',
        help='Path to directory containing input files',
        required=True,
    )
    parser.add_argument(
        '--save-intermediate',
        help='Keep intermediate files.',
        action='store_true'
    )
    parser.add_argument(
        '--warmup',
        help='Number of warmup estimates to burn in STAN model.',
        default=1000,
        type=int
    ),
    parser.add_argument(
        '--keeper',
        help='Number of estimates to keep from STAN model.',
        default=1000,
        type=int,
    )
    parser.add_argument(
        '--output-dir',
        help='Path to directory where output will be placed.',
        required=True
    )
    parser.add_argument(
        '--ldlink-cache-dir',
        help='Path to directory to save ldlink cache database.',
        default='~/.beastie',
    )

    return parser.parse_args()

def load_config_from_args(args):
    return ConfigurationData(
        prefix=args.prefix,
        vcf_file_name=args.vcf_filename,
        vcf_sample_name=args.vcf_sample_name,
        pileup_file_name=args.pileup_filename,
        shapeit2=args.shapeit2_phasing_filename,
        simulation_pileup_file_name=args.simulation_pileup_filename,
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
        input_dir=args.input_dir,
        SAVE_INT=args.save_intermediate,
        WARMUP=args.warmup,
        KEEPER=args.keeper,
        output_dir=args.output_dir,
        ldlink_cache_dir=os.path.expanduser(
            args.ldlink_cache_dir
        ),
    )

def load_configuration(config_file):
    # If you do not change names of parameters config file, you do not need to modify this run_config.py
    # ====================================================
    config = configparser.ConfigParser()
    config.read(
        config_file
    )  # read in parameters defined in the config you want: parameters_HG_chr20.cfg

    # ==================================================== Nothing has to be changed below this line
    inputs = config["inputs"]
    outputs = config["outputs"]
    config = ConfigurationData(
        prefix=inputs["prefix"],
        vcf_file_name=inputs["vcf_file_name"],
        vcf_sample_name=inputs["vcf_sample_name"],
        pileup_file_name=inputs["pileup_file_name"],
        shapeit2=inputs.get("shapeit2", None),
        simulation_pileup_file_name=inputs.get("simulation_pileup_file_name", None),
        gencode_path=inputs.get("gencode_path", None),
        ancestry=inputs["ancestry"],
        min_total_cov=int(inputs.get("min_total_cov", 1)),
        min_single_cov=int(inputs.get("min_single_cov", 0)),
        sigma=float(inputs.get("sigma", 0.5)),
        cutoff=float(inputs.get("cutoff", 0.5)),
        alpha=float(inputs.get("alpha", 0.05)),
        chr_start=inputs["chr_start"],
        chr_end=inputs["chr_end"],
        read_length=inputs["read_length"],
        LD_token=inputs["LD_token"],
        modelName=inputs["modelName"],
        STAN=inputs["STAN"],
        ref_dir=inputs["ref_dir"],
        input_dir=inputs["input_dir"],
        SAVE_INT=bool(inputs.get("SAVE_INT", False)),
        WARMUP=int(inputs.get("WARMUP", 1000)),
        KEEPER=int(inputs.get("KEEPER", 1000)),
        output_dir=outputs["out_path"],
        ldlink_cache_dir=os.path.expanduser(
            inputs.get("ldlink_cache_dir", "~/.beastie")
        ),
    )

    # ======== Pre-requisite: pre-defined directories and input files and reference files/directories

    return config


###############################################
# BEASTIE
###############################################


def run(config):
    in_path = os.path.join(config.input_dir, config.prefix)
    out_dir = os.path.join(config.output_dir, config.prefix)
    vcf_file = os.path.join(in_path, config.vcf_file_name)
    pileup_file = os.path.join(in_path, config.pileup_file_name)
    shapeit2_file = (
        os.path.join(in_path, config.shapeit2) if config.shapeit2 is not None else None
    )
    simulation_pileup_file = (
        os.path.join(in_path, config.simulation_pileup_file_name)
        if config.simulation_pileup_file_name is not None
        else None
    )
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
        in_path,
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