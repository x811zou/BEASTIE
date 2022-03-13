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
        "-c",
        "--config",
        help="File containing parameters for running BEASTIE",
        required=True,
    )

    return parser.parse_args()


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


if __name__ == "__main__":
    args = check_arguments()
    config = load_configuration(args.config)
    run(config)
