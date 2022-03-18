#!/usr/bin/env python
# =========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
# =========================================================================
import os
import sys
import logging
from unittest.mock import NonCallableMagicMock
import pandas as pd
from pathlib import Path
from pkg_resources import resource_filename
from .annotation import annotateAF
from .extractHets import count_all_het_sites
from .helpers import runhelper
from .intersect_hets import Intersect_exonicHetSnps
from .parse_mpileup import Parse_mpileup_allChr


def check_file_existence(model, vcfgz, ref_dir, pileup, simulation_pileup, shapeit2):
    ##### STAN model
    if False and not os.path.exists(model):
        logging.error(
            "Oops! STAN model {0} doesn't exist. Please try again ...".format(model)
        )
        sys.exit(1)
    else:
        logging.info("Great! STAN model {0} exists.".format(model))
    vcfgztbi = "{0}.tbi".format(vcfgz)
    if not os.path.isfile(vcfgz):
        logging.error(
            "Oops! vcfgz file {0} doesn't exist. Please try again ...".format(vcfgz)
        )
        exit(1)
    if not os.path.isfile(vcfgztbi):
        logging.warning("We will generate vcfgz index for you ...")
        cmd = "tabix -fp vcf %s" % (vcfgz)
        runhelper(cmd)

    ##### reference directory : AF file and gencode directory
    ref_dir = resource_filename("BEASTIE", "reference/")
    if not os.path.exists(ref_dir):
        logging.error(
            "Oops! REFERENCE path {0} doesn't exist. Please try again ...".format(
                ref_dir
            )
        )
        sys.exit(1)
    else:
        logging.info("Great! REFERENCE path {0} exists.".format(ref_dir))
    ##### pileup_file
    if not os.path.isfile(pileup):
        logging.error(
            "Oops!  pileup file {0} doesn't exist. Please try again ...".format(pileup)
        )
    else:
        logging.info("Great! pileup file {0} exists.".format(pileup))
    ##### simulation pileup_file
    if simulation_pileup is not None:
        if not os.path.isfile(simulation_pileup):
            logging.error(
                "Oops!  simulation pileup file {0} doesn't exist. Please try again ...".format(
                    simulation_pileup
                )
            )
        else:
            logging.info("Great! simulation pileup file {0} exists.".format(pileup))

    ##### shapeit2_file
    if shapeit2 is not None:
        if not os.path.isfile(shapeit2):
            logging.error(
                "Oops!  shapeit2 file {0} doesn't exist. Please try again ...".format(
                    shapeit2
                )
            )
        else:
            logging.info("Great! shapeit2 file {0} exists.".format(shapeit2))


def run(
    prefix,
    vcf_sample_name,
    output_path,
    tmp_path,
    gencode_path,
    model,
    vcfgz,
    ref_dir,
    ancestry,
    chr_start,
    chr_end,
    min_total_cov,
    min_single_cov,
    read_length,
    pileup,
    simulation_pileup,
    shapeit2_file,
):
    #####
    ##### 0.0 Check input file existence
    #####
    check_file_existence(
        model, vcfgz, ref_dir, pileup, simulation_pileup, shapeit2_file
    )

    #####
    ##### 1.1 Generate hetSNP file: extract heterozygous bi-allelic SNPs for specific chromosomes from all gencode transcripts
    #####
    hetSNP = os.path.join(output_path, f"{prefix}_hetSNP_chr{chr_start}-{chr_end}.tsv")
    logging.info("=================")
    logging.info("================= Starting common step 1.1")
    if not os.path.exists(hetSNP):
        logging.info("....... start extracting heterozygous bi-allelic SNPs from VCF")
        count_all_het_sites(
            vcfgz,
            hetSNP,
            int(chr_start),
            int(chr_end),
            gencode_path,
        )
    else:
        logging.info("================= Skipping common step 1.1")
        logging.info("=================")
    data11 = pd.read_csv(hetSNP, sep="\t", header=0, index_col=False)
    logging.debug(
        "output {0} has {1} het SNPs from VCF file".format(
            os.path.basename(hetSNP), data11.shape[0]
        )
    )
    if data11.shape[0] < 2:
        os.remove(hetSNP)
        logging.error("..... existed hetSNP file is empty, please try again!")
        sys.exit(1)
    else:
        logging.info(
            "....... {0} save to {1}".format(
                os.path.basename(hetSNP), os.path.dirname(hetSNP)
            )
        )

    #####
    ##### 1.2 Annotation: AF
    #####
    logging.info("=================")
    logging.info("================= Starting common step 1.2")
    hetSNP_AF = f"{os.path.splitext(hetSNP)[0]}_AF.tsv"
    if os.path.isfile(hetSNP_AF):
        logging.info("================= Skipping common step 1.2")
        logging.info("=================")
    else:
        logging.info("..... start annotating hetSNP from 1.1 with AF from 1000 Genome")
        annotateAF(ancestry, hetSNP, hetSNP_AF)

    data13 = pd.read_csv(hetSNP_AF, sep="\t", header=0, index_col=False)
    logging.debug(
        "output {0} annotates {1} het SNPs with AF annotation".format(
            os.path.basename(hetSNP_AF), data13.shape[0]
        )
    )
    if data13.shape[0] < 2:
        os.remove(hetSNP_AF)
        logging.error(
            "....... existing annotated hetSNP with AF is empty, please try again!"
        )
        sys.exit(1)
    else:
        logging.info(
            "....... {0} save to {1}".format(
                os.path.basename(hetSNP_AF), os.path.dirname(hetSNP_AF)
            )
        )

    #####
    ##### 1.3 Generate parsed pileup file : parse samtools mpile up output files
    #####
    parsed_pileup = os.path.join(
        output_path, f"{prefix}_parsed_pileup_chr{chr_start}-{chr_end}.tsv"
    )
    logging.info("=================")
    logging.info("================= Starting common step 1.3")
    # required data
    if simulation_pileup is None:
        simulation_parsed_pileup = None
        if not os.path.exists(parsed_pileup):
            logging.info("....... start parsing samtools pileup result for input data")
            Parse_mpileup_allChr(
                vcf_sample_name,
                vcfgz,
                pileup,
                min_total_cov,
                min_single_cov,
                parsed_pileup,
            )
        else:
            logging.info("================= Skipping common step 1.3")
            logging.info("=================")
    # required data + simulation data
    else:
        simulation_parsed_pileup = os.path.join(
            output_path,
            f"{prefix}_parsed_pileup_chr{chr_start}-{chr_end}.simulation.tsv",
        )
        if not os.path.exists(parsed_pileup):
            logging.info("....... start parsing samtools pileup result for input data")
            Parse_mpileup_allChr(
                vcf_sample_name,
                vcfgz,
                pileup,
                min_total_cov,
                min_single_cov,
                parsed_pileup,
            )
        if not os.path.exists(simulation_parsed_pileup):
            logging.info(
                "....... start parsing samtools pileup result for simulation data"
            )
            Parse_mpileup_allChr(
                vcf_sample_name,
                vcfgz,
                simulation_pileup,
                min_total_cov,
                min_single_cov,
                simulation_parsed_pileup,
            )
        if os.path.exists(parsed_pileup) and os.path.exists(simulation_parsed_pileup):
            logging.info("================= Skipping common step 1.3")
            logging.info("=================")

    data12 = pd.read_csv(parsed_pileup, sep="\t", header=0, index_col=False)
    logging.debug(
        "output {0} has {1} variants from parsing pileup file".format(
            os.path.basename(parsed_pileup), data12.shape[0]
        )
    )
    if data12.shape[1] < 2:
        os.remove(parsed_pileup)
        logging.error("....... existed parsed pileup file is empty, please try again!")
        sys.exit(1)
    else:
        logging.info(
            "....... {0} save to {1}".format(
                os.path.basename(parsed_pileup), os.path.dirname(parsed_pileup)
            )
        )

    if simulation_pileup is not None:
        data122 = pd.read_csv(
            simulation_parsed_pileup, sep="\t", header=0, index_col=False
        )
        logging.debug(
            "output {0} has {1} variants from parsing pileup file".format(
                os.path.basename(simulation_parsed_pileup), data122.shape[0]
            )
        )
        if data122.shape[1] < 2:
            os.remove(simulation_parsed_pileup)
            logging.error(
                "....... existed parsed pileup file for simulation data is empty, please try again!"
            )
            sys.exit(1)
        else:
            logging.info(
                "....... {0} save to {1}".format(
                    os.path.basename(simulation_parsed_pileup),
                    os.path.dirname(simulation_parsed_pileup),
                )
            )

    #####
    ##### 1.4 Combine hetSNPs and parsed mpileup & thinning reads: one reads only count once
    #####
    hetSNP_intersect_unique = os.path.join(
        tmp_path, f"TEMP.{prefix}_hetSNP_intersected_filtered.tsv"
    )
    logging.info("=================")
    logging.info("================= Starting specific step 1.4")
    if simulation_pileup is not None:
        hetSNP_intersect_unique_sim = os.path.join(
            tmp_path, f"TEMP.{prefix}_hetSNP_intersected_filtered.simulation.tsv"
        )
    else:
        hetSNP_intersect_unique_sim = None

    if not os.path.exists(hetSNP_intersect_unique):
        logging.info(
            "....... start filtering variants with min total counts {0}, and min single allele counts {1}, each read of read length {2} only count for each variant".format(
                min_total_cov, min_single_cov, read_length
            )
        )
        Intersect_exonicHetSnps(
            parsed_pileup,
            simulation_parsed_pileup,
            hetSNP_AF,
            read_length,
            min_total_cov,
            min_single_cov,
            hetSNP_intersect_unique,
            hetSNP_intersect_unique_sim,
        )
    else:
        logging.info("================= Skipping specific step 1.4")
        logging.info("=================")
    data14 = pd.read_csv(hetSNP_intersect_unique, sep="\t", header=0, index_col=False)
    logging.debug(
        "output {0} has {1} het SNPs with AF annotation after taking intersection between output from 1.2 and 1.3".format(
            os.path.basename(hetSNP_intersect_unique), data14.shape[0]
        )
    )
    if data14.shape[1] < 2:
        os.remove(hetSNP_intersect_unique)
        logging.error(
            "....... existed {0} with filtered sites file is empty, please try again!".format(
                os.path.basename(hetSNP_intersect_unique)
            )
        )
        sys.exit(1)
    else:
        logging.info(
            "....... {0} save to {1}".format(
                os.path.basename(hetSNP_intersect_unique),
                os.path.dirname(hetSNP_intersect_unique),
            )
        )
    if hetSNP_intersect_unique_sim is not None:
        data144 = pd.read_csv(
            hetSNP_intersect_unique_sim, sep="\t", header=0, index_col=False
        )
        logging.debug(
            "output {0} has {1} het SNPs with AF annotation after taking intersection between output from 1.2 and 1.3".format(
                os.path.basename(hetSNP_intersect_unique_sim), data144.shape[0]
            )
        )
        if data144.shape[1] < 2:
            os.remove(hetSNP_intersect_unique_sim)
            logging.error(
                "....... existed {0} with filtered sites file is empty, please try again!".format(
                    os.path.basename(hetSNP_intersect_unique_sim)
                )
            )
            sys.exit(1)
        else:
            logging.info(
                "....... {0} save to {1}".format(
                    os.path.basename(hetSNP_intersect_unique_sim),
                    os.path.dirname(hetSNP_intersect_unique_sim),
                )
            )
    logging.info("================= finish step1! ")

    return hetSNP_intersect_unique, hetSNP_intersect_unique_sim
