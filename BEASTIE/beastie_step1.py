#!/usr/bin/env python
# =========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
# =========================================================================
import os
import sys
import logging
import pandas as pd
from pathlib import Path

from pkg_resources import resource_filename
from . import annotation
from .extractHets import count_all_het_sites
from .helpers import runhelper
from .intersect_hets import Intersect_exonicHetSnps
from .parse_mpileup import Parse_mpileup_allChr


def check_file_existence(
    prefix,
    in_path,
    output_path,
    tmp_path,
    model,
    vcf,
    ref_dir,
    pileup,
    hetSNP,
    parsed_pileup,
    chr_start,
    chr_end,
):
    ##### STAN model
    if not os.path.exists(model):
        logging.error(
            "Oops! STAN model {0} doesn't exist. Please try again ...".format(model)
        )
        sys.exit(1)
    ##### vcf & vcfgz
    vcfgz = "{0}.gz".format(vcf)
    vcfgztbi = "{0}.tbi".format(vcfgz)
    if not os.path.isfile(vcf) and not os.path.isfile(vcfgz):
        logging.error(
            "Oops! vcf file {0} or vcf.gz file {1} doesn't exist. Please try again ...".format(
                vcf, vcfgz
            )
        )
        exit(1)
    elif os.path.isfile(vcf) and (not os.path.isfile(vcfgz)):
        logging.warning("We will generate vcfgz for you ...".format(vcfgz))
        cmd = "bgzip -c %s > %s" % (vcf, vcfgz)
        runhelper(cmd)
        cmd = "tabix -fp vcf %s" % (vcfgz)
        runhelper(cmd)
    elif (os.path.isfile(vcfgz)) and (not os.path.isfile(vcf)):
        logging.warning(
            "Oops! VCF file {0} not found. We will generate that for you ...".format(
                vcf
            )
        )
        cmd = "gunzip -k %s" % (vcfgz)
        runhelper(cmd)
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
    ##### pileup_file
    pileup_file = os.path.join(in_path, pileup)
    if not os.path.isfile(pileup_file):
        logging.error(
            "Oops!  pileup file {0} doesn't exist in {1}. Please try again ...".format(
                pileup, in_path
            )
        )
    ##### hetSNP_file
    if hetSNP is not None:
        hetSNP_file = in_path + hetSNP
        if not os.path.isfile(hetSNP_file):
            logging.warning(
                "Alright, hetSNP file {0} doesn't exist in {1}. We will generate that for you ...".format(
                    hetSNP, in_path
                )
            )
        else:
            logging.warning(
                "Found existed hetSNP file {0} doesn't exist in {1}.".format(
                    hetSNP, in_path
                )
            )
    else:
        hetSNP_file = os.path.join(
            output_path, f"{prefix}_hetSNP_chr{chr_start}-{chr_end}.tsv"
        )
        logging.info("We will generate {0} for you ...".format(hetSNP_file))
    ##### TEMP output generation: hetSNP_file
    hetSNP_intersect_unique_file = os.path.join(
        tmp_path, f"TEMP.{prefix}_hetSNP_intersected_filtered.tsv"
    )

    logging.info(
        "We will generate intermediate file {0} ...".format(
            hetSNP_intersect_unique_file
        )
    )

    ##### parsed_pileup_file
    if parsed_pileup is not None:
        parsed_pileup_file = os.path.join(in_path, parsed_pileup)
        if not os.path.isfile(parsed_pileup_file):
            logging.warning(
                "Alright, parsed_pileup file {0} doesn't exist in {1}. We will generate that for you ...".format(
                    parsed_pileup, in_path
                )
            )
    else:
        parsed_pileup_file = os.path.join(
            output_path, f"{prefix}_parsed_pileup_chr{chr_start}-{chr_end}.tsv"
        )
        logging.info("We will generate {0} for you ...".format(parsed_pileup_file))
    return (
        vcfgz,
        pileup_file,
        hetSNP_file,
        hetSNP_intersect_unique_file,
        parsed_pileup_file,
    )


def run(
    prefix,
    vcf_sample_name,
    in_path,
    output_path,
    tmp_path,
    gencode_file,
    model,
    vcf,
    ref_dir,
    ancestry,
    chr_start,
    chr_end,
    min_total_cov,
    min_single_cov,
    read_length,
    pileup,
    hetSNP,
    parsed_pileup,
):
    (
        vcfgz,
        pileup,
        hetSNP,
        hetSNP_intersect_unique,
        parsed_pileup,
    ) = check_file_existence(
        prefix,
        in_path,
        output_path,
        tmp_path,
        model,
        vcf,
        ref_dir,
        pileup,
        hetSNP,
        parsed_pileup,
        chr_start,
        chr_end,
    )
    #####
    ##### 0.0 simulate data from BAM file
    #####

    # if not os.path.exists(simulated_folder) and os.path.isfile(fastq1) and os.path.isfile(fastq2):
    #     logging.info("=================")
    #     logging.info("================= Starting common step 0.0")
    #     logging.info("....... start simulating reads")
    #     simulate_reads(simulated_folder,fastq1,fastq2,depth,config)
    # else:
    #     logging.info("=================")
    #     logging.info("================= Skipping common step 0.0")
    # samtools view -h -o $SAM $BAM

    #####
    ##### 1.1 Generate hetSNP file: extract heterozygous bi-allelic SNPs for specific chromosomes from all gencode transcripts
    #####

    if not os.path.exists(hetSNP):
        logging.info("=================")
        logging.info("================= Starting common step 1.1")
        logging.info("....... start extracting heterozygous bi-allelic SNPs from VCF")
        count_all_het_sites(
            prefix,
            vcfgz,
            hetSNP,
            int(chr_start),
            int(chr_end),
            gencode_file,
        )
    else:
        logging.info("=================")
        logging.info("================= Skipping common step 1.1")
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
    ##### 1.2 Generate parsed pileup file: parse samtools mpile up output files
    #####
    if not os.path.exists(parsed_pileup):
        logging.info("=================")
        logging.info("================= Starting common step 1.2")
        logging.info("....... start parsing samtools pileup result")
        Parse_mpileup_allChr(
            vcf_sample_name, vcfgz, pileup, min_total_cov, min_single_cov, parsed_pileup
        )
    else:
        logging.info("=================")
        logging.info("================= Skipping common step 1.2")
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

    #####
    ##### 1.3 Annotation: AF
    #####
    hetSNP_AF = f"{os.path.splitext(hetSNP)[0]}_AF.tsv"
    if os.path.isfile(hetSNP_AF):
        logging.info("=================")
        logging.info("================= Skipping common step 1.3")

    else:
        logging.info("=================")
        logging.info("================= Starting common step 1.3")
        logging.info("..... start annotating hetSNP from 1.1 with AF from 1000 Genome")
        annotation.annotateAF(ancestry, hetSNP, hetSNP_AF)

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
    ##### 1.4 Combine hetSNPs and parsed mpileup & thinning reads: one reads only count once
    #####
    if not os.path.exists(hetSNP_intersect_unique):
        logging.info("=================")
        logging.info("================= Starting specific step 1.4")
        logging.info(
            "....... start filtering variants with min total counts {0}, and min single allele counts {1}, each read of read length {2} only count for each variant".format(
                min_total_cov, min_single_cov, read_length
            )
        )
        Intersect_exonicHetSnps(
            parsed_pileup,
            hetSNP_AF,
            read_length,
            min_total_cov,
            min_single_cov,
            hetSNP_intersect_unique,
        )
    else:
        logging.info("=================")
        logging.info("================= Skipping specific step 1.4")
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
    logging.info("================= finish step1! ")
    sys.exit()
    return hetSNP_intersect_unique
