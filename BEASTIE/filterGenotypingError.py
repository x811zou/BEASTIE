#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
from faulthandler import is_enabled
import dataclasses
import multiprocessing
import os
import re
import sys
import logging
from unittest.mock import NonCallableMagicMock
import pandas as pd
import numpy as np
from pathlib import Path
from pkg_resources import resource_filename
from .annotationAF import annotateAF
from .helpers import runhelper
from .intersect_hets import Intersect_exonicHetSnps
from .parse_mpileup import Parse_mpileup_allChr
from scipy.stats import fisher_exact


def parse_mpileup(
    input_file, output_file, vcf_sample_name, vcfgz, min_total_cov, min_single_cov
):
    if not is_valid_parsed_pileup(output_file):
        logging.info(
            "....... start parsing samtools pileup result for {0}".format(output_file)
        )
        Parse_mpileup_allChr(
            vcf_sample_name,
            vcfgz,
            input_file,
            min_total_cov,
            min_single_cov,
            output_file,
        )
        parsed_pileup_data = pd.read_csv(
            output_file, sep="\t", header=0, index_col=False
        )
        logging.debug(
            "output {0} has {1} variants from parsing simulation mpileup file".format(
                output_file, parsed_pileup_data.shape[0]
            )
        )
    else:
        logging.info(f"...... skipped parsing mpileup for {input_file}")
    return "hello"


def check_file_existence(vcfgz, pileup, het_snp_file):
    ##### VCFGZ file
    vcfgztbi = "{0}.tbi".format(vcfgz)
    if not os.path.isfile(vcfgz):
        logging.error(
            "Oops! vcfgz file {0} doesn't exist. Please try again ...".format(vcfgz)
        )
        exit(1)
    if os.path.isfile(vcfgz) and (
        not os.path.isfile(vcfgztbi)
        or os.path.getmtime(vcfgztbi) < os.path.getmtime(vcfgz)
    ):
        logging.warning("We will generate the latest vcfgz index for you ...")
        cmd = "tabix -fp vcf %s" % (vcfgz)
        runhelper(cmd)
    if os.path.isfile(vcfgz) and os.path.isfile(vcfgztbi):
        logging.info(
            "Great! VCFGZ file {0} and index {1} exists.".format(vcfgz, vcfgztbi)
        )
    ##### pileup_file
    if not os.path.isfile(pileup):
        logging.error(
            "Oops!  pileup file {0} doesn't exist. Please try again ...".format(pileup)
        )
    else:
        logging.info("Great! pileup file {0} exists.".format(pileup))

    ##### hetSNP file
    if not os.path.exists(het_snp_file):
        logging.error(
            "Oops! hetSNP file {0} doesn't exist. Please try again ...".format(
                het_snp_file
            )
        )
        exit(1)
    else:
        hetsnp = pd.read_csv(het_snp_file, sep="\t", header=0, index_col=False)
        logging.debug(
            "output {0} has {1} het SNPs from VCF file".format(
                os.path.basename(het_snp_file), hetsnp.shape[0]
            )
        )
        if hetsnp.shape[0] < 2:
            os.remove(het_snp_file)
            logging.error("..... existed hetSNP file is empty, please try again!")
            sys.exit(1)


def is_valid_parsed_pileup(filepath):
    if not os.path.exists(filepath):
        return False
    parsed_pileup_data = pd.read_csv(filepath, sep="\t", header=0, index_col=False)
    if parsed_pileup_data.shape[0] < 2:
        os.remove(filepath)
        logging.info("....... existed parsed pileup file is empty, removing")
        return False
    return True


def process_gene(data):
    data_sub = data[["chrN", "pos", "refCount", "altCount"]]
    data_rest = data_sub[(data_sub["refCount"] != 0) & (data_sub["altCount"] != 0)]
    size = data_rest.shape[0]
    if size == 0:
        min_rest_data = "NA"
        max_rest_data = "NA"
    else:
        data_rest["max_count"] = data_rest[["refCount", "altCount"]].max(axis=1)
        data_rest["min_count"] = data_rest[["refCount", "altCount"]].min(axis=1)
        min_rest_data = data_rest["min_count"].sum()
        max_rest_data = data_rest["max_count"].sum()

    def process_row(row):
        fishTest = 100
        if ((row["refCount"] == 0) or (row["altCount"] == 0)) and (size > 0):
            zero_counts = row[["refCount", "altCount"]].min()
            nonzero_counts = row[["refCount", "altCount"]].max()
            table = np.array(
                [[nonzero_counts, zero_counts], [max_rest_data, min_rest_data]]
            )
            oddsr, p = fisher_exact(table, alternative="two-sided")
            fishTest = p
        return (
            row["chrN"],
            row["pos"],
            row["refCount"],
            row["altCount"],
            max_rest_data,
            min_rest_data,
            fishTest,
        )

    return data.apply(process_row, axis=1)


def filter_genotypeEr(
    genotypeEr_cutoff,
    hetSNP_intersect_unique_filename,
    filtered_hetSNP_intersec_pileup,
    genotypeEr_filename,
):

    hetSNP_intersect_unique = pd.read_csv(
        hetSNP_intersect_unique_filename, sep="\t", header=0, index_col=False
    )
    hetSNP_intersect_unique = hetSNP_intersect_unique[
        [
            "chr",
            "chrN",
            "pos",
            "rsid",
            "AF",
            "geneID",
            "genotype",
            "refCount",
            "altCount",
            "totalCount",
            "altRatio",
        ]
    ]
    base_out = os.path.splitext(hetSNP_intersect_unique_filename)[0]
    applyFilter_filename = f"{base_out}_underGenotypingErTesting.tsv"
    beforeFilter_filename = f"{base_out}_beforeGenotypingErFiltered.tsv"
    # genotyping error fisher exact test score
    grouped_df = (
        hetSNP_intersect_unique.groupby("geneID").apply(process_gene).reset_index()
    )
    grouped_df[
        [
            "chrN",
            "pos",
            "refCount",
            "altCount",
            "max_rest_data",
            "min_reast_data",
            "FishTest",
        ]
    ] = pd.DataFrame(grouped_df[0].tolist(), index=hetSNP_intersect_unique.index)
    grouped_df_sub = grouped_df.drop(grouped_df.columns[[1, 2]], axis=1)
    grouped_df_sub.to_csv(applyFilter_filename, index=False, sep="\t", header=True)
    new_df = pd.merge(
        hetSNP_intersect_unique,
        grouped_df_sub,
        how="inner",
        on=["chrN", "pos", "geneID", "refCount", "altCount"],
    )
    new_df.to_csv(beforeFilter_filename, index=False, sep="\t", header=True)
    debiased_df = new_df[new_df["FishTest"] > genotypeEr_cutoff]
    logging.debug(
        f"{debiased_df.shape[0]} out of {new_df.shape[0]} ({round(debiased_df.shape[0] / new_df.shape[0] * 100,2,)}%) het SNPs pass genotyping error with Fisher exact test p-val > {genotypeEr_cutoff}"
    )
    debiased_df.to_csv(
        filtered_hetSNP_intersec_pileup, index=False, sep="\t", header=True
    )
    biased_df = new_df[new_df["FishTest"] <= genotypeEr_cutoff]
    biased_df = biased_df[["chrN", "pos"]]
    biased_df.to_csv(genotypeEr_filename, index=False, sep="\t", header=True)
    logging.debug(
        f"{biased_df.shape[0]} het SNPs with genotyping error will be filtered out in shapeit2 phasing and simulation steps"
    )


def run(
    prefix,
    vcf_sample_name,
    output_path,
    af_path,
    vcfgz,
    ancestry,
    chr_start,
    chr_end,
    read_length,
    min_single_cov,
    min_total_cov,
    het_snp_file,
    genotypeErfiltered_file,
    pileup,
    genotypeEr_cutoff,
    filtered_hetSNP_filename,
):
    #####
    ##### 1.1 Check input file existence
    #####
    tmp_path = output_path + "/" + "filterGenotypingError"
    Path(output_path).mkdir(parents=True, exist_ok=True)
    Path(tmp_path).mkdir(parents=True, exist_ok=True)
    check_file_existence(vcfgz, pileup, het_snp_file)
    chr_suffix = f"_chr{chr_start}-{chr_end}"

    #####
    ##### 1.2 Annotation: AF
    #####
    logging.info("=================")
    logging.info("================= Starting common step 1.2")
    tmp_path = output_path + "/" + "filterGenotypingError"
    hetSNP_AF = os.path.join(tmp_path, f"{prefix}_hetSNP_AF{chr_suffix}.tsv")
    if os.path.isfile(hetSNP_AF):
        logging.info("================= Skipping common step 1.2")
        logging.info("=================")
    else:
        logging.info("..... start annotating hetSNP from 1.1 with AF from 1000 Genome")
        annotateAF(af_path, ancestry, het_snp_file, hetSNP_AF)

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
    logging.info("=================")
    logging.info("================= Starting common step 1.3")
    # required data

    with multiprocessing.Pool(2) as pool:
        handles = []
        parsed_pileup = os.path.join(
            tmp_path, f"{prefix}_parsed_pileup{chr_suffix}.tsv"
        )
        handle = pool.apply_async(
            parse_mpileup,
            (
                pileup,
                parsed_pileup,
                vcf_sample_name,
                vcfgz,
                min_total_cov,
                min_single_cov,
            ),
        )
        handles.append(handle)

        pool.close()
        pool.join()

        for handle in handles:
            handle.get()

    #####
    ##### 1.4 Combine hetSNPs and parsed mpileup & thinning reads: one reads only count once
    #####
    hetSNP_intersect_pileup = os.path.join(
        tmp_path, f"TEMP.{prefix}_hetSNP_intersected_filtered.tsv"
    )
    logging.info("=================")
    logging.info("================= Starting specific step 1.4")

    if not os.path.exists(hetSNP_intersect_pileup):
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
            hetSNP_intersect_pileup,
        )
    else:
        logging.info("================= Skipping specific step 1.4")
        logging.info("=================")
    data14 = pd.read_csv(hetSNP_intersect_pileup, sep="\t", header=0, index_col=False)

    logging.debug(
        "output {0} has {1} het SNPs with AF annotation after taking intersection between output from 1.2 and 1.3".format(
            os.path.basename(hetSNP_intersect_pileup), data14.shape[0]
        )
    )
    if data14.shape[1] < 2:
        os.remove(hetSNP_intersect_pileup)
        logging.error(
            "....... existed {0} with filtered sites file is empty, please try again!".format(
                os.path.basename(hetSNP_intersect_pileup)
            )
        )
        sys.exit(1)
    else:
        logging.info(
            "....... {0} save to {1}".format(
                os.path.basename(hetSNP_intersect_pileup),
                os.path.dirname(hetSNP_intersect_pileup),
            )
        )

    #####
    ##### 1.5 Use fisher exact test to filter variants with genotyping error
    #####

    logging.info("=================")
    logging.info("================= Starting specific step 1.5")
    logging.info("....... start filtering variants with genotyping error")

    filter_genotypeEr(
        genotypeEr_cutoff,
        hetSNP_intersect_pileup,
        filtered_hetSNP_filename,
        genotypeErfiltered_file,
    )

    # checking
    data15 = pd.read_csv(genotypeErfiltered_file, sep="\t", header=0, index_col=False)
    if data15.shape[0] < 2:
        os.remove(genotypeErfiltered_file)
        logging.error(
            "....... existed {0} is empty, please try again!".format(
                os.path.basename(genotypeErfiltered_file)
            )
        )
        sys.exit(1)
    else:
        logging.info(
            "....... {0} save to {1}".format(
                os.path.basename(genotypeErfiltered_file),
                output_path,
            )
        )
