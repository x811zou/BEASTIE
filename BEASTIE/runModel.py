#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import logging
import os
import sys
import shutil
import multiprocessing
from pathlib import Path
import pandas as pd
from . import binomial_for_real_data
from pkg_resources import resource_filename
from .helpers import runhelper
from .parse_mpileup import Parse_mpileup_allChr
from .annotateLD import annotateLD
from .prepare_model import (
    filter_alignBias,
    re_allocateReads,
    add_simulationData,
    generate_modelCount,
    significant_genes,
    update_model_input_lambda_phasing,
)
from statsmodels.stats.multitest import multipletests
current_script_path = os.path.abspath(__file__)
current_script_dir = os.path.dirname(current_script_path)
quickBEAST_path = os.path.join(current_script_dir, '../QuickBEAST')
sys.path.append(quickBEAST_path)

import calculate_p_value_from_qb_mode


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


def is_valid_parsed_pileup(filepath):
    if not os.path.exists(filepath):
        return False
    parsed_pileup_data = pd.read_csv(filepath, sep="\t", header=0, index_col=False)
    if parsed_pileup_data.shape[0] < 2:
        os.remove(filepath)
        logging.info("....... existed parsed pileup file is empty, removing")
        return False
    return True


def check_file_existence2(
    vcfgz,
    simulation_pileup,
    filtered_hetSNP_intersect_pileup,
    shapeit2_file,
    collected_alignmentBias_file,
):
    ##### VCFGZ file
    vcfgztbi = "{0}.tbi".format(vcfgz)
    if not os.path.isfile(vcfgz):
        logging.error(
            "Oops! filtered vcfgz file {0} doesn't exist. Please try again ...".format(
                vcfgz
            )
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
    ##### simulation pileup_file
    if simulation_pileup is not None:
        if not os.path.isfile(simulation_pileup):
            logging.error(
                "Oops! simulation pileup file {0} doesn't exist. Please try again ...".format(
                    simulation_pileup
                )
            )
        else:
            logging.info(
                "Great! simulation pileup file {0} exists.".format(simulation_pileup)
            )

    ##### real data filtered hetSNP pileup file
    if not os.path.exists(filtered_hetSNP_intersect_pileup):
        logging.error(
            "Oops! hetSNP file {0} doesn't exist. Please try again ...".format(
                filtered_hetSNP_intersect_pileup
            )
        )
        exit(1)
    else:
        hetsnp = pd.read_csv(
            filtered_hetSNP_intersect_pileup, sep="\t", header=0, index_col=False
        )
        logging.debug(
            "output {0} has {1} het SNPs from VCF file".format(
                os.path.basename(filtered_hetSNP_intersect_pileup), hetsnp.shape[0]
            )
        )
        if hetsnp.shape[0] < 2:
            os.remove(filtered_hetSNP_intersect_pileup)
            logging.error("..... existed hetSNP file is empty, please try again!")
            sys.exit(1)

    ##### shapeit2
    if shapeit2_file is not None:
        if not os.path.isfile(shapeit2_file):
            logging.error(
                "Oops!  shapeit2 phased file {0} doesn't exist. Please try again ...".format(
                    shapeit2_file
                )
            )
        else:
            logging.info(
                "Great! shapeit2 phased file {0} exists.".format(shapeit2_file)
            )
    ##### collected_alignmentBias_file
    if collected_alignmentBias_file is not None:
        if not os.path.isfile(collected_alignmentBias_file):
            logging.error(
                "Oops!  File with alignment Bias SNPs {0} doesn't exist. Please try again ...".format(
                    collected_alignmentBias_file
                )
            )
        else:
            logging.info(
                "Great! File with alignment Bias SNPs {0} exists.".format(
                    collected_alignmentBias_file
                )
            )


def create_file_name(prefix, tmp_path, shapeit2_input):
    ##### TEMP output generation: meta_file
    meta_file = os.path.join(
        tmp_path,
        f"{prefix}.meta.tsv",
    )
    meta_error = os.path.join(
        tmp_path,
        f"{prefix}.meta.w_error.tsv",
    )
    if shapeit2_input is not None:
        hetSNP_intersect_unique_forlambda_file = os.path.join(
            tmp_path,
            f"TEMP.{prefix}_hetSNP_intersected_filtered.shapeit2.dropNA.forLambda.tsv",
        )
        hetSNP_intersect_unique_lambdaPredicted_file = os.path.join(
            tmp_path,
            f"TEMP.{prefix}_hetSNP_intersected_filtered.shapeit2.dropNA.lambdaPredicted.tsv",
        )
    else:
        hetSNP_intersect_unique_forlambda_file = os.path.join(
            tmp_path, f"TEMP.{prefix}_hetSNP_intersected_filtered.forLambda.tsv"
        )
        hetSNP_intersect_unique_lambdaPredicted_file = os.path.join(
            tmp_path, f"TEMP.{prefix}_hetSNP_intersected_filtered.lambdaPredicted.tsv"
        )
    logging.info(
        "We will generate intermedate file {0} ...".format(os.path.basename(meta_file))
    )
    logging.info(
        "We will generate intermedate file {0} ...".format(os.path.basename(meta_error))
    )
    logging.info(
        "We will generate intermediate file {0} ...".format(
            os.path.basename(hetSNP_intersect_unique_forlambda_file)
        )
    )
    logging.info(
        "We will generate intermediate file {0} ...".format(
            os.path.basename(hetSNP_intersect_unique_lambdaPredicted_file)
        )
    )
    return (meta_file, meta_error)


def run(
    prefix,
    vcfgz,
    vcf_sample_name,
    collected_alignmentBias_file,
    simulation_pileup,
    filtered_hetSNP_intersect_pileup,
    output_path,
    tmp_path,
    result_path,
    shapeit2_file,
    binomialp_cutoff,
    ase_cutoff,
    ancestry,
    chr_start,
    chr_end,
    min_total_cov,
    min_single_cov,
    read_length,
    SAVE_INT,
    nophasing,
    atacseq,
    LD_token,
    ldlink_cache_dir,
    ldlink_token_db,
):
    #####
    ##### 2.1 Check input file existence
    #####
    logging.info("=================")
    logging.info("================= Starting specific step 2.1")
    logging.info("....... start checking file existence")
    Path(output_path).mkdir(parents=True, exist_ok=True)
    check_file_existence2(
        vcfgz,
        simulation_pileup,
        filtered_hetSNP_intersect_pileup,
        shapeit2_file,
        collected_alignmentBias_file,
    )
    (meta, meta_error) = create_file_name(prefix, tmp_path, shapeit2_file)
    chr_suffix = f"_chr{chr_start}-{chr_end}"
    #####
    ##### 2.2 parse simulation data
    #####
    logging.info("=================")
    logging.info("================= Starting specific step 2.2")
    logging.info("....... start parsing simulation pileup data")
    if simulation_pileup is not None:
        with multiprocessing.Pool(2) as pool:
            handles = []
            simulation_parsed_pileup = os.path.join(
                tmp_path, f"{prefix}_parsed_pileup{chr_suffix}.simulation.tsv"
            )
            handle = pool.apply_async(
                parse_mpileup,
                (
                    simulation_pileup,
                    simulation_parsed_pileup,
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
    else:
        logging.info("....... simulator data is NOT provided, skip step")

    #####
    ##### 2.3 use simulation data to filter variants with alignment bias
    #####
    logging.info("=================")
    logging.info(
        "================= Starting specific step 2.3 aligmment bias filtering using simulation data"
    )
    logging.info("....... start filtering variants with alignment bias")
    biased_variant = None
    if collected_alignmentBias_file is not None:
        logging.info(
            "....... alignment Bias SNP list is provided {0}".format(
                os.path.basename(collected_alignmentBias_file)
            )
        )
    else:
        logging.info("....... alignment Bias SNP list is NOT provided, skip")
        if simulation_pileup is None:
            logging.info("....... simulator data is NOT provided, skip step")
        else:
            logging.info(
                "....... simulator data is provided {0}".format(
                    os.path.basename(simulation_parsed_pileup)
                )
            )
            biased_variant = add_simulationData(simulation_parsed_pileup)
            # running
            alignBiasfiltered_filename = filter_alignBias(
                prefix,
                tmp_path,
                binomialp_cutoff,
                filtered_hetSNP_intersect_pileup,
                biased_variant,
                collected_alignmentBias_file,
            )
            # checking
            data23 = pd.read_csv(
                alignBiasfiltered_filename, sep="\t", header=0, index_col=False
            )
            if data23.shape[0] < 2:
                os.remove(alignBiasfiltered_filename)
                logging.error(
                    "....... existed {0} is empty, please try again!".format(
                        os.path.basename(alignBiasfiltered_filename)
                    )
                )
                sys.exit(1)
            else:
                logging.info(
                    "....... {0} save to {1}".format(
                        os.path.basename(alignBiasfiltered_filename),
                        os.path.dirname(alignBiasfiltered_filename),
                    )
                )

    #####
    ##### 2.4 phase data with shapeit2 or VCF
    #####
    logging.info("=================")
    logging.info("================= Starting specific step 2.4 phase data")
    logging.info(
        "....... start phasing with shapeit2 or VCF, convert data into model input format"
    )
    if atacseq is True or nophasing is True:
        print(
            "....... Phasing not provided: skip using phasing information"
        )
        phasing_method = "nophasing"
        if simulation_pileup is None:
            alignBiasfiltered_filename = filtered_hetSNP_intersect_pileup
        phased_filename = os.path.join(tmp_path, f"{os.path.splitext(os.path.basename(alignBiasfiltered_filename))[0]}.nophasing.tsv")
        phase_difference_filename = None
    else:
        if shapeit2_file is None:
            logging.info(
                "....... Phasing provided: shapeit2 phasing is NOT provided, we use VCF phasing information"
            )
            phasing_method = "VCF"
            phased_filename = (
                f"{os.path.splitext(alignBiasfiltered_filename)[0]}.phasedByVCF.tsv"
            )
            phase_difference_filename = None
        else:
            logging.info(
                "....... Phasing provided: shapeit2 phasing is provided {0}".format(
                    os.path.basename(shapeit2_file)
                )
            )
            phasing_method = "shapeit2"
            phased_filename = f"{os.path.splitext(alignBiasfiltered_filename)[0]}.phasedByshapeit2.tsv"
            phase_difference_filename = f"{os.path.splitext(alignBiasfiltered_filename)[0]}.phasingDifference.tsv"
    
    # running
    phased_clean_filename = f"{os.path.splitext(phased_filename)[0]}.cleaned.tsv"
    re_allocateReads(
        alignBiasfiltered_filename,
        shapeit2_file,
        phasing_method,
        phased_filename,
        phased_clean_filename,
        collected_alignmentBias_file,
        simulation_pileup,
        phase_difference_filename,
    )

    # checking
    data24_1 = pd.read_csv(phased_filename, sep="\t", header=0, index_col=False)
    if data24_1.shape[0] < 2:
        os.remove(phased_filename)
        logging.error(
            "....... existed {0} is empty, please try again!".format(
                os.path.basename(phased_filename)
            )
        )
        sys.exit(1)
    else:
        logging.info(
            "....... {0} save to {1}".format(
                os.path.basename(phased_filename),
                os.path.dirname(phased_filename),
            )
        )

    logging.info("....... start converting data into model input format")
    logging.info(
        "....... start converting {0} for model input".format(
            os.path.basename(phased_clean_filename)
        )
    )
    (
        base_modelin,
        base_modelin_error,
    ) = generate_modelCount(phased_clean_filename,atacseq)

    #####
    ##### 2.5 Annotation LD
    #####
    if atacseq is not True:
        logging.info("=================")
        logging.info("================= Starting specific step 2.5 annotate LD")
        if phasing_method != "nophasing":
            FORCE_ANNOTATE_LD = False
            if FORCE_ANNOTATE_LD or not os.path.isfile(meta):
                logging.info("....... Phasing provided: start annotating LD information")
                logging.debug("input {0} ".format(os.path.basename(phased_clean_filename)))
                annotateLD(
                    ancestry,
                    phased_clean_filename,
                    LD_token,
                    meta,
                    ldlink_cache_dir,
                    ldlink_token_db,
                )
            else:
                logging.info("=================")
                logging.info("================= Skipping specific step 2.5")
            data25 = pd.read_csv(meta, sep="\t", header=0, index_col=False)
            logging.debug(
                "output {0} has {1} het SNPs for LD (d',r2) annotation".format(
                    os.path.basename(meta), data23.shape[0]
                )
            )
            if data25.shape[0] < 2:
                os.remove(meta)
                logging.error(
                    "....... existed {0} is empty, please try again!".format(
                        os.path.basename(meta)
                    )
                )
                sys.exit(1)
            else:
                for files in os.listdir(os.path.dirname(meta)):
                    if "TEMP_chr" in files:
                        logging.info("..... remove created TEMP files: {0}".format(files))
                        os.remove(os.path.dirname(meta) + "/" + files)
                logging.info(
                    "....... {0} save to {1}".format(
                        os.path.basename(meta),
                        os.path.dirname(meta),
                    )
                )
        else:
            logging.info("....... Phasing not provided: skip annotating LD information")

    #####
    ##### 2.6 logistic regression model predict switching phasing error, linear regression model predicts lambda
    #####
    logging.info("=================")
    logging.info("================= Starting step specific 2.6 regression model")

    if not os.path.isfile(meta_error):
        if phasing_method != "nophasing":
            logging.info(
                "....... start predicting switching eror for {0}".format(
                    os.path.basename(meta)
                )
            )
        predict_lambda_phasing_error = resource_filename("BEASTIE", "predict_lambda_phasingError.R")
        beastie_wd = resource_filename("BEASTIE", ".")
        cmd = f"Rscript --vanilla {predict_lambda_phasing_error} {tmp_path} {prefix} {phased_clean_filename} {meta} {meta_error} {beastie_wd} {phasing_method}"
        runhelper(cmd)
        if os.path.isfile(meta_error):
            logging.info(f"....... finish predicting error")

    if phasing_method != "nophasing":
        data26_2 = pd.read_csv(
            meta_error,
            sep="\t",
            header=None,
            index_col=False,
        )
        logging.info(
            "switching error prediction model: input {0}".format(os.path.basename(meta))
        )
        logging.debug(
            "output {0} has {1} het SNPs with switching error prediction".format(
                os.path.basename(meta_error), data26_2.shape[0] - 1
            )
        )

    if phasing_method != "nophasing":
        #####
        ##### 2.7 adding switching error information to model input
        #####
        logging.info("=================")
        logging.info("================= Starting specific step 2.7")
        logging.info(
            "....... start adding model input with predicted switching error information"
        )

        if not os.path.isfile(base_modelin_error):
            logging.info(
                "output {0} is generated for BEASTIE stan model".format(
                    os.path.basename(base_modelin_error)
                )
            )
        else:
            logging.info(
                "....... {0} exists, overwrites".format(
                    os.path.basename(base_modelin_error)
                )
            )
        default_phasing_error = 0.05
        update_model_input_lambda_phasing(
            "pred_error_GIAB", base_modelin, base_modelin_error, meta_error, default_phasing_error
        )
        model_input = base_modelin_error
    else:
        model_input = base_modelin

    #####
    ##### 2.8 run quickbeast
    #####
    logging.info("=================")
    logging.info("================= Starting specific step 2.8 running quickBEAST")

    qb_executable = '/usr/local/bin/QuickBEAST'

    logging.info("....... running quickBEAST")
    input_genes = calculate_p_value_from_qb_mode.parse_gene_input_file(model_input)
    gene_df = calculate_p_value_from_qb_mode.run_qb_parallel(qb_executable, input_genes)

    logging.info("....... running null simulations for p-value calculation")
    null_simulation_df = calculate_p_value_from_qb_mode.run_null_simulations(qb_executable, input_genes, False)
    gene_df = pd.merge(gene_df, null_simulation_df, on="geneID")

    logging.info("....... calculating p-values")
    gene_df['mode_st_p_value'] = gene_df.apply(lambda row: calculate_p_value_from_qb_mode.calculate_st_2sided_pvalue(row['st_df'], row['st_loc'], row['st_scale'], row['qb_mode']), axis=1)
    gene_df.rename(columns={'mode_st_p_value': 'qb_p_value'}, inplace=True)
    gene_df = gene_df[['geneID', 'n_hets', 'total_count', 'qb_mode', 'qb_p_value']]

    outfilename = os.path.join(result_path, f"{prefix}_quickBEAST.tsv")
    gene_df.to_csv(outfilename, sep="\t", header=True, index=False)

    logging.info("....... saved quickBeast to {0}".format(outfilename))


    logging.info("....... running binomial")
    df_binomial = binomial_for_real_data.run(base_modelin,atacseq)
    logging.info("....... done with running binomial")
    binomial_out_path = os.path.join(result_path, f"{prefix}_binomial.tsv")
    df_binomial.to_csv(binomial_out_path, sep="\t", header=True, index=False)
    logging.info("....... saved binomial to {0}".format(binomial_out_path))

    # join gene_df and df_binomial files
    merged_df = pd.merge(gene_df, df_binomial, on="geneID")
    # FDR adjustment using multipletests function on the "qb" column from merged_df to produce "qb_fdr" column
    qb_p = merged_df["qb_p_value"].values
    _, qb_fdr, _, _ = multipletests(qb_p, method='fdr_bh')
    merged_df["qb_fdr"] = qb_fdr
    ns_p = merged_df["NaiveSum_pval"].values
    _, ns_fdr, _, _ = multipletests(ns_p, method='fdr_bh')
    merged_df["NS_fdr"] = ns_fdr
    ms_p = merged_df["MajorSite_pval"].values
    _, ms_fdr, _, _ = multipletests(ms_p, method='fdr_bh')
    merged_df["MS_fdr"] = ms_fdr

    columns_to_drop = ['NaiveSum_esti','MajorSite_esti']
    merged_df = merged_df.drop(columns=columns_to_drop)
    merged_out_path = os.path.join(result_path, f"{prefix}_ASE.tsv")
    merged_df.to_csv(merged_out_path, sep="\t", header=True, index=False)

    # report number and percentage significant ASE genes found in each method
    significant_qb = merged_df[merged_df["qb_fdr"] <= ase_cutoff]
    significant_ns = merged_df[merged_df["NS_fdr"] <= ase_cutoff]
    significant_ms = merged_df[merged_df["MS_fdr"] <= ase_cutoff]
    logging.info(f"Alpha (ASE cutoff): {ase_cutoff}")
    logging.info(f"Number/Percentage of significant genes found by QuickBEAST: {significant_qb.shape[0]} / {significant_qb.shape[0] / merged_df.shape[0]*100}%")
    logging.info(f"Number of significant genes found by NaiveSum: {significant_ns.shape[0]} / {significant_ns.shape[0] / merged_df.shape[0]*100}%")
    logging.info(f"Number of significant genes found by MajorSite: {significant_ms.shape[0]} / {significant_ms.shape[0] / merged_df.shape[0]*100}%")

    if not SAVE_INT:
        logging.info("....... removing TEMP folder {0}".format(tmp_path))
        shutil.rmtree(tmp_path)

    logging.info("=================")
    logging.info(">>  Yep! You are done running BEASTIE!")
