#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import logging
import os
import sys
import shutil
import pandas as pd
import BEASTIE.ADM_for_real_data as ADM_for_real_data
import BEASTIE.binomial_for_real_data as binomial_for_real_data
import BEASTIE.run_model_stan_wrapper as run_model_stan_wrapper
from pkg_resources import resource_filename
from .helpers import runhelper
from .annotateLD import annotateLD
from .prepare_model import (
    re_allocateReads,
    add_simulationData,
    generate_modelCount,
    significant_genes,
    update_model_input_lambda_phasing,
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
    shapeit2_input,
    hetSNP_intersect_unique,
    hetSNP_intersect_unique_sim,
    prefix,
    alpha,
    model,
    sigma,
    tmp_path,
    result_path,
    cutoff,
    SAVE_INT,
    WARMUP,
    KEEPER,
    total_cov,
    either_cov,
    ancestry,
    LD_token,
    chr_start,
    chr_end,
    ldlink_cache_dir,
    ldlink_token_db,
):
    (meta, meta_error) = create_file_name(prefix, tmp_path, shapeit2_input)

    #####
    ##### 2.1 phase data with shapeit2 or VCF, use simulation data to filter biased variants, and convert data into model input format
    #####
    logging.info("=================")
    logging.info("================= Starting specific step 2.1")
    logging.info(
        "....... start converting {0} for model input".format(
            os.path.basename(hetSNP_intersect_unique)
        )
    )

    p_cutoff = 0.05
    if shapeit2_input is None:
        logging.info(
            "....... shapeit2 phasing is NOT provided, we use VCF phasing information"
        )
        phasing_method = "VCF"
        phased_filename = (
            f"{os.path.splitext(hetSNP_intersect_unique)[0]}.phasedByVCF.tsv"
        )
    else:
        logging.info(
            "....... shapeit2 phasing is provided {0}".format(
                os.path.basename(shapeit2_input)
            )
        )
        phasing_method = "shapeit2"
        phased_filename = (
            f"{os.path.splitext(hetSNP_intersect_unique)[0]}.phasedByshapeit2.tsv"
        )
    filename_cleaned = f"{os.path.splitext(phased_filename)[0]}.cleaned.tsv"
    re_allocateReads(
        shapeit2_input,
        hetSNP_intersect_unique,
        phasing_method,
        phased_filename,
        filename_cleaned,
    )
    data21 = pd.read_csv(phased_filename, sep="\t", header=0, index_col=False)
    if data21.shape[0] < 2:
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

    #####
    ##### 2.2 use simulation data to filter biased variants, and convert data into model input format
    #####

    if hetSNP_intersect_unique_sim is None:
        logging.info("....... simulator data is NOT provided")
        biased_variant = None
    else:
        logging.info(
            f"....... simulator data is provided {hetSNP_intersect_unique_sim}"
        )
        biased_variant = add_simulationData(hetSNP_intersect_unique_sim)

    (
        file_for_LDannotation,
        file_for_lambda,
        lambdaPredicted_file,
        base_modelin,
        base_modelin_error,
    ) = generate_modelCount(filename_cleaned, biased_variant)

    data22 = pd.read_csv(file_for_lambda, sep="\t", header=0, index_col=False)
    logging.debug(
        "output {0} has {1} genes for lambda prediction".format(
            os.path.basename(file_for_lambda), data22.shape[0]
        )
    )
    if data22.shape[0] < 2:
        os.remove(file_for_lambda)
        logging.error(
            "....... existed {0} with filtered sites prepared for lambda model file is empty, please try again!".format(
                os.path.basename(file_for_lambda)
            )
        )
        sys.exit(1)
    else:
        logging.info(
            "....... {0} save to {1}".format(
                os.path.basename(file_for_lambda),
                os.path.dirname(file_for_lambda),
            )
        )

    #####
    ##### 2.3 Annotation LD
    #####
    FORCE_ANNOTATE_LD = False
    if FORCE_ANNOTATE_LD or not os.path.isfile(meta):
        logging.info("=================")
        logging.info("================= Starting specific step 2.3")
        logging.info("....... start annotating LD information")
        logging.debug("input {0} ".format(os.path.basename(file_for_LDannotation)))
        annotateLD(
            ancestry,
            file_for_LDannotation,
            LD_token,
            meta,
            ldlink_cache_dir,
            ldlink_token_db,
        )
    else:
        logging.info("=================")
        logging.info("================= Skipping specific step 2.3")
    data23 = pd.read_csv(meta, sep="\t", header=0, index_col=False)
    logging.debug(
        "output {0} has {1} het SNPs for LD (d',r2) annotation".format(
            os.path.basename(meta), data22.shape[0]
        )
    )
    if data23.shape[0] < 2:
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

    #####
    ##### 2.4 logistic regression model predict switching phasing error, linear regression model predicts lambda
    #####
    logging.info("=================")
    logging.info("================= Starting step specific 2.4")
    logging.info(
        "....... start predicting lambda for {0}".format(
            os.path.basename(file_for_lambda)
        )
    )
    logging.info(
        "....... start predicting switching eror for {0}".format(os.path.basename(meta))
    )

    predict_lambda_phasing_error = resource_filename(
        "BEASTIE", "predict_lambda_phasingError.R"
    )
    beastie_wd = resource_filename("BEASTIE", ".")
    cmd = f"Rscript --vanilla {predict_lambda_phasing_error} {alpha} {tmp_path} {prefix} {model} {file_for_LDannotation} {file_for_lambda} {lambdaPredicted_file} {meta} {meta_error} {beastie_wd}"
    runhelper(cmd)
    data24_1 = pd.read_csv(
        lambdaPredicted_file,
        sep="\t",
        header=None,
        index_col=False,
    )
    data24_2 = pd.read_csv(
        meta_error,
        sep="\t",
        header=None,
        index_col=False,
    )
    logging.info(
        "lambda prediction model: input {0}".format(os.path.basename(file_for_lambda))
    )
    logging.info(
        "lambda prediction model: input alpha (family wise error rate) is {0}, adjusted after size of input {1} is {2}".format(
            alpha, data24_1.shape[0], alpha / data24_1.shape[0]
        )
    )
    logging.debug(
        "output {0} has {1} genes with lambda prediction".format(
            os.path.basename(lambdaPredicted_file),
            data24_1.shape[0],
        )
    )
    logging.info(
        "switching error prediction model: input {0}".format(os.path.basename(meta))
    )
    logging.debug(
        "output {0} has {1} het SNPs with switching error prediction".format(
            os.path.basename(meta_error), data24_2.shape[0] - 1
        )
    )

    #####
    ##### 2.5 adding switching error information to model input
    #####
    logging.info("=================")
    logging.info("================= Starting specific step 2.5")
    logging.info(
        "....... start adding model input with predicted swhitching error information"
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
    update_model_input_lambda_phasing(
        "pred_error_GIAB", base_modelin, base_modelin_error, meta_error
    )

    #####
    ##### 2.6 running stan model
    #####
    logging.info("=================")
    logging.info("================= Starting specific step 2.5")
    logging.info("....... start running iBEASTIE model")
    df_ibeastie, picklename = run_model_stan_wrapper.run(
        prefix,
        base_modelin_error,
        sigma,
        alpha,
        model,
        result_path,
        lambdaPredicted_file,
        WARMUP,
        KEEPER,
        either_cov,
        total_cov,
    )
    if df_ibeastie.shape[0] > 2:
        logging.info("....... done with running model!")
    else:
        logging.error("....... model output is empty, please try again!")
        sys.exit(1)

    df_adm = ADM_for_real_data.run(prefix, base_modelin, result_path, picklename)
    logging.info("....... done with running ADM method")

    df_binomial = binomial_for_real_data.run(
        prefix, base_modelin, result_path, picklename
    )
    logging.info("....... done with running binomial")

    #####
    ##### 2.7 generating output
    #####
    logging.info("=================")
    logging.info("================= Starting specific step 2.6")
    logging.info("....... start generating gene list")
    outfilename = os.path.join(result_path, f"{prefix}_ASE_all.tsv")
    outfilename_ase = os.path.join(
        result_path, f"{prefix}_ASE_cutoff_{cutoff}_filtered.tsv"
    )
    if (os.path.isfile(outfilename)) and (os.path.isfile(outfilename_ase)):
        logging.info(
            "....... {0} exists, but we will overwrites".format(
                os.path.basename(outfilename)
            )
        )
        logging.info(
            "....... {0} exists, but we will overwrites".format(
                os.path.basename(outfilename_ase)
            )
        )
    significant_genes(
        df_ibeastie,
        df_binomial,
        df_adm,
        outfilename,
        outfilename_ase,
        cutoff,
        lambdaPredicted_file,
    )
    logging.info("....... done with significant_gene")
    if not SAVE_INT:
        shutil.rmtree(tmp_path)
        logging.info("....... remove TEMP folder {0}".format(tmp_path))

    data26_1 = pd.read_csv(outfilename, sep="\t", header=0, index_col=False)
    data26_2 = pd.read_csv(outfilename_ase, sep="\t", header=0, index_col=False)
    if data26_1.shape[0] <= 2 or data26_2.shape[0] <= 2:
        logging.error(
            "....... existed {0} or {1} is empty, please try again!".format(
                os.path.basename(outfilename),
                os.path.basename(outfilename_ase),
            )
        )
        sys.exit(1)
    logging.info(
        "output {0} has all {1} genes".format(
            os.path.basename(outfilename),
            data26_1.shape[0],
        )
    )
    logging.info(
        "output {0} has filtered {1} genes passing ASE cutoff {2}".format(
            os.path.basename(outfilename_ase),
            data26_2.shape[0],
            cutoff,
        )
    )
    logging.info("=================")
    logging.info(">>  Yep! You are done running BEASTIE!")
