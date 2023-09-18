#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import multiprocessing
import os
import pickle
import logging
import statistics
import subprocess
import tempfile
import time
import numpy as np
import pandas as pd
from math import floor, log10, log2
from .misc_tools.StanParser import StanParser
from .helpers import runhelper
from scipy import stats
import statistics
import sys
import math


def writeInitializationFile(filename):
    OUT = open(filename, "wt")
    print("theta <- 1", file=OUT)
    OUT.close()


def writeOutputsFile(filename):
    OUT = open(filename, "wt")
    OUT.close()


def writeReadCounts(fields, start, numReps, varName, OUT):
    print(varName, "<- c(", file=OUT, end="")
    for rep in range(numReps):
        print(fields[start + rep * 2], file=OUT, end="")
        if rep + 1 < numReps:
            print(",", file=OUT, end="")
    print(")", file=OUT)


def writePi(fields, numReps, varName, OUT):
    print(varName, "<- c(", file=OUT, end="")
    start = numReps * 2 + 3
    for rep in range(start, len(fields)):
        print(fields[rep], file=OUT, end="")
        if rep + 1 < len(fields):
            print(",", file=OUT, end="")
    print(")", file=OUT)


def writeInputsFile(fields, filename, sigma):
    Mreps = int(fields[1])
    OUT = open(filename, "wt")
    print("M <-", Mreps, file=OUT)
    writeReadCounts(fields, 2, Mreps, "A", OUT)  # alt
    writeReadCounts(fields, 3, Mreps, "R", OUT)  # ref
    print("sigma <-", sigma, file=OUT)
    OUT.close()


def writeInputsFile_i(fields, filename, sigma):
    Mreps = int(fields[1])
    OUT = open(filename, "wt")
    print("M <-", Mreps, file=OUT)
    writeReadCounts(fields, 2, Mreps, "A", OUT)  # alt
    writeReadCounts(fields, 3, Mreps, "R", OUT)  # ref
    print("sigma <-", sigma, file=OUT)
    writePi(fields, Mreps, "pi", OUT)  # ref
    print("N_MISSING_PI <-", fields[Mreps * 2 + 2], file=OUT)
    OUT.close()


def getBaseline(fields):
    if len(fields) >= 5:
        base_thetas = []
        Mreps = int(fields[1])
        for rep in range(Mreps):
            A = float(fields[2 + rep * 2])
            R = float(fields[3 + (rep) * 2])
            base = (A + 1) / (R + 1)
            base_thetas.append(base)
        med_base_theta = statistics.median(base_thetas)
    return med_base_theta


def getFieldIndex(label, fields):
    numFields = len(fields)
    index = None
    for i in range(7, numFields):
        if fields[i] == label:
            index = i
    return index


def computeBeastieScoreLog2(log2_thetas, l):
    assert l >= 1

    # 1. no transformation
    min_l = 1 / l
    max_l = l
    min_l_log2 = math.log2(min_l)
    max_l_log2 = math.log2(max_l)

    n_less_log2 = np.count_nonzero(log2_thetas < min_l_log2)
    n_more_log2 = np.count_nonzero(log2_thetas > max_l_log2)

    n_total = len(log2_thetas)
    max_log2_score = max(n_less_log2, n_more_log2) / n_total
    sum_log2_score = (n_less_log2 + n_more_log2) / n_total
    return max_log2_score, sum_log2_score


def runModel(
    model,
    fields,
    tmp_output_file,
    stan_output_file,
    init_file,
    sigma,
    WARMUP,
    KEEPER,
    phasing_method,
):
    if len(fields) >= 4:
        geneID = str(fields[0])
        writeOutputsFile(stan_output_file)
        # logging.debug(geneID)
        if phasing_method == "nophasing":
            writeInputsFile(fields, tmp_output_file, sigma)
        else:
            writeInputsFile_i(fields, tmp_output_file, sigma)
        writeInitializationFile(init_file)
        cmd = (
            "%s sample num_samples=%s num_warmup=%s data file=%s init=%s output file=%s refresh=0"
            % (model, KEEPER, WARMUP, tmp_output_file, init_file, stan_output_file)
        )
        # print(cmd)
        # logging.debug(cmd)
        try:
            subprocess.run(
                cmd,
                shell=True,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except Exception as e:
            logging.error("STAN MODEL FAILED")
            logging.error(fields)
            raise e
        parser = StanParser(stan_output_file)
        thetas = parser.getVariable("theta")
        return geneID, thetas
    else:
        logging.error("lines with no enough elements")


def getCredibleInterval(thetas, alpha, n):
    halfAlpha = alpha / 2.0
    leftIndex = int(halfAlpha * n)
    rightIndex = n - leftIndex
    left = thetas[leftIndex + 1]
    right = thetas[rightIndex - 1]
    return left, right


def summarize(thetas, alpha):
    mean = statistics.mean(thetas)
    median = statistics.median(thetas)
    thetas_log2 = [log2(x) for x in thetas]
    thetas_log2_abs = [abs(log2(x)) for x in thetas]
    log2_mean = statistics.mean(thetas_log2)
    log2_median = statistics.median(thetas_log2)
    abslog2_mean = statistics.mean(thetas_log2_abs)
    abslog2_median = statistics.median(thetas_log2_abs)
    # print(f"median {median}")
    variance = np.var(thetas)
    log2_variance = np.var(thetas_log2)
    abslog2_variance = np.var(thetas_log2_abs)
    # print(f"variance {variance}")
    # print(f"length of thetas {len(thetas)}")
    n = len(thetas)
    thetas.sort()
    # print(f"length of sorted thetas {len(thetas)}")
    CI_left, CI_right = getCredibleInterval(thetas, alpha, n)
    # mad = stats.median_abs_deviation(thetas, scale="normal")
    # mad = stats.median_abs_deviation(thetas, scale=1/1.4826)
    mad = stats.median_absolute_deviation(thetas)
    return (
        mean,
        median,
        variance,
        CI_left,
        CI_right,
        mad,
        log2_mean,
        log2_median,
        log2_variance,
        abslog2_mean,
        abslog2_median,
        abslog2_variance,
    )


def parse_stan_output_initializer(thetas, lambdas):
    global g_thetas, g_lambdas
    g_thetas = thetas
    g_lambdas = lambdas

def parse_stan_output_worker_atacseq(line):
    global g_thetas, g_lambdas

    fields = line.rstrip().split()
    gene_id = fields[0]
    gene_thetas = g_thetas.get(gene_id)
    if not gene_thetas:
        return None

    lambdas_choice_gam = g_lambdas.loc[
        g_lambdas["peakID"] == gene_id, "gam_lambda"
    ].iloc[0]

    (
        mean,
        median,
        variance,
        left_CI,
        right_CI,
        mad,
        log2_mean,
        log2_median,
        log2_variance,
        abslog2_mean,
        abslog2_median,
        abslog2_variance,
    ) = summarize(gene_thetas, 0.05)
    log2_thetas = np.log2(np.array(gene_thetas))
    _, sum_prob_lambda_gam = computeBeastieScoreLog2(log2_thetas, lambdas_choice_gam)

    return (
        gene_id,
        round(mad, 3),
        round(median, 3),
        round(mean, 3),
        round(variance, 3),
        round(left_CI, 3),
        round(right_CI, 3),
        sum_prob_lambda_gam,
        log2_median,
        log2_mean,
        log2_variance,
        abslog2_median,
        abslog2_mean,
        abslog2_variance,
    )

def parse_stan_output_worker(line):
    global g_thetas, g_lambdas

    fields = line.rstrip().split()
    gene_id = fields[0]
    gene_thetas = g_thetas.get(gene_id)
    if not gene_thetas:
        return None
    lambdas_choice_gam = g_lambdas.loc[
        g_lambdas["peakID"] == gene_id, "gam_lambda"
    ].iloc[0]
    (
        mean,
        median,
        variance,
        left_CI,
        right_CI,
        mad,
        log2_mean,
        log2_median,
        log2_variance,
        abslog2_mean,
        abslog2_median,
        abslog2_variance,
    ) = summarize(gene_thetas, 0.05)
    log2_thetas = np.log2(np.array(gene_thetas))
    _, sum_prob_lambda_gam = computeBeastieScoreLog2(log2_thetas, lambdas_choice_gam)

    return (
        gene_id,
        round(mad, 3),
        round(median, 3),
        round(mean, 3),
        round(variance, 3),
        round(left_CI, 3),
        round(right_CI, 3),
        sum_prob_lambda_gam,
        log2_median,
        log2_mean,
        log2_variance,
        abslog2_median,
        abslog2_mean,
        abslog2_variance,
    )


def parse_stan_output_new(input_file, thetas_file, lambdas_file,atacseq):
    thetas = pickle.load(open(thetas_file, "rb"))
    # names = ['geneID','median_altratio','num_hets','totalRef','totalAlt','total_reads','predicted_lambda']
    lambdas = pd.read_csv(lambdas_file, delimiter="\t", header=0)

    with open(input_file, "rt") as IN, multiprocessing.Pool(
        initializer=parse_stan_output_initializer,
        initargs=(
            thetas,
            lambdas,
        ),
    ) as pool:
        if atacseq is not True:
            rows = pool.map(parse_stan_output_worker, IN, chunksize=1)
        else:
            rows = pool.map(parse_stan_output_worker_atacseq, IN, chunksize=1)
    if atacseq is True:
        df = pd.DataFrame(
            rows,
            columns=[
                "peakID",
                "median_abs_deviation",
                "posterior_median",
                "posterior_mean",
                "posterior_variance",
                "CI_left",
                "CI_right",
                "posterior_mass_support_ALT_gam",
                "log2_posterior_median",
                "log2_posterior_mean",
                "log2_posterior_variance",
                "abslog2_posterior_median",
                "abslog2_posterior_mean",
                "abslog2_posterior_variance",
            ],
        )
    else:
        df = pd.DataFrame(
            rows,
            columns=[
                "geneID",
                "median_abs_deviation",
                "posterior_median",
                "posterior_mean",
                "posterior_variance",
                "CI_left",
                "CI_right",
                "posterior_mass_support_ALT_gam",
                "log2_posterior_median",
                "log2_posterior_mean",
                "log2_posterior_variance",
                "abslog2_posterior_median",
                "abslog2_posterior_mean",
                "abslog2_posterior_variance",
            ],
        )

    return df


def run_model_worker_initializer(
    tmp_dir,
    model,
    sigma,
    WARMUP,
    KEEPER,
    phasing_method,
):
    global g_model, g_sigma, g_WARMUP, g_KEEPER, g_phasing_method
    g_model = model
    g_sigma = sigma
    g_WARMUP = WARMUP
    g_KEEPER = KEEPER
    g_phasing_method = phasing_method

    global g_tmp_output_file, g_stan_output_file, g_init_file
    g_tmp_output_file = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False).name
    g_stan_output_file = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False).name
    g_init_file = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False).name


def run_model_worker(line):
    global g_model, g_tmp_output_file, g_stan_output_file, g_init_file, g_sigma, g_WARMUP, g_KEEPER, g_phasing_method

    fields = line.rstrip().split()
    geneID, thetas = runModel(
        g_model,
        fields,
        g_tmp_output_file,
        g_stan_output_file,
        g_init_file,
        g_sigma,
        g_WARMUP,
        g_KEEPER,
        g_phasing_method,
    )
    return geneID, thetas


def save_raw_theta_parallel(
    input_filepath,
    out_filepath,
    model,
    sigma,
    WARMUP,
    KEEPER,
    phasing_method,
):
    processes = os.cpu_count()
    with tempfile.TemporaryDirectory() as tmp_dir, multiprocessing.Pool(
        processes=processes,
        initializer=run_model_worker_initializer,
        initargs=(tmp_dir, model, sigma, WARMUP, KEEPER, phasing_method),
    ) as pool, open(input_filepath, "rt") as input_file:
        items = []
        for item in pool.imap_unordered(run_model_worker, input_file):
            items.append(item)
            if len(items) % 500 == 0:
                logging.info(f".... Processed thetas for {len(items)} genes")
        logging.info(".... Processed thetas for all genes")

    model_theta = dict(items)
    pickle.dump(model_theta, open(out_filepath, "wb"))


def run(
    prefix,
    inFile,
    sigma,
    modelpath,
    out0,
    lambdas_file,
    WARMUP,
    KEEPER,
    phasing_method,
    ancestry,
    atacseq,
):
    if phasing_method != "nophasing":
        out_BEASTIE = "iBEASTIE"
    else:
        out_BEASTIE = "BEASTIE-fix-uniform"
    logging.debug(
        "Number of WARMUP samples is {0}, Number of posterior estimates is {1}".format(
            WARMUP, KEEPER
        )
    )
    out = os.path.join(out0, "output_pkl")
    out_path = os.path.join(out, out_BEASTIE)
    if not os.path.exists(out):
        os.makedirs(out)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    ###########################################################################################
    outname1 = "stan.pickle"
    out_path = os.path.join(out_path, "theta")
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    thetas_file = os.path.join(out_path, outname1)
    # step1
    if not os.path.isfile(thetas_file):
        save_raw_theta_parallel(
            inFile,
            thetas_file,
            modelpath,
            sigma,
            WARMUP,
            KEEPER,
            phasing_method,
        )
    logging.info(
        "...... Finshed running {0} and saved raw theta at : {1}".format(
            os.path.basename(modelpath), os.path.dirname(thetas_file)
        )
    )
    logging.info("...... Start parse_stan_output")
    # df = parse_stan_output(
    #     out0, prefix, inFile, thetas_file, lambdas_file, os.path.basename(modelpath)
    # )
    df = parse_stan_output_new(inFile, thetas_file, lambdas_file,atacseq)
    if atacseq is not True:
        df["ancestry"] = ancestry
    logging.info("...... Finish parse_stan_output")

    if "iBEASTIE" in os.path.basename(modelpath):
        modelname = "iBEASTIE"
    else:
        modelname = "BEASTIE_fix_uniform"
    parsed_output_path = os.path.join(out0, f"{prefix}_ASE_{modelname}.tsv")
    df.to_csv(
        parsed_output_path,
        sep="\t",
        header=True,
        index=False,
    )
    logging.info(f"...... Saved parsed output at {parsed_output_path}")

    return df, outname1
