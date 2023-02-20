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


def getMaxProb_RMSE(thetas):
    p_less1 = len([i for i in thetas if i < 1]) / len(thetas)
    p_more1 = 1 - p_less1
    max_prob1 = max(p_less1, p_more1)
    # 2. transform thetas, and then calculate proportion
    thetas_log2 = [log2(x) for x in thetas]
    p_less2 = len([i for i in thetas_log2 if i < 0]) / len(thetas_log2)
    p_more2 = 1 - p_less2
    max_prob2 = max(p_less2, p_more2)
    return max_prob1  # , RMSE


def getMaxProb_lambda(thetas, Lambda):
    # 1. no transformation
    one_over_Lambda = float(1 / float(Lambda))
    # changes maded in 07/22
    min_l = one_over_Lambda  # min(Lambda,one_over_Lambda)
    max_l = Lambda  # max(Lambda,one_over_Lambda)
    #
    p_less1 = len([i for i in thetas if i < min_l]) / len(thetas)
    p_more1 = len([i for i in thetas if i > max_l]) / len(thetas)
    lambda_prob1 = max(p_less1, p_more1)
    # 2. transform thetas, and then calculate proportion
    thetas_log2 = [math.log2(x) for x in thetas]
    p_less2 = len([i for i in thetas_log2 if i < math.log2(one_over_Lambda)]) / len(
        thetas
    )
    p_more2 = len([i for i in thetas_log2 if i > math.log2(float(Lambda))]) / len(
        thetas
    )
    lambda_prob2 = max(p_less2, p_more2)
    # 3. sum tail
    lambda_sum1 = p_less1 + p_more1
    # 4. sum tail  transform thetas, and then calculate proportion
    lambda_sum2 = p_less2 + p_more2
    return round(lambda_prob2, 3), round(lambda_sum2, 3)


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
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        parser = StanParser(stan_output_file)
        thetas = parser.getVariable("theta")
        return geneID, thetas
    else:
        logging.error("lines with no enough elements")


def parse_lambda_validation_simulation(thetas, alphas, lambdas_file, lm):
    with open(lambdas_file, "rt") as IN:
        for idx, line in enumerate(IN):
            if idx == lm - 1:
                fields = line.rstrip().split()
                # print(fields)
                print("model %s - alpha %s - lambda: %s" % (lm, alphas[0], fields[0]))
                prob_lambda1, _, sum_lambda1, _ = getMaxProb_lambda(
                    thetas, float(fields[0])
                )
                print("model %s - alpha %s - lambda: %s" % (lm, alphas[1], fields[1]))
                prob_lambda2, _, sum_lambda2, _ = getMaxProb_lambda(
                    thetas, float(fields[1])
                )
                print("model %s - alpha %s - lambda: %s" % (lm, alphas[2], fields[2]))
                prob_lambda3, _, sum_lambda3, _ = getMaxProb_lambda(
                    thetas, float(fields[2])
                )
                return (
                    prob_lambda1,
                    sum_lambda1,
                    prob_lambda2,
                    sum_lambda2,
                    prob_lambda3,
                    sum_lambda3,
                )


# def getMedian(thetas):
#     # Precondition: thetas is already sorted
#     thetas.sort()
#     n = len(thetas)
#     mid = int(n / 2)
#     if n % 2 == 0:
#         return (thetas[mid - 1] + thetas[mid]) / 2.0
#     return thetas[mid]


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


def parse_stan_output(out, prefix, input_file, out1, KEEPER, lambdas_file, model):
    thetas = pickle.load(
        open(
            out1,
            "rb",
        )
    )
    lambdas = pd.read_csv(
        lambdas_file, delimiter="\t", header=0
    )  # names = ['geneID','median_altratio','num_hets','totalRef','totalAlt','total_reads','predicted_lambda']

    prob_sum_lambda_gam4 = []
    prob_sum_lambda_gam3 = []
    model_theta_med = []  # 150
    model_theta_mean = []  # 150
    model_theta_var = []  # 150
    model_log2_theta_med = []  # 150
    model_log2_theta_mean = []  # 150
    model_log2_theta_var = []  # 150
    model_abslog2_theta_med = []  # 150
    model_abslog2_theta_mean = []  # 150
    model_abslog2_theta_var = []  # 150
    model_mad = []
    CI_left = []
    CI_right = []
    geneID = []
    with open(input_file, "rt") as IN:
        # i=0
        for line in IN:
            fields = line.rstrip().split()
            ID = fields[0]
            # read the ith geneID
            # j=i+int(KEEPER)-1
            gene_thetas = thetas.get(ID)
            # print(gene_thetas)
            # print(len(gene_thetas))
            if len(gene_thetas) > 1:
                # print(">>>> record")
                geneID.append(ID)
                lambdas_choice_gam4 = lambdas.loc[
                    lambdas["geneID"] == ID, "gam4_lambda"
                ].iloc[0]
                lambdas_choice_gam3 = lambdas.loc[
                    lambdas["geneID"] == ID, "gam3_lambda"
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
                # print(f"mad {mad}")
                max_prob = getMaxProb_RMSE(gene_thetas)
                max_prob_lambda, sum_prob_lambda_gam4 = getMaxProb_lambda(
                    gene_thetas, lambdas_choice_gam4
                )
                _, sum_prob_lambda_gam3 = getMaxProb_lambda(
                    gene_thetas, lambdas_choice_gam3
                )
                # i=i+int(KEEPER)
                prob_sum_lambda_gam4.append(sum_prob_lambda_gam4)
                prob_sum_lambda_gam3.append(sum_prob_lambda_gam3)
                CI_left.append(round(left_CI, 3))
                CI_right.append(round(right_CI, 3))
                model_theta_mean.append(mean)
                model_theta_med.append(median)
                model_theta_var.append(variance)
                model_mad.append(round(mad, 3))
                model_log2_theta_mean.append(log2_mean)
                model_log2_theta_med.append(log2_median)
                model_log2_theta_var.append(log2_variance)
                model_abslog2_theta_mean.append(abslog2_mean)
                model_abslog2_theta_med.append(abslog2_median)
                model_abslog2_theta_var.append(abslog2_variance)
    df = {
        "geneID": geneID,
        "median_abs_deviation": model_mad,
        "posterior_median": model_theta_med,
        "posterior_mean": model_theta_mean,
        "posterior_variance": model_theta_var,
        "CI_left": CI_left,
        "CI_right": CI_right,
        "posterior_mass_support_ALT_gam4": prob_sum_lambda_gam4,
        "posterior_mass_support_ALT_gam3": prob_sum_lambda_gam3,
        "log2_posterior_median": model_log2_theta_med,
        "log2_posterior_mean": model_log2_theta_mean,
        "log2_posterior_variance": model_log2_theta_var,
        "abslog2_posterior_median": model_abslog2_theta_med,
        "abslog2_posterior_mean": model_abslog2_theta_mean,
        "abslog2_posterior_variance": model_abslog2_theta_var,
    }
    df = pd.DataFrame(df)
    # df["posterior_mean"] = df["posterior_mean"].apply(
    #     lambda x: round(x, 3 - int(floor(log10(abs(x)))))
    # )
    # df["posterior_median"] = df["posterior_median"].apply(
    #     lambda x: round(x, 3 - int(floor(log10(abs(x)))))
    # )
    # df["posterior_variance"] = df["posterior_variance"].apply(
    #     lambda x: round(x, 3 - int(floor(log10(abs(x)))))
    # )
    df["posterior_mean"] = df["posterior_mean"].apply(lambda x: round(x, 3))
    df["posterior_median"] = df["posterior_median"].apply(lambda x: round(x, 3))
    df["posterior_variance"] = df["posterior_variance"].apply(lambda x: round(x, 3))
    if "iBEASTIE" in model:
        modelname = "iBEASTIE"
    else:
        modelname = "BEASTIE_fix_uniform"
    df.to_csv(
        out + "/" + prefix + "_ASE_" + modelname + ".tsv",
        sep="\t",
        header=True,
        index=False,
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
    prefix, inFile, sigma, modelpath, out0, lambdas_file, WARMUP, KEEPER, phasing_method
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
    out1 = os.path.join(out_path, outname1)
    # step1
    if not os.path.isfile(out1):
        save_raw_theta_parallel(
            inFile,
            out1,
            modelpath,
            sigma,
            WARMUP,
            KEEPER,
            phasing_method,
        )
    logging.info(
        "...... Finshed running {0} and saved raw theta at : {1}".format(
            os.path.basename(modelpath), os.path.dirname(out1)
        )
    )
    logging.info("...... Start parse_stan_output")
    df = parse_stan_output(
        out0, prefix, inFile, out1, KEEPER, lambdas_file, os.path.basename(modelpath)
    )
    logging.info("...... Finish parse_stan_output")
    # step2
    return df, outname1
