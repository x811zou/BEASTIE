#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import os
import pickle
import logging
import statistics
import numpy as np
import pandas as pd
from math import floor, log10, log2
from .misc_tools.StanParser import StanParser
from .helpers import runhelper
from scipy import stats
import statistics
import sys


def writeInitializationFile(filename):
    OUT = open(filename, "wt")
    print("theta <- 1", file=OUT)
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
    # thetas_log2 = [math.log2(x) for x in thetas]
    # p_less2 = len([i for i in thetas_log2 if i < math.log2(one_over_Lambda)])/len(thetas)
    # p_more2 = len([i for i in thetas_log2 if i > math.log2(float(Lambda))])/len(thetas)
    # lambda_prob2 = max(p_less2,p_more2)
    # 3. sum tail
    lambda_sum1 = p_less1 + p_more1
    # 4. sum tail  transform thetas, and then calculate proportion
    # lambda_sum2 = p_less2 + p_more2
    return round(lambda_prob1, 3), round(lambda_sum1, 3)


def runModel(
    model, fields, tmp_output_file, stan_output_file, init_file, sigma, WARMUP, KEEPER
):
    if len(fields) >= 4:
        geneID = str(fields[0])
        # logging.debug(geneID)
        writeInputsFile_i(fields, tmp_output_file, sigma)
        writeInitializationFile(init_file)
        cmd = (
            "%s sample num_samples=%s num_warmup=%s data file=%s init=%s output file=%s refresh=0"
            % (model, KEEPER, WARMUP, tmp_output_file, init_file, stan_output_file)
        )
        # print(cmd)
        # logging.debug(cmd)
        runhelper(cmd)  # Parse MCMC output
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
    median = statistics.median(thetas)
    # print(f"median {median}")
    variance = np.var(thetas)
    # print(f"variance {variance}")
    # print(f"length of thetas {len(thetas)}")
    n = len(thetas)
    thetas.sort()
    # print(f"length of sorted thetas {len(thetas)}")
    CI_left, CI_right = getCredibleInterval(thetas, alpha, n)
    # mad = stats.median_abs_deviation(thetas, scale="normal")
    # mad = stats.median_abs_deviation(thetas, scale=1/1.4826)
    mad = stats.median_absolute_deviation(thetas)
    return median, variance, CI_left, CI_right, mad


def parse_stan_output(out, prefix, input_file, out1, KEEPER, lambdas_file):
    thetas = pickle.load(
        open(
            out1,
            "rb",
        )
    )
    lambdas = pd.read_csv(
        lambdas_file, delimiter="\t", header=0
    )  # names = ['geneID','median_altratio','num_hets','totalRef','totalAlt','total_reads','predicted_lambda']

    prob_sum_lambda = []
    model_theta_med = []  # 150
    model_theta_var = []  # 150
    model_mad = []  # 150
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
                lambdas_choice = lambdas.loc[lambdas["geneID"] == ID].iloc[
                    0, 6
                ]  # has to change here
                # log_lambda=(log(alpha/(1-alpha)) -(as.numeric(model$coefficients[1])+as.numeric(model$coefficients[3])*as.integer(totalCount)))/as.numeric(model$coefficients[2]))
                # predicted_lambda = exp(log_lambda)
                # predicted_lambda_plus1 = predicted_lambda_1+1)

                median, variance, left_CI, right_CI, mad = summarize(gene_thetas, 0.05)
                # print(f"mad {mad}")
                max_prob = getMaxProb_RMSE(gene_thetas)
                max_prob_lambda, sum_prob_lambda = getMaxProb_lambda(
                    gene_thetas, lambdas_choice
                )
                # i=i+int(KEEPER)
                prob_sum_lambda.append(sum_prob_lambda)
                CI_left.append(round(left_CI, 3))
                CI_right.append(round(right_CI, 3))
                model_theta_med.append(median)
                model_theta_var.append(variance)
                model_mad.append(round(mad, 3))

    df = {
        "geneID": geneID,
        "median_abs_deviation": model_mad,
        "posterior_median": model_theta_med,
        "posterior_variance": model_theta_var,
        "CI_left": CI_left,
        "CI_right": CI_right,
        "posterior_mass_support_ALT": prob_sum_lambda,
    }
    df = pd.DataFrame(df)
    df["posterior_median"] = df["posterior_median"].apply(
        lambda x: round(x, 3 - int(floor(log10(abs(x)))))
    )
    df["posterior_variance"] = df["posterior_variance"].apply(
        lambda x: round(x, 3 - int(floor(log10(abs(x)))))
    )
    df.to_csv(
        out + "/" + prefix + "_ASE_ibeastie.tsv", sep="\t", header=True, index=False
    )
    return df


def save_raw_theta(
    out0,
    models,
    input_file,
    tmp_output_file,
    stan_output_file,
    init_file,
    sigma,
    WARMUP,
    KEEPER,
):
    model_theta = {}  # 150
    with open(input_file, "rt") as IN:
        i = 0
        for line in IN:
            i += 1
            # logging.debug(line)
            fields = line.rstrip().split()
            # logging.debug(fields)
            geneID, thetas = runModel(
                models,
                fields,
                tmp_output_file,
                stan_output_file,
                init_file,
                sigma,
                WARMUP,
                KEEPER,
            )
            # model_theta.extend(thetas)
            model_theta[geneID] = thetas
    logging.debug(
        ".... Number of genes : {0}, and length of thetas: {1}".format(
            i, len(model_theta)
        )
    )
    pickle.dump(model_theta, open(out0, "wb"))
    # print(model_theta)


def run(
    prefix,
    inFile,
    sigma,
    alpha,
    models,
    out0,
    lambdas_file,
    WARMUP,
    KEEPER,
    either_cov,
    total_cov,
):
    logging.debug(
        "Number of WARMUP samples is {0}, Number of posterior estimates is {1}".format(
            WARMUP, KEEPER
        )
    )
    tmpFile = "tmp_output.txt"
    initFile = "initialization_stan.txt"
    outFile = "stan_output.txt"
    out = os.path.join(out0, "output_pkl")
    out_path = os.path.join(out, "ibeastie")
    if not os.path.exists(out):
        os.makedirs(out)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    tmp_output_file = os.path.join(out_path, tmpFile)
    init_file = os.path.join(out_path, initFile)
    stan_output_file = os.path.join(out_path, outFile)
    ###########################################################################################
    outname1 = "stan.pickle"
    out_path = os.path.join(out_path, "theta")
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    out1 = os.path.join(out_path, outname1)
    # step1
    if not os.path.isfile(out1):
        save_raw_theta(
            out1,
            models,
            inFile,
            tmp_output_file,
            stan_output_file,
            init_file,
            sigma,
            WARMUP,
            KEEPER,
        )
    logging.info(
        "...... Finshed running {0} and saved raw theta at : {1}".format(
            os.path.basename(models), os.path.dirname(out1)
        )
    )
    df = parse_stan_output(out0, prefix, inFile, out1, KEEPER, lambdas_file)
    # step2
    return df, outname1
