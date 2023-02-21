#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import logging
import multiprocessing
import os
import os.path
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# def create_binomial_library(depth):
#     count_p = {}
#     for i in range(int(depth)+1):
#         count_p[i] = binom.cdf(i, int(depth), 0.5)
#     return count_p


def getBaseline(fields, depth):
    if len(fields) >= 4:
        ############
        # count_p = create_binomial_library(depth)
        Mreps = int(fields[1])
        total_AR = 0
        for rep in range(Mreps):
            A = float(fields[2 + rep * 2])  # +1
            R = float(fields[3 + rep * 2])  # +1
            # logging.info('... {} site, A-{} R-{}'.format(rep,A,R))
            if A <= R:
                min_AR = A
                max_AR = R
            else:
                min_AR = R
                max_AR = A
            if total_AR < A + R:
                max_A = A
                max_R = R
                total_AR = A + R
            if str(rep) == "0":
                FS_AAR = A / (A + R)
                if FS_AAR == 1:
                    FS_AAR = FS_AAR - 0.001
                FS_esti = FS_AAR / (1 - FS_AAR)
                # logging.info('... FS_esti - {}'.format(FS_esti))
                FS_prob = stats.binom_test(A, A + R, p=0.5, alternative="two-sided")
                # logging.info('... FS_prob - {}'.format(FS_prob))
            # esti3
            base3 = abs(0.5 - max_AR / (A + R))
            # esti4
            diff_AR = float(abs(0.5 - max_AR / (A + R)))
        MS_AAR = max_A / (max_A + max_R)
        if MS_AAR == 1:
            MS_AAR = MS_AAR - 0.001
        MS_esti = MS_AAR / (1 - MS_AAR)
        # logging.info('... MS_esti - {}'.format(MS_esti))
        MS_prob = stats.binom_test(max_A, max_A + max_R, p=0.5, alternative="two-sided")
        # logging.info('... MS_prob - {}'.format(MS_prob))
        return (
            round(FS_esti, 3),
            round(FS_prob, 3),
            round(MS_esti, 3),
            round(MS_prob, 3),
        )
    else:
        return (None, None, None, None)


def getBaseline_pooled(fields, depth, hets):
    if len(fields) >= 4:
        base_thetas = []
        Mreps = int(fields[1])
        pooled_A = 0
        pooled_R = 0
        pooled_min = 0
        for rep in range(Mreps):
            A = float(fields[2 + rep * 2])  # +1
            R = float(fields[3 + rep * 2])  # +1
            # logging.info('... {0} site: A-{1} R-{2}'.format(rep,A,R))
            pooled_A = pooled_A + A
            pooled_R = pooled_R + R
            pooled_min = pooled_min + min(A, R)
        sum_AR = pooled_A + pooled_R
        # logging.info('... {0} sites pooled: A-{1} R-{2}'.format(rep,pooled_A,pooled_R))
        # logging.info('... {0} sites pooled: minAR-{1} maxAR-{2}'.format(rep,pooled_min,sum_AR-pooled_min))
        # Naive Sum (NS): summing up the counts for A assuming phasing is correct
        NS_AAR = pooled_A / sum_AR
        if NS_AAR == 1:  # if REF has 0 count
            NS_AAR = NS_AAR - 0.001
        NS_esti = NS_AAR / (1 - NS_AAR)
        NS_prob = stats.binom_test(
            pooled_A, pooled_A + pooled_R, p=0.5, alternative="two-sided"
        )
        # Pseudo phasing: lower counts assuming for ALT allele, high counts assuming for REF allele
        pseudo_esti = pooled_min / sum_AR
        # pseudo_esti = pseudo_AAR/(1-pseudo_AAR)
        pseudo_p = stats.binom_test(
            pooled_min, pooled_A + pooled_R, p=0.5, alternative="less"
        )
        # logging.info('... NS_esti: {0}, NS_prob: {1}'.format(NS_esti,NS_prob))
        # logging.info('... pseudo_esti: {0}, pseudo_prob: {1}'.format(pseudo_esti,pseudo_p))
        return (
            round(NS_esti, 3),
            round(NS_prob, 3),
            round(pseudo_esti, 3),
            round(pseudo_p, 3),
        )
    else:
        return (None, None, None, None)


def worker(line):
    fields = line.rstrip().split()
    geneID = fields[0]
    h = int(fields[1])
    d = int(fields[2]) + int(fields[3])
    FS_esti, FS_prob, MS_esti, MS_prob = getBaseline(fields, d)
    NS_esti, NS_prob, pseudo_esti, pseudo_p = getBaseline_pooled(fields, d, h)
    return (
        geneID,
        FS_esti,
        FS_prob,
        NS_esti,
        NS_prob,
        pseudo_esti,
        pseudo_p,
        MS_esti,
        MS_prob,
    )


def run(inFile):
    rows = []
    with open(inFile, "rt") as IN, multiprocessing.Pool() as pool:
        rows = pool.map(worker, IN)

    binomial_df = pd.DataFrame(
        rows,
        columns=[
            "geneID",
            "FirstSite_esti",
            "FirstSite_pval",
            "NaiveSum_esti",
            "NaiveSum_pval",
            "Pseudo_esti",
            "Pseudo_pval",
            "MajorSite_esti",
            "MajorSite_pval",
        ],
    )
    return binomial_df
