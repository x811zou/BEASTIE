#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import logging
import multiprocessing
import os
import os.path
import pickle
import statistics

import numpy as np
import pandas as pd
from scipy.stats import percentileofscore


def AA_estimate(A, R):
    AA = abs(A - R) / (A + R)
    return AA


def fixed_simulator(D, M, N, estimate):
    esti = []
    for k in range(N):
        AR = []
        for i in range(M):
            A = np.random.binomial(D[i], 0.5)
            R = D[i] - A
            AR.append(AA_estimate(A, R))

        esti.append(statistics.mean(AR))
    # find the tail
    pval = 1 - percentileofscore(esti, estimate) / 100
    # print("pvalue: "+str(pval))
    return pval


def getAA(fields):
    if len(fields) >= 4:
        esti = []
        Mreps = int(fields[1])
        # Mreps = hets
        # p = float(theta)/(float(theta)+1)
        depth = []
        for rep in range(Mreps):
            A = float(fields[2 + rep * 2])
            R = float(fields[3 + rep * 2])
            estimate = abs(R - A) / (A + R)
            esti.append(estimate)
            total_counts = A + R
            depth.append(int(total_counts))
        avg_esti = statistics.mean(esti)
        hets = Mreps
        pval = fixed_simulator(depth, hets, 1000, avg_esti)
    return round(avg_esti, 3), round(pval, 3)


def process_line(line):
    fields = line.rstrip().split()
    gene = fields[0]
    esti, pval = getAA(fields)
    return gene, esti, pval


def run(in_path, out_path):
    with open(in_path, "rt") as IN, multiprocessing.Pool() as pool:
        rows = pool.map(process_line, IN)

    logging.info("....... start saving ADM file")
    ADM_df = pd.DataFrame(
        rows,
        columns=["geneID", "ADM_esti", "ADM_pval"],
    )
    ADM_df.to_csv(out_path, sep="\t", header=True, index=False)

    return ADM_df
