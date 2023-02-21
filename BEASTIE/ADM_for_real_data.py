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


def fixed_simulator(total_counts, num_iterations, estimate):
    alternate_counts = np.random.binomial(
        total_counts, 0.5, size=(num_iterations, len(total_counts))
    )
    reference_counts = total_counts - alternate_counts
    AA = np.abs(alternate_counts - reference_counts) / (
        alternate_counts + reference_counts
    )

    means = np.mean(AA, axis=1)
    pval = 1 - percentileofscore(means, estimate) / 100
    return pval


def getAA(fields):
    if len(fields) >= 4:
        esti = []
        Mreps = int(fields[1])
        # Mreps = hets
        # p = float(theta)/(float(theta)+1)
        total_counts = []
        for rep in range(Mreps):
            A = float(fields[2 + rep * 2])
            R = float(fields[3 + rep * 2])
            estimate = abs(R - A) / (A + R)
            esti.append(estimate)
            total_count = A + R
            total_counts.append(int(total_count))
        avg_esti = statistics.mean(esti)
        pval = fixed_simulator(total_counts, 1000, avg_esti)
    return round(avg_esti, 3), round(pval, 3)


def process_line(line):
    fields = line.rstrip().split()
    gene = fields[0]
    esti, pval = getAA(fields)
    return gene, esti, pval


def run(in_path, out_path):
    # parallel
    with open(in_path, "rt") as IN, multiprocessing.Pool() as pool:
        rows = pool.map(process_line, IN)

    # serial
    # with open(in_path, "rt") as IN:
    #     rows = [process_line(line) for line in IN]

    logging.info("....... start saving ADM file")
    ADM_df = pd.DataFrame(
        rows,
        columns=["geneID", "ADM_esti", "ADM_pval"],
    )
    ADM_df.to_csv(out_path, sep="\t", header=True, index=False)

    return ADM_df
