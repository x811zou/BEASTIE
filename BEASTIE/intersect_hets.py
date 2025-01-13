#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import os
import logging
import statistics
from xxlimited import Xxo
import numpy as np
import pandas as pd
import sys


def Intersect_exonicHetSnps(
    parsed_mpileup_file,
    hetSNP_AF,
    read_length,
    min_totalcounts,
    min_singlecounts,
    hetSNP_intersect_pileup,
    atacseq
):
    df = pd.read_csv(parsed_mpileup_file, index_col=False, sep="\t")
    df_sub = df.drop(
        [
            "refAllele",
            "altAllele",
            "lowMAPQDepth",
            "lowBaseQDepth",
            "rawDepth",
            "otherCount",
        ],
        axis=1,
    )
    hetSNP_data = pd.read_csv(hetSNP_AF, index_col=False, sep="\t")
    hetSNP_data["contig"] = hetSNP_data["chrN"]
    hetSNP_data["position"] = hetSNP_data["pos"]
    out_base = os.path.dirname(hetSNP_intersect_pileup)
    out1 = "{0}/TEMP.filtered_beforeThinning.tsv".format(out_base)
    out11 = "{0}/TEMP.filtered_beforeThinning.simulation.tsv".format(out_base)
    out2 = "{0}/TEMP.filtered_afterThinning_beforeDropTrans.tsv".format(out_base)
    if len(df_sub):
        df2 = df_sub[df_sub["refCount"] >= int(min_singlecounts)]
        df3 = df2[df2["altCount"] >= int(min_singlecounts)]
        df4 = df3[df3["totalCount"] >= int(min_totalcounts)]
        df5 = df4[df4["if_SNP"] == "Y"]
        df6 = df5[df5["if_biallelic"] == "Y"]
        def custom_sort(chrom):
            """ Custom sort: numeric chromosomes first, followed by non-numeric ones. """
            try:
                # Try to convert to integer for numeric sorting
                return int(chrom), ''
            except ValueError:
                # If not numeric, return a high integer and the chromosome string for sorting
                return float('inf'), chrom

        df6["contig"] = df6["contig"].astype(str)  # Ensure all values are strings
        df6 = df6.sort_values(by="contig", key=lambda col: col.map(custom_sort))

        df6["position"] = df6["position"].astype(int)
        df6 = df6.drop_duplicates()
        hetSNP_data["contig"] = hetSNP_data["contig"].astype(str)
        hetSNP_data = hetSNP_data.sort_values(by="contig", key=lambda col: col.map(custom_sort))
        
        hetSNP_data["position"] = hetSNP_data["position"].astype(int)
        # take intersection between these two dataframes
        df_overlapped = pd.merge(
            hetSNP_data, df6, on=["contig", "position"], how="inner"
        )
        df_overlapped = df_overlapped.drop_duplicates()
        if atacseq is True:
            print("ATACseq data skipping double counting reads within a reads length")
            df_overlapped.sort_values(
                ["chrN", "pos","geneID", "peakID"],
                ascending=[True, True,True, True],
                inplace=True,
            )
            df_overlapped = df_overlapped.drop_duplicates()
            df_overlapped.sort_values(
                ["chrN", "pos", "peakID"], ascending=[True, True, True], inplace=True
            )
            df_overlapped_uni = df_overlapped.groupby(["chr", "pos"]).first().reset_index()
            df_overlapped_uni.to_csv(
                hetSNP_intersect_pileup, index=False, sep="\t", header=True
            )
        else:
            df_overlapped.sort_values(
                ["chrN", "geneID", "transcript_pos", "transcriptID"],
                ascending=[True, True, True, True],
                inplace=True,
            )
            df_overlapped.to_csv(out1, index=False, sep="\t", header=True)
            df_overlapped["snp_window"] = (
                df_overlapped["transcript_pos"]
                - df_overlapped.groupby(["geneID", "transcriptID"])[
                    "transcript_pos"
                ].transform(np.min)
            ) // int(read_length)
            df_overlapped = (
                df_overlapped.sort_values("totalCount", ascending=False)
                .groupby(["geneID", "transcriptID", "snp_window"])
                .nth(0)
            )
            df_overlapped.reset_index(inplace=True)
            df_overlapped.to_csv(out2, index=False, sep="\t", header=True)

            df_overlapped = df_overlapped.drop(
                ["snp_window", "transcript_pos", "transcriptID"], axis=1
            )
            df_overlapped = df_overlapped.drop_duplicates()
            df_overlapped.sort_values(
                ["chrN", "pos", "geneID"], ascending=[True, True, True], inplace=True
            )
            df_overlapped_uni = df_overlapped.groupby(["chr", "pos"]).first().reset_index()
            df_overlapped_uni.to_csv(
                hetSNP_intersect_pileup, index=False, sep="\t", header=True
            )


def summary_statistics(data, title):
    print(title + " statistics:")
    if isinstance(data, list):
        print(
            "  #: {:-2} \n  Max: {:-2} \n  Min: {:-2} \n  Mean: {:-2} \n  Median: {:-2} \n  Variance: {:-2} \n  Std: {:-2} \n  25% quantile: {:-2} \n  75% quantile: {:-2}  \n  75% quantile + IQR*1.5: {:-2} \n  90% quantile: {:-2}  \n  91% quantile: {:-2}  \n  95% quantile: {:-2} \n  98% quantile: {:-2} \n  99% quantile: {:-2}".format(
                len(data),
                round(max(data), 3),
                round(min(data), 3),
                round(statistics.mean(data), 3),
                round(statistics.median(data), 3),
                round(statistics.variance(data), 3),
                round(statistics.stdev(data), 3),
                round(np.quantile(data, 0.25), 3),
                round(np.quantile(data, 0.75), 3),
                round(
                    np.quantile(data, 0.75)
                    + (np.quantile(data, 0.75) - np.quantile(data, 0.25)) * 1.5,
                    3,
                ),
                round(np.quantile(data, 0.90), 3),
                round(np.quantile(data, 0.91), 3),
                round(np.quantile(data, 0.95), 3),
                round(np.quantile(data, 0.98), 3),
                round(np.quantile(data, 0.99), 3),
            )
        )
    else:
        print("Not a list")
