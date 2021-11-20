import os
import logging
import statistics
import numpy as np
import pandas as pd


def Intersect_exonicHetSnps(
    parsed_mpileup_file,
    hetSNP_AF,
    read_length,
    min_totalcounts,
    min_singlecounts,
    hetSNP_intersect_unique,
    hetSNP_intersect_unique_forlambda_file,
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

    out_base = os.path.splitext(hetSNP_intersect_unique)[0]
    out1 = "{0}_filtered_beforeThinning.TEMP.tsv".format(out_base)
    out2 = "{0}_filtered_afterThinning_beforeDropTrans.TEMP.tsv".format(out_base)
    if len(df_sub):
        df2 = df_sub[df_sub["refCount"] >= int(min_singlecounts)]
        df3 = df2[df2["altCount"] >= int(min_singlecounts)]
        df4 = df3[df3["totalCount"] >= int(min_totalcounts)]
        df5 = df4[df4["if_SNP"] == "Y"]
        df6 = df5[df5["if_biallelic"] == "Y"]
        df6["contig"] = df6["contig"].astype(int)
        df6["position"] = df6["position"].astype(int)
        df6 = df6.drop_duplicates()
        hetSNP_data["contig"] = hetSNP_data["contig"].astype(int)
        hetSNP_data["position"] = hetSNP_data["position"].astype(int)
        # take intersection between these two dataframes
        df_overlapped = pd.merge(
            hetSNP_data, df6, on=["contig", "position"], how="inner"
        )
        df_overlapped = df_overlapped.drop_duplicates()
        df_overlapped.sort_values(
            ["chrN", "transcript_pos", "geneID", "transcriptID"],
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
            hetSNP_intersect_unique, index=False, sep="\t", header=True
        )

        df6_med = (
            df_overlapped_uni.groupby(["geneID"])
            .agg({"altRatio": "median"})
            .reset_index()
            .rename(columns={"altRatio": "median.altRatio"})
        )
        df6_nhets = (
            df_overlapped_uni.groupby(by="geneID")
            .agg({"position": pd.Series.nunique})
            .rename(columns={"position": "number.of.hets"})
        )
        df6_totalref = (
            df_overlapped_uni.groupby(by="geneID")
            .agg({"refCount": "sum"})
            .reset_index()
            .rename(columns={"totalRef": "total refAllele"})
        )
        df6_totalalt = (
            df_overlapped_uni.groupby(by="geneID")
            .agg({"altCount": "sum"})
            .reset_index()
            .rename(columns={"totalAlt": "total altAllele"})
        )
        df_summary_1 = pd.merge(df6_med, df6_nhets, on=["geneID"], how="inner")
        df_summary_2 = pd.merge(df_summary_1, df6_totalref, on=["geneID"], how="inner")
        df_summary_3 = pd.merge(df_summary_2, df6_totalalt, on=["geneID"], how="inner")
        df_summary_3["totalCount"] = df_summary_3["refCount"] + df_summary_3["altCount"]
    else:
        logging.error("file no lines")
    logging.info(
        "..... Aggregated pasred mpile up data dim {0} saved to {1}".format(
            df_overlapped_uni.size, hetSNP_intersect_unique
        )
    )
    df_summary_3.to_csv(
        hetSNP_intersect_unique_forlambda_file, index=False, sep="\t", header=True
    )
    logging.info(
        "..... Number of genes in this file for lambda prediction is {0} saved to {1}".format(
            df_summary_3.size, hetSNP_intersect_unique_forlambda_file
        )
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
