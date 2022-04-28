#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

import dataclasses
from pickle import NONE
import re
import os
import sys
import logging
import scipy.stats
import pandas as pd
import numpy as np
from math import floor, log10, pi, isnan
from scipy.stats import binom_test, fisher_exact


def change_phasing(data):
    row_n = data.shape[0]
    data["patCount"] = data["refCount"]
    data["matCount"] = data["altCount"]
    reference_row = data.iloc[0]
    if not isnan(int(reference_row["e_paternal"])):
        for x in range(1, row_n):
            row = data.iloc[x]
            if int(row["e_paternal"]) != int(reference_row["e_paternal"]) and not isnan(
                int(row["e_paternal"])
            ):
                data.at[x, "patCount"] = row["altCount"]
                data.at[x, "matCount"] = row["refCount"]
    return data


def filter_alignBias(p_cutoff, snp_input, simulator_df=None):
    base_out = os.path.splitext(snp_input)[0]
    hetSNP_intersect_unique = pd.read_csv(
        snp_input, sep="\t", header=0, index_col=False
    )
    # genotyping error fisher exact test score
    grouped_df = (
        hetSNP_intersect_unique.groupby("geneID").apply(process_gene).reset_index()
    )
    grouped_df[
        [
            "chrN",
            "pos",
            "refCount",
            "altCount",
            "max_rest_data",
            "min_reast_data",
            "FishTest",
        ]
    ] = pd.DataFrame(grouped_df[0].tolist(), index=hetSNP_intersect_unique.index)
    grouped_df_sub = grouped_df.drop(grouped_df.columns[[1, 2]], axis=1)
    new_df = pd.merge(
        hetSNP_intersect_unique,
        grouped_df_sub,
        how="inner",
        on=["chrN", "pos", "geneID", "refCount", "altCount"],
    )
    ################# 1. filter out variants with alignment bias first
    if simulator_df is not None:
        beforeFilter_filename = f"{base_out}_alignBiasbeforeFilter.tsv"
        afterFilter_filename = f"{base_out}_alignBiasFiltered.tsv"
        # alignment bias binomial test score
        overlapped_variants = new_df.merge(simulator_df, on=["chr", "pos"], how="inner")
        # saving
        overlapped_variants.to_csv(
            beforeFilter_filename, index=False, sep="\t", header=True
        )
        Notoverlapped_variants = hetSNP_intersect_unique.drop_duplicates().merge(
            simulator_df.drop_duplicates(),
            on=["chr", "pos"],
            how="left",
            indicator=True,
        )
        Notoverlapped_variants = Notoverlapped_variants[
            Notoverlapped_variants["_merge"] == "left_only"
        ]
        Notoverlapped_variants_filename = f"{base_out}_notoverlapped.tsv"
        Notoverlapped_variants.to_csv(
            Notoverlapped_variants_filename, index=False, sep="\t", header=True
        )
        logging.debug(
            "{0} het SNPs in input data, {1} het SNPs in simulation data, overlapped het SNPs {2} in real data, not overlapped {3}".format(
                hetSNP_intersect_unique.shape[0],
                simulator_df.shape[0],
                overlapped_variants.shape[0],
                Notoverlapped_variants.shape[0],
            )
        )
        simulator_df_biased = overlapped_variants[
            overlapped_variants["alt_binomial_p"] > p_cutoff
        ]
        debiased_df = simulator_df_biased[simulator_df_biased["FishTest"] > p_cutoff]
        logging.debug(
            "{0} ({1}%) overlapped het SNPs pass alignment bias filtering p-val > {2} \n{3} ({4}%) het SNPs pass genotyping error filtering p-val > {5}".format(
                simulator_df_biased.shape[0],
                round(
                    simulator_df_biased.shape[0]
                    / hetSNP_intersect_unique.shape[0]
                    * 100,
                    2,
                ),
                p_cutoff,
                debiased_df.shape[0],
                round(
                    debiased_df.shape[0] / simulator_df_biased.shape[0] * 100,
                    2,
                ),
                p_cutoff,
            )
        )
        gene_df_filtered = debiased_df[
            [
                "chr",
                "chrN",
                "pos",
                "rsid",
                "AF",
                "geneID",
                "genotype",
                "refCount",
                "altCount",
                "totalCount",
                "altRatio",
                "alt_binomial_p",
                "FishTest",
            ]
        ]
        gene_df_filtered.to_csv(
            afterFilter_filename, index=False, sep="\t", header=True
        )
        filename = afterFilter_filename
    else:
        debiased_df = simulator_df_biased[simulator_df_biased["FishTest"] > p_cutoff]
        hetSNP_intersect_unique_subset = debiased_df[
            [
                "chr",
                "chrN",
                "pos",
                "rsid",
                "AF",
                "geneID",
                "genotype",
                "refCount",
                "altCount",
                "totalCount",
                "altRatio",
                "alt_binomial_p",
                "FishTest",
            ]
        ]
        afterFilter_filename = f"{base_out}_NoalignBiasFiltered.tsv"
        hetSNP_intersect_unique_subset.to_csv(
            afterFilter_filename, index=False, sep="\t", header=True
        )
        filename = afterFilter_filename
        logging.debug(
            "{0} has {1} het SNPs, after genotyping error filtering, {2} has {3} het SNPs, no alignment bias filtering".format(
                os.path.basename(snp_input),
                hetSNP_intersect_unique.shape[0],
                os.path.basename(afterFilter_filename),
                hetSNP_intersect_unique_subset.shape[0],
            )
        )
    return filename


def process_gene(data):
    data_sub = data[["chrN", "pos", "refCount", "altCount"]]
    data_rest = data_sub[(data_sub["refCount"] != 0) & (data_sub["altCount"] != 0)]
    size = data_rest.shape[0]
    if size == 0:
        min_rest_data = "NA"
        max_rest_data = "NA"
    else:
        data_rest["max_count"] = data_rest[["refCount", "altCount"]].max(axis=1)
        data_rest["min_count"] = data_rest[["refCount", "altCount"]].min(axis=1)
        min_rest_data = data_rest["min_count"].sum()
        max_rest_data = data_rest["max_count"].sum()

    def process_row(row):
        fishTest = 100
        if ((row["refCount"] == 0) or (row["altCount"] == 0)) and (size > 0):
            zero_counts = row[["refCount", "altCount"]].min()
            nonzero_counts = row[["refCount", "altCount"]].max()
            table = np.array(
                [[nonzero_counts, zero_counts], [max_rest_data, min_rest_data]]
            )
            oddsr, p = fisher_exact(table, alternative="two-sided")
            fishTest = p
        return (
            row["chrN"],
            row["pos"],
            row["refCount"],
            row["altCount"],
            max_rest_data,
            min_rest_data,
            fishTest,
        )

    return data.apply(process_row, axis=1)


def re_allocateReads(
    alignBiasfiltered_filename,
    shapeit2_input,
    version,
    filename,
    filename_cleaned,
    phase_difference_file=None,
):
    # read data after alignment bias filtering
    hetSNP_intersect_unique_filtered = pd.read_csv(
        alignBiasfiltered_filename, sep="\t", header=0, index_col=False
    )
    # phasing
    if os.path.isfile(filename_cleaned):
        logging.info(
            "....... phased data is pre-existed! {0}".format(os.path.basename(filename))
        )
    else:
        if version == "shapeit2":
            #### version1: with shapeit phasing
            shapeit2 = pd.read_csv(shapeit2_input, sep="\t", header=0, index_col=False)
            hetSNP_shapeit2 = hetSNP_intersect_unique_filtered.merge(
                shapeit2, how="inner", on=["chr", "pos"]
            )
            logging.debug(
                "size of filtered hetSNP unique is {0} ; size of shapeit2 phased variants is {1}; innerjoin size is {2}".format(
                    len(hetSNP_intersect_unique_filtered),
                    len(shapeit2),
                    len(hetSNP_shapeit2),
                )
            )
            phasing_data = hetSNP_shapeit2
        #### version2: without shapeit phasing
        else:
            hetSNP_intersect_unique_filtered[
                ["e_paternal", "e_maternal"]
            ] = hetSNP_intersect_unique_filtered.genotype.str.split("|", expand=True)
            phasing_data = hetSNP_intersect_unique_filtered
        phasing_data = phasing_data[
            [
                "chr",
                "chrN",
                "pos",
                "rsid",
                "AF",
                "geneID",
                "genotype",
                "refCount",
                "altCount",
                "totalCount",
                "altRatio",
                "FishTest",
                "alt_binomial_p",
                "e_paternal",
                "e_maternal",
            ]
        ]
        geneIDs = phasing_data["geneID"].unique()
        new_df = pd.DataFrame()
        for gene in geneIDs:
            selected_gene = phasing_data[phasing_data["geneID"] == gene]
            selected_gene = selected_gene.reset_index(drop=True)
            edited_selected_gene = change_phasing(selected_gene)
            new_df = new_df.append(edited_selected_gene)
        new_df.to_csv(filename, sep="\t", header=True, index=False)
        ###################
        # checking phasing difference within a gene
        if version == "shapeit2" and phase_difference_file is not None:
            new_df_diff = new_df.rename(
                columns={
                    "e_paternal": "shapeit2_paternal",
                    "e_maternal": "shapeit2_maternal",
                    "patCount": "shapeit2_patCount",
                    "matCount": "shapeit2_matCount",
                }
            )
            new_df_diff[["e_paternal", "e_maternal"]] = new_df_diff.genotype.str.split(
                "|", expand=True
            )
            new_df_phasingDiff = pd.DataFrame()
            for gene in geneIDs:
                selected_gene = new_df_diff[new_df_diff["geneID"] == gene]
                selected_gene = selected_gene.reset_index(drop=True)
                edited_selected_gene = change_phasing(selected_gene)
                new_df_phasingDiff = new_df_phasingDiff.append(edited_selected_gene)
            new_df_phasingDiff = new_df_phasingDiff.rename(
                columns={
                    "e_paternal": "vcf_paternal",
                    "e_maternal": "vcf_maternal",
                    "patCount": "vcf_patCount",
                    "matCount": "vcf_matCount",
                }
            )
            new_df_phasingDiff["diff"] = np.where(
                (
                    new_df_phasingDiff["shapeit2_patCount"]
                    == new_df_phasingDiff["vcf_patCount"]
                )
                & (
                    new_df_phasingDiff["shapeit2_matCount"]
                    == new_df_phasingDiff["vcf_matCount"]
                ),
                0,
                1,
            )
            new_df_phasingDiff.to_csv(
                phase_difference_file, sep="\t", header=True, index=False
            )
        ###################
        # drop NAN
        # new_df_dropNA = new_df.dropna()
        # print(new_df_dropNA)
        # new_df = new_df.rename(columns={"SNP_id": "rsid"})
        new_df_dropNA = new_df.drop(
            ["refCount", "altCount", "e_paternal", "e_maternal"], axis=1
        )
        logging.debug(
            "size of filtered hetSNP unique is {0} ; phased data size is {1}, phased data after cleaning size is {2}".format(
                len(hetSNP_intersect_unique_filtered),
                len(new_df),
                len(new_df_dropNA),
            )
        )
        new_df_dropNA.to_csv(filename_cleaned, sep="\t", header=True, index=False)


def add_simulationData(sim_filename):
    simulator_df = pd.read_csv(sim_filename, sep="\t", header=0, index_col=False)
    simulator_df["alt_binomial_p"] = simulator_df.apply(
        lambda x: binom_test(
            x["altCount"],
            x["totalCount"],
            p=0.5,
            alternative="two-sided",
        ),
        axis=1,
    )
    simulator_df.to_csv(sim_filename, index=False, sep="\t", header=True)
    simulator_df = simulator_df[["chr", "pos", "alt_binomial_p"]]
    return simulator_df


def generate_modelCount(phased_filename):
    base_out = os.path.splitext(phased_filename)[0]
    data = pd.read_csv(phased_filename, sep="\t", header=0, index_col=False)
    data = data[
        [
            "chr",
            "chrN",
            "pos",
            "patCount",
            "matCount",
            "totalCount",
            "altRatio",
            "rsid",
            "AF",
            "geneID",
            "genotype",
        ]
    ]
    file_for_lambda = "{0}.forlambda.tsv".format(base_out)
    lambdaPredicted_file = "{0}.lambdaPredicted.tsv".format(base_out)
    out_modelinput = "{0}.modelinput.tsv".format(base_out)
    out_modelinput_error = "{0}.modelinput.w_error.tsv".format(base_out)

    df_med = (
        data.groupby(["geneID"])
        .agg({"altRatio": "median"})
        .reset_index()
        .rename(columns={"altRatio": "median.altRatio"})
    )

    df_nhets = (
        data.groupby(by="geneID")
        .agg({"pos": pd.Series.nunique})
        .rename(columns={"pos": "number.of.hets"})
    )

    df_totalref = (
        data.groupby(by="geneID")
        .agg({"patCount": "sum"})
        .reset_index()
        .rename(columns={"patCount": "total.patCount"})
    )

    df_totalalt = (
        data.groupby(by="geneID")
        .agg({"matCount": "sum"})
        .reset_index()
        .rename(columns={"matCount": "total.matCount"})
    )
    df_summary_1 = pd.merge(df_med, df_nhets, on=["geneID"], how="inner")
    df_summary_2 = pd.merge(df_summary_1, df_totalref, on=["geneID"], how="inner")
    df_summary_3 = pd.merge(df_summary_2, df_totalalt, on=["geneID"], how="inner")
    df_summary_3["totalCount"] = (
        df_summary_3["total.patCount"] + df_summary_3["total.matCount"]
    )

    df_summary_3 = df_summary_3.drop(["total.patCount", "total.matCount"], axis=1)
    df_summary_3.to_csv(file_for_lambda, index=False, sep="\t", header=True)

    counter = 0
    if not os.path.isfile(out_modelinput):
        unique_gene = data["geneID"].unique()
        rst = ""
        for each in unique_gene:
            counter += 1
            idx = each
            gene_lst = [idx, sum(data["geneID"] == each)] + list(
                np.ravel(data[data["geneID"] == each][["matCount", "patCount"]].values)
            )
            rst += "\t".join([str(x) for x in gene_lst]) + "\n"
        file1 = open(out_modelinput, "w")
        file1.write(rst)
        file1.close()
        logging.debug(
            "output {0} has {1} genes as input for stan model".format(
                os.path.basename(out_modelinput), counter
            )
        )
    else:
        logging.info(
            "....... {0} exists at {1}".format(
                os.path.basename(out_modelinput),
                os.path.dirname(out_modelinput),
            )
        )
    return (
        file_for_lambda,
        lambdaPredicted_file,
        out_modelinput,
        out_modelinput_error,
    )


def count_reads(fields):
    if len(fields) >= 4:
        base_thetas = []
        Mreps = int(fields[1])
        total_A = 0
        total_R = 0
        for rep in range(Mreps):
            A = float(fields[2 + rep * 2])
            R = float(fields[3 + (rep) * 2])
            base = (A + 1) / (R + 1)
            base_thetas.append(base)
            total_A += A
            total_R += R
        total_reads = total_A + total_R
    return total_reads


def update_model_input_lambda_phasing(
    pred_prob_column, base_modelin, base_modelin_error, meta_error
):
    pred_prob_column = "pred_error_GIAB"
    outfile = base_modelin_error
    model_input = base_modelin
    phasing_error = pd.read_csv(meta_error, header=0, sep="\t")
    updated_line = ""
    counter = 0

    with open(model_input, "r") as stream_in:
        for i, line in enumerate(
            stream_in
        ):  # start reading in my pileup results line by line
            line = re.sub("\n", "", line)
            # logging.info("{0}".format(line))
            counter += 1
            if counter % 1000 == 0:
                logging.info("{0} line processed".format(counter))
            # logging.debug(line)
            geneID = line.split("\t")[0]
            # logging.debug("geneID :{0}".format(geneID))
            nhet = int(line.split("\t")[1])
            # logging.debug("num of hets :{0}".format(nhet))
            phasing_error_selected = phasing_error[(phasing_error["geneID"] == geneID)]
            # logging.debug(phasing_error_selected)
            pred_prob = phasing_error_selected[str(pred_prob_column)].tolist()

            if len(pred_prob) != 0:
                del pred_prob[0]
            # logging.debug("original predicted prob :{0}".format(pred_prob))
            # count number of NAs in the predicted switching error for each gene
            nNAs = np.count_nonzero(np.isnan(pred_prob))
            # logging.debug("num of NAs :{0}".format(nNAs))
            # replace NA with -1 for predicted switching error
            pred_prob = [-1 if x != x else x for x in pred_prob]
            # logging.debug("corrected predicted prob :{0}".format(pred_prob))
            PI_pred = []
            # for gene with 1 het, number of NA is 1, put 0 there
            if nhet == 1:
                updated_line += "%s\t%d\n" % (line, 0)
            # for gene with more than 1 hets, append predicted switching error after nNAs
            else:
                for m in range(int(nhet) - 1):  # 0,1,2
                    if m == 0:  # the first het site
                        PI_pred = "%d\t%s" % (nNAs, pred_prob[m])
                    else:
                        PI_pred = "%s\t%s" % (PI_pred, pred_prob[m])
                updated_line += "%s\t%s\n" % (line, PI_pred)
            # if counter == 20:
            #     sys.exit()
    file1 = open(outfile, "w")
    file1.write(updated_line)
    file1.close()


def significant_genes(
    df_ibeastie,
    df_binomial,
    df_adm,
    outfilename,
    outfilename_ase,
    cutoff,
    hetSNP_intersect_unique_lambdaPredicted_file,
):
    data_modeloutput = pd.read_csv(
        hetSNP_intersect_unique_lambdaPredicted_file, header=0, sep="\t"
    )
    logging.debug(
        "size of hetSNP_intersect_unique_lambdaPredicted_file file {}".format(
            len(data_modeloutput)
        )
    )
    logging.debug("size of df_beastie {}".format(len(df_ibeastie)))

    df_output = pd.merge(data_modeloutput, df_ibeastie, on=["geneID"], how="inner")
    logging.debug(
        "size of model output is {0} ; size of hetSNP_intersect_unique_lambdaPredicted_file is {1}; intersection size is {2}".format(
            len(df_ibeastie), len(data_modeloutput), len(df_output)
        )
    )
    ncount = df_output[df_output["posterior_mass_support_ALT"] > cutoff].count()[8]
    logging.info(
        "{} genes with ASE out of total genes {} ({}%) at @ {} > ASE cutoff {}".format(
            ncount,
            len(df_output),
            round((ncount / len(df_output)) * 100, 3),
            "posterior_mass_support_ALT",
            0.5,
        )
    )
    df_output["foldLog2MedSq_over_std"] = (
        abs(np.log2(df_output["posterior_median"]))
    ) ** 2 / abs(np.log2(df_output["posterior_variance"]))
    # df_output["foldLog2MedSq_over_var"]=df_output["foldLog2MedSq_over_var"].apply(lambda x: round(x, 3 - int(floor(log10(abs(x))))))
    df_output["median.altRatio"] = df_output["median.altRatio"].round(3)
    # df_output["extreme_val"] = df_output["foldLog2MedSq_over_var"].apply(lambda x: round(x, 3 - int(floor(log10(abs(x))))))
    extremes = scipy.stats.chi2.ppf(0.95, 1)
    ncount2 = df_output[df_output["foldLog2MedSq_over_std"] > extremes].count()[8]
    logging.debug(
        "{} genes with ASE out of total genes {} ({}%) at @ {} > extreme values: {} = QCHISQ(p=0.95,df=1)".format(
            ncount2,
            len(df_output),
            round((ncount2 / len(data_modeloutput)) * 100, 3),
            "foldLog2MedSq_over_std",
            extremes,
        )
    )
    df_output = df_output.assign(
        extreme_val=lambda dataframe: dataframe["foldLog2MedSq_over_std"].map(
            lambda foldLog2MedSq_over_std: "Y"
            if foldLog2MedSq_over_std > extremes
            else "N"
        )
    )
    df_output["ASE"] = df_output["posterior_mass_support_ALT"]
    # print(">> merged V2")
    # print(df_output.head(n=5))
    df_output = df_output.assign(
        ASE=lambda dataframe: dataframe["posterior_mass_support_ALT"].map(
            lambda posterior_mass_support_ALT: "Y"
            if posterior_mass_support_ALT > 0.5
            else "N"
        )
    )
    df_output_bi = pd.merge(df_output, df_binomial, on=["geneID"], how="inner")

    if df_output_bi.shape[0] == 0:
        print(">> merged binomial with beastie is empty")
        df_output_bi = df_output
    else:
        df_output_bi = df_output_bi.drop(
            [
                "FirstSite_pval",
                "FirstSite_esti",
                "NaiveSum_esti",
                "Pseudo_pval",
                "Pseudo_esti",
                "MajorSite_esti",
            ],
            axis=1,
        )
    df_output_bi_adm = pd.merge(df_output_bi, df_adm, on=["geneID"], how="inner")
    # print(">> merged binomial with beastie and ADM")
    if df_output_bi_adm.shape[0] == 0:
        print(">> merged binomial with beastie and ADM is empty")
        df_output_bi_adm = df_output
    else:
        df_output_bi_adm = df_output_bi_adm.drop(["ADM_esti"], axis=1)
    # print(df_output_bi_adm.head(n=5))

    logging.info("....... done with dropping")
    df_output = df_output_bi_adm
    # outfilename="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/other_example/HG00096/output/s-0.5_a-0.05_sinCov0_totCov1_W1000K1000/HG00096_ASE_all.tsv"
    # outfilename_ase="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/other_example/HG00096/output/s-0.5_a-0.05_sinCov0_totCov1_W1000K1000/HG00096_ASE_cutoff_0.5_filtered.tsv"
    df_output_bi.to_csv(outfilename, sep="\t", header=True, index=False)
    df_output_ase = df_output[df_output["ASE"] == "Y"]
    df_output_ase.to_csv(outfilename_ase, sep="\t", header=True, index=False)
