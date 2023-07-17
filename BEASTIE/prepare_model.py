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
from scipy.stats import binom_test


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


def filter_alignBias(
    prefix,
    tmp_path,
    binomialp_cutoff,
    genotypeErfiltered_hetSNP_intersect_pileup,
    simulator_df=None,
    collected_alignmentBias_file=None,
):
    new_df = pd.read_csv(
        genotypeErfiltered_hetSNP_intersect_pileup, sep="\t", header=0, index_col=False
    )

    ################# 1. w/ alignment bias SNP file --> filter out variants with alignment bias from given SNP list
    if collected_alignmentBias_file is not None:
        SNPs_needtobe_filtered = pd.read_csv(
            collected_alignmentBias_file, sep="\t", header=None, index_col=False
        )
        SNPs_needtobe_filtered.columns = ["chrN_pos"]
        SNP_list = SNPs_needtobe_filtered.chrN_pos.values.tolist()
        new_df["chrN_pos"] = new_df.apply(
            lambda x: "%s_%s" % (x["chrN"], x["pos"]), axis=1
        )
        # SNP_list2 = ["9_140483339", "9_140507688"]
        new_df_SNPlistfiltered = new_df[~new_df["chrN_pos"].isin(SNP_list)]
        new_df_SNPlistfiltered = new_df_SNPlistfiltered[
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
                "genotypeTest",
            ]
        ]
        afterFilter_filename = os.path.join(
            tmp_path, f"{prefix}_givenBiasSNPlistFiltered.tsv"
        )
        new_df_SNPlistfiltered.to_csv(
            afterFilter_filename, index=False, sep="\t", header=True
        )
        filename = afterFilter_filename
        logging.debug(
            "{0} has {1} het SNPs after genotyping error filtering. After given alignment Bias SNP list filtering {2} has {3} het SNPs".format(
                os.path.basename(genotypeErfiltered_hetSNP_intersect_pileup),
                new_df.shape[0],
                os.path.basename(afterFilter_filename),
                new_df_SNPlistfiltered.shape[0],
            )
        )
    ################# 2. w/ simulation data --> filter out variants with alignment bias from simulation data
    elif simulator_df is not None:
        beforeFilter_filename = os.path.join(
            tmp_path, f"{prefix}_real_alignBiasbeforeFilter.tsv"
        )
        afterFilter_filename = os.path.join(
            tmp_path, f"{prefix}_real_alignBiasafterFilter.tsv"
        )
        # alignment bias binomial test score
        overlapped_variants = new_df.merge(
            simulator_df, on=["chrN", "pos"], how="inner"
        )
        # saving
        overlapped_variants.to_csv(
            beforeFilter_filename, index=False, sep="\t", header=True
        )
        Notoverlapped_variants = new_df.drop_duplicates().merge(
            simulator_df.drop_duplicates(),
            on=["chrN", "pos"],
            how="left",
            indicator=True,
        )
        Notoverlapped_variants = Notoverlapped_variants[
            Notoverlapped_variants["_merge"] == "left_only"
        ]
        Notoverlapped_variants_filename = os.path.join(
            tmp_path, f"{prefix}_real_notoverlapped.tsv"
        )
        Notoverlapped_variants.to_csv(
            Notoverlapped_variants_filename, index=False, sep="\t", header=True
        )
        logging.debug(
            "{0} het SNPs in real data after genotyping error filtering, {1} het SNPs in simulation data, overlapped het SNPs {2} in real data, not overlapped {3}".format(
                new_df.shape[0],
                simulator_df.shape[0],
                overlapped_variants.shape[0],
                Notoverlapped_variants.shape[0],
            )
        )
        no_alignment_bias = overlapped_variants[
            overlapped_variants["alt_binomial_p"] > binomialp_cutoff
        ]

        logging.debug(
            "{0} ({1}%) overlapped het SNPs pass alignment bias filtering p-val > {2}".format(
                no_alignment_bias.shape[0],
                round(
                    no_alignment_bias.shape[0] / overlapped_variants.shape[0] * 100,
                    10,
                ),
                binomialp_cutoff,
            )
        )
        gene_df_filtered = no_alignment_bias[
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
                "genotypeTest",
            ]
        ]
        gene_df_filtered.to_csv(
            afterFilter_filename, index=False, sep="\t", header=True
        )
        filename = afterFilter_filename

    ################# 3. no simulation data --> do not filter out variants with alignment bias
    else:
        hetSNP_intersect_unique_subset = new_df[
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
                "genotypeTest",
            ]
        ]
        afterFilter_filename = os.path.join(
            tmp_path, f"{prefix}_NoalignBiasFiltered.tsv"
        )
        hetSNP_intersect_unique_subset.to_csv(
            afterFilter_filename, index=False, sep="\t", header=True
        )
        filename = afterFilter_filename
        logging.debug(
            "{0} has {1} het SNPs, after genotyping error filtering, {2} has {3} het SNPs, no alignment bias filtering".format(
                os.path.basename(genotypeErfiltered_hetSNP_intersect_pileup),
                new_df.shape[0],
                os.path.basename(afterFilter_filename),
                hetSNP_intersect_unique_subset.shape[0],
            )
        )
    return filename


def re_allocateReads(
    alignBiasfiltered_filename,
    shapeit2_input,
    version,
    filename,
    filename_cleaned,
    collected_alignmentBias_file=None,
    simulation_pileup=None,
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
        # use everything can be phased by shapeit2, no matter whether VCF have phased genotype
        if version != "nophasing":
            if version == "shapeit2":
                #### version1: with shapeit phasing
                shapeit2 = pd.read_csv(
                    shapeit2_input,
                    sep="\t",
                    header=0,
                    index_col=False,
                )
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
                #### version2: VCF phasing
            elif version == "VCF":
                hetSNP_intersect_unique_filtered = hetSNP_intersect_unique_filtered[
                    (hetSNP_intersect_unique_filtered["genotype"] != "1/0")
                    & (hetSNP_intersect_unique_filtered["genotype"] != "0/1")
                ]
                hetSNP_intersect_unique_filtered[
                    ["e_paternal", "e_maternal"]
                ] = hetSNP_intersect_unique_filtered.genotype.str.split(
                    "|", expand=True
                )
                phasing_data = hetSNP_intersect_unique_filtered
            ### cleaning
            geneIDs = phasing_data["geneID"].unique()
            new_df = pd.DataFrame()
            for gene in geneIDs:
                selected_gene = phasing_data[phasing_data["geneID"] == gene]
                selected_gene = selected_gene.reset_index(drop=True)
                edited_selected_gene = change_phasing(selected_gene)
                new_df = pd.concat([new_df, edited_selected_gene], ignore_index=True)
        #### version3: no phasing
        else:
            new_df = hetSNP_intersect_unique_filtered
            new_df["patCount"] = new_df["refCount"]
            new_df["matCount"] = new_df["altCount"]
        new_df.to_csv(filename, sep="\t", header=True, index=False)
        # collected alignment bias SNP list filtering
        # if (collected_alignmentBias_file is not None) or (simulation_pileup is None):
        #     phasing_data = phasing_data[
        #         [
        #             "chr",
        #             "chrN",
        #             "pos",
        #             "rsid",
        #             "AF",
        #             "geneID",
        #             "genotype",
        #             "refCount",
        #             "altCount",
        #             "totalCount",
        #             "altRatio",
        #             "e_paternal",
        #             "e_maternal",
        #         ]
        #     ]
        # elif simulation_pileup is not None:
        #     phasing_data = phasing_data[
        #         [
        #             "chr",
        #             "chrN",
        #             "pos",
        #             "rsid",
        #             "AF",
        #             "geneID",
        #             "genotype",
        #             "refCount",
        #             "altCount",
        #             "totalCount",
        #             "altRatio",
        #             "genotypeTest",
        #             "alt_binomial_p",
        #             "e_paternal",
        #             "e_maternal",
        #         ]
        #     ]
        ###################
        # checking phasing difference within a gene (for sites with VCF phased genotype and shapeit2 phasing)
        if version == "shapeit2" and phase_difference_file is not None:
            new_df_filtered = new_df[
                (new_df["genotype"] != "1/0") & (new_df["genotype"] != "0/1")
            ]
            logging.debug(
                "size of shapeit2 phased sites inner join vcf het sites is {0} ; overlapped with only VCF phased sites is {1}".format(
                    len(new_df), len(new_df_filtered)
                )
            )
            if len(new_df) == len(new_df_filtered):
                new_df_diff = new_df_filtered.rename(
                    columns={
                        "e_paternal": "shapeit2_paternal",
                        "e_maternal": "shapeit2_maternal",
                        "patCount": "shapeit2_patCount",
                        "matCount": "shapeit2_matCount",
                    }
                )
                new_df_diff[
                    ["e_paternal", "e_maternal"]
                ] = new_df_diff.genotype.str.split("|", expand=True)
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
            else:
                logging.debug(
                    "there are sites without VCF phasing, skip phasing difference file generation!"
                )
        ###################
        # drop NAN
        # print(new_df)
        new_df_dropNA = new_df
        # new_df_dropNA = new_df.dropna()
        # new_df = new_df.rename(columns={"SNP_id": "rsid"})
        # new_df_dropNA = new_df.drop(
        #     ["refCount", "altCount", "e_paternal", "e_maternal"], axis=1
        # )
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
    simulator_df["chrN"] = simulator_df["contig"]
    simulator_df["pos"] = simulator_df["position"]
    simulator_df = simulator_df[["chrN", "pos", "alt_binomial_p"]]
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

    # df_summary_3 = df_summary_3.drop(["total.patCount", "total.matCount"], axis=1)
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
    df_beastie,
    df_binomial,
    df_adm,
    outfilename,
    outfilename_sub,
    outfilename_ase,
    ase_cutoff,
    hetSNP_intersect_unique_lambdaPredicted_file,
    adjusted_alpha,
):
    data_modeloutput = pd.read_csv(
        hetSNP_intersect_unique_lambdaPredicted_file, header=0, sep="\t"
    )
    logging.debug(
        "size of hetSNP_intersect_unique_lambdaPredicted_file file {}".format(
            len(data_modeloutput)
        )
    )
    logging.debug("size of df_beastie {}".format(len(df_beastie)))

    df_output = pd.merge(data_modeloutput, df_beastie, on=["geneID"], how="inner")
    logging.debug(
        "size of model output is {0} ; size of hetSNP_intersect_unique_lambdaPredicted_file is {1}; intersection size is {2}".format(
            len(df_beastie), len(data_modeloutput), len(df_output)
        )
    )

    # df_output["foldLog2MedSq_over_std"] = (
    #     abs(np.log2(df_output["posterior_median"]))
    # ) ** 2 / abs(np.log2(df_output["posterior_variance"]))
    # df_output["foldLog2MedSq_over_var"]=df_output["foldLog2MedSq_over_var"].apply(lambda x: round(x, 3 - int(floor(log10(abs(x))))))
    df_output["median.altRatio"] = df_output["median.altRatio"].round(3)
    # df_output["extreme_val"] = df_output["foldLog2MedSq_over_var"].apply(lambda x: round(x, 3 - int(floor(log10(abs(x))))))
    extremes = scipy.stats.chi2.ppf(0.95, 1)
    # ncount2 = df_output[df_output["foldLog2MedSq_over_std"] > extremes].count()[8]
    # logging.debug(
    #     "{} genes with ASE out of total genes {} ({}%) at @ {} > extreme values: {} = QCHISQ(p=0.95,df=1)".format(
    #         ncount2,
    #         len(df_output),
    #         round((ncount2 / len(data_modeloutput)) * 100, 3),
    #         "foldLog2MedSq_over_std",
    #         extremes,
    #     )
    # )
    # df_output = df_output.assign(
    #     extreme_val=lambda dataframe: dataframe["foldLog2MedSq_over_std"].map(
    #         lambda foldLog2MedSq_over_std: "Y"
    #         if foldLog2MedSq_over_std > extremes
    #         else "N"
    #     )
    # )
    # df_output["ASE"] = df_output["posterior_mass_support_ALT"]
    # output that contain everything
    # df_output = df_output.assign(
    #     ASE=lambda dataframe: dataframe["posterior_mass_support_ALT"].map(
    #         lambda posterior_mass_support_ALT: "Y"
    #         if posterior_mass_support_ALT > ase_cutoff
    #         else "N"
    #     )
    # )

    # subset output
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
    # df_output = df_output_bi_adm
    df_output = df_output_bi

    def beastie_gam(row):
        if float(row["posterior_mass_support_ALT_gam"]) > ase_cutoff:
            return 1
        else:
            return 0

    def NS(row, adjusted_alpha):
        if float(row["NaiveSum_pval"]) <= float(adjusted_alpha):
            return 1
        else:
            return 0

    def MS(row, adjusted_alpha):
        if float(row["MajorSite_pval"]) <= float(adjusted_alpha):
            return 1
        else:
            return 0
    def beta11(row, adjusted_alpha):
        if float(row["beta_1_1_pval"]) <= float(adjusted_alpha):
            return 1
        else:
            return 0
    def beta1010(row, adjusted_alpha):
        if float(row["beta_10_10_pval"]) <= float(adjusted_alpha):
            return 1
        else:
            return 0
    def beta2020(row, adjusted_alpha):
        if float(row["beta_20_20_pval"]) <= float(adjusted_alpha):
            return 1
        else:
            return 0
    def beta5050(row, adjusted_alpha):
        if float(row["beta_50_50_pval"]) <= float(adjusted_alpha):
            return 1
        else:
            return 0
    def beta100100(row, adjusted_alpha):
        if float(row["beta_100_100_pval"]) <= float(adjusted_alpha):
            return 1
        else:
            return 0

    df_output["beastie_ASE_gam"] = df_output.apply(
        lambda row: beastie_gam(row), axis=1
    )
    df_output["NS_ASE"] = df_output.apply(lambda row: NS(row, adjusted_alpha), axis=1)
    df_output["MS_ASE"] = df_output.apply(lambda row: MS(row, adjusted_alpha), axis=1)
    df_output["beta11_ASE"] = df_output.apply(lambda row: beta11(row, adjusted_alpha), axis=1)
    df_output["beta1010_ASE"] = df_output.apply(lambda row: beta1010(row, adjusted_alpha), axis=1)
    df_output["beta2020_ASE"] = df_output.apply(lambda row: beta2020(row, adjusted_alpha), axis=1)
    df_output["beta5050_ASE"] = df_output.apply(lambda row: beta5050(row, adjusted_alpha), axis=1)
    df_output["beta100100_ASE"] = df_output.apply(lambda row: beta100100(row, adjusted_alpha), axis=1)

    ncount1 = df_output["beastie_ASE_gam"].sum()
    logging.info(
        "{} genes with ASE out of total genes {} ({}%) at @ {} > ASE cutoff {}".format(
            ncount1,
            len(df_output),
            round((ncount1 / len(df_output)) * 100, 3),
            "posterior_mass_support_ALT_gam",
            ase_cutoff,
        )
    )
    ncount2 = df_output["NS_ASE"].sum()
    logging.info(
        "{} genes with ASE out of total genes {} ({}%) at @ {} <= adjusted alpha {}".format(
            ncount2,
            len(df_output),
            round((ncount2 / len(df_output)) * 100, 3),
            "Naive Sum",
            adjusted_alpha,
        )
    )
    ncount3 = df_output["MS_ASE"].sum()
    logging.info(
        "{} genes with ASE out of total genes {} ({}%) at @ {} <= adjusted alpha {}".format(
            ncount3,
            len(df_output),
            round((ncount3 / len(df_output)) * 100, 3),
            "Major Site",
            adjusted_alpha,
        )
    )
    ncount4 = df_output["beta11_ASE"].sum()
    logging.info(
        "{} genes with ASE out of total genes {} ({}%) at @ {} <= adjusted alpha {}".format(
            ncount4,
            len(df_output),
            round((ncount4 / len(df_output)) * 100, 3),
            "beta_1_1",
            adjusted_alpha,
        )
    )
    ncount5 = df_output["beta1010_ASE"].sum()
    logging.info(
        "{} genes with ASE out of total genes {} ({}%) at @ {} <= adjusted alpha {}".format(
            ncount5,
            len(df_output),
            round((ncount5 / len(df_output)) * 100, 3),
            "beta_10_10",
            adjusted_alpha,
        )
    )
    ncount6 = df_output["beta2020_ASE"].sum()
    logging.info(
        "{} genes with ASE out of total genes {} ({}%) at @ {} <= adjusted alpha {}".format(
            ncount6,
            len(df_output),
            round((ncount6 / len(df_output)) * 100, 3),
            "beta_20_20",
            adjusted_alpha,
        )
    )
    ncount7 = df_output["beta5050_ASE"].sum()
    logging.info(
        "{} genes with ASE out of total genes {} ({}%) at @ {} <= adjusted alpha {}".format(
            ncount7,
            len(df_output),
            round((ncount7 / len(df_output)) * 100, 3),
            "beta_50_50",
            adjusted_alpha,
        )
    )
    ncount8 = df_output["beta100100_ASE"].sum()
    logging.info(
        "{} genes with ASE out of total genes {} ({}%) at @ {} <= adjusted alpha {}".format(
            ncount8,
            len(df_output),
            round((ncount8 / len(df_output)) * 100, 3),
            "beta_100_100",
            adjusted_alpha,
        )
    )
    df_output_sub = df_output.drop(
        [
            "median_abs_deviation",
            "CI_left",
            "CI_right",
            "abslog2_posterior_variance",
            "abslog2_posterior_mean",
            "abslog2_posterior_median",
            "log2_posterior_variance",
            "log2_posterior_mean",
            "log2_posterior_median",
            "median_abs_deviation",
            "Pseudo_pval",
        ],
        axis=1,
    )

    # outfilename="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/other_example/HG00096/output/s-0.5_a-0.05_sinCov0_totCov1_W1000K1000/HG00096_ASE_all.tsv"
    # outfilename_ase="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/other_example/HG00096/output/s-0.5_a-0.05_sinCov0_totCov1_W1000K1000/HG00096_ASE_cutoff_0.5_filtered.tsv"
    df_output.to_csv(outfilename, sep="\t", header=True, index=False)
    df_output_sub.to_csv(outfilename_sub, sep="\t", header=True, index=False)
    # df_output_ase = df_output_sub[df_output_sub["ASE"] == "Y"]
    # df_output_ase.to_csv(outfilename_ase, sep="\t", header=True, index=False)
