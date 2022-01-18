#!/usr/bin/env python
# =========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
# =========================================================================

import re
import os
import sys
import logging
import scipy.stats
import pandas as pd
import numpy as np
from math import floor, log10, pi, isnan


def change_phasing(data):
    row_n = data.shape[0]
    data["patCount"] = data["refCount"]
    data["matCount"] = data["altCount"]
    reference_row = data.iloc[0]
    if not isnan(reference_row["e_paternal"]):
        for x in range(1, row_n):
            row = data.iloc[x]
            if row["e_paternal"] != reference_row["e_paternal"] and not isnan(
                row["e_paternal"]
            ):
                data.at[x, "patCount"] = row["altCount"]
                data.at[x, "matCount"] = row["refCount"]
    return data


def incorporate_shapeit2(
    shapeit2,
    hetSNP_intersect_unique,
    filename,
    version="V2",
):
    hetSNP_shapeit2 = pd.merge(
        hetSNP_intersect_unique, shapeit2, how="left", on=["chr", "pos"]
    )
    hetSNP_shapeit2 = hetSNP_shapeit2[
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
            "e_paternal",
            "e_maternal",
        ]
    ]
    base_out = os.path.splitext(filename)[0]
    hetSNP_intersect_unique_shapeit2 = f"{base_out}.shapeit2.tsv"
    hetSNP_intersect_unique_shapeit2_dropNA = f"{base_out}.shapeit2.dropNA.tsv"
    #### version1
    if version == "V1":
        geneIDs = hetSNP_shapeit2["geneID"].unique()
        new_df = pd.DataFrame()
        for gene in geneIDs:
            selected_gene = hetSNP_shapeit2[hetSNP_shapeit2["geneID"] == gene]
            selected_gene = selected_gene.reset_index(drop=True)
            edited_selected_gene = change_phasing(selected_gene)
            new_df = new_df.append(edited_selected_gene)

    #### version2
    else:
        new_df = hetSNP_shapeit2
        row_n = new_df.shape[0]
        s = new_df["e_paternal"] == 1
        new_df["patCount"] = new_df["refCount"]
        new_df["matCount"] = new_df["altCount"]
        new_df.loc[s, ["patCount", "matCount"]] = new_df.loc[
            s, ["altCount", "refCount"]
        ].values
    new_df.to_csv(hetSNP_intersect_unique_shapeit2, sep="\t", header=True, index=False)
    # drop NAN
    new_df_dropNA = new_df.dropna()
    new_df_dropNA = new_df_dropNA.drop(
        ["refCount", "altCount", "e_paternal", "e_maternal"], axis=1
    )
    logging.debug(
        "size of hetSNP unique is {0} ; size of shapeit2 phased variants is {1}; leftjoin size is {2}; leftjoin dropNA size is {3}".format(
            len(hetSNP_intersect_unique), len(shapeit2), len(new_df), len(new_df_dropNA)
        )
    )
    new_df_dropNA.to_csv(
        hetSNP_intersect_unique_shapeit2_dropNA, sep="\t", header=True, index=False
    )
    # return new_df_dropNA


def add_shapepit2(
    filename,
    shapeit2_input,
):
    hetSNP_intersect_unique_shapeit2 = f"{os.path.splitext(filename)[0]}.shapeit2.tsv"
    hetSNP_intersect_unique_shapeit2_dropNA = (
        f"{os.path.splitext(filename)[0]}.shapeit2.dropNA.tsv"
    )
    if not os.path.isfile(hetSNP_intersect_unique_shapeit2_dropNA):
        gene_df = pd.read_csv(filename, sep="\t", header=0, index_col=False)
        shapeit2_df = pd.read_csv(shapeit2_input, sep="\t", header=0, index_col=False)
        incorporate_shapeit2(shapeit2_df, gene_df, filename)

        logging.info(
            "output {0} has hetSNP data with shapeit2 phasing".format(
                os.path.basename(hetSNP_intersect_unique_shapeit2_dropNA)
            )
        )
    else:
        logging.info(
            "....... {0} existed has hetSNP data with shapeit2 phasing".format(
                os.path.basename(hetSNP_intersect_unique_shapeit2_dropNA)
            )
        )


def generate_modelCount(filename, shapeit2_input):
    if shapeit2_input is None:
        gene_df = pd.read_csv(filename, sep="\t", header=0, index_col=False)
        data = gene_df
        data["patCount"] = data["refCount"]
        data["matCount"] = data["altCount"]
        data = data[
            [
                "chr",
                "chrN",
                "pos",
                "rsid",
                "AF",
                "geneID",
                "genotype",
                "patCount",
                "matCount",
                "totalCount",
                "altRatio",
            ]
        ]
    else:
        base_out = os.path.splitext(filename)[0]
        filename = f"{base_out}.shapeit2.dropNA.tsv"
        gene_df = pd.read_csv(filename, sep="\t", header=0, index_col=False)
        data = gene_df
    file_for_LDannotation = filename
    logging.debug(
        "input {0} has {1} het SNPs".format(
            os.path.basename(filename), gene_df.shape[0]
        )
    )
    base_out = os.path.splitext(filename)[0]
    file_for_lambda = "{0}.forlambda.tsv".format(base_out)
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
        unique_gene = gene_df["geneID"].unique()
        rst = ""
        for each in unique_gene:
            counter += 1
            idx = each
            gene_lst = [idx, sum(gene_df["geneID"] == each)] + list(
                np.ravel(
                    gene_df[gene_df["geneID"] == each][["matCount", "patCount"]].values
                )
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
    return file_for_LDannotation, file_for_lambda, out_modelinput, out_modelinput_error


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
            nNAs = np.count_nonzero(np.isnan(pred_prob))
            # logging.debug("num of NAs :{0}".format(nNAs))
            pred_prob = [-1 if x != x else x for x in pred_prob]  # -1 if missing
            # logging.debug("corrected predicted prob :{0}".format(pred_prob))
            PI_pred = []
            if nhet == 1:
                updated_line += "%s\t%d\n" % (line, 0)  # 0 for 1 het gene
            else:
                for m in range(int(nhet) - 1):  # 0,1,2
                    if m == 0:
                        PI_pred = "%d\t%s" % (nNAs, pred_prob[m])
                    else:
                        PI_pred = "%s\t%s" % (PI_pred, pred_prob[m])
                updated_line += "%s\t%s\n" % (line, PI_pred)
            # if counter == 10:
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
            cutoff,
        )
    )
    df_output["foldLog2MedSq_over_var"] = (
        abs(np.log2(df_output["posterior_median"]))
    ) ** 2 / df_output["posterior_variance"]
    # df_output["foldLog2MedSq_over_var"]=df_output["foldLog2MedSq_over_var"].apply(lambda x: round(x, 3 - int(floor(log10(abs(x))))))
    df_output["median.altRatio"] = df_output["median.altRatio"].round(3)
    # df_output["extreme_val"] = df_output["foldLog2MedSq_over_var"].apply(lambda x: round(x, 3 - int(floor(log10(abs(x))))))
    extremes = scipy.stats.chi2.ppf(0.95, 1)
    ncount2 = df_output[df_output["foldLog2MedSq_over_var"] > extremes].count()[8]
    logging.debug(
        "{} genes with ASE out of total genes {} ({}%) at @ {} > extreme values: {} = QCHISQ(p=0.95,df=1)".format(
            ncount2,
            len(df_output),
            round((ncount2 / len(data_modeloutput)) * 100, 3),
            "foldLog2MedSq_over_var",
            extremes,
        )
    )
    df_output = df_output.assign(
        extreme_val=lambda dataframe: dataframe["foldLog2MedSq_over_var"].map(
            lambda foldLog2MedSq_over_var: "Y"
            if foldLog2MedSq_over_var > extremes
            else "N"
        )
    )
    df_output["ASE"] = df_output["posterior_mass_support_ALT"]
    # print(">> merged V2")
    # print(df_output.head(n=5))
    df_output = df_output.assign(
        ASE=lambda dataframe: dataframe["posterior_mass_support_ALT"].map(
            lambda posterior_mass_support_ALT: "Y"
            if posterior_mass_support_ALT > cutoff
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
