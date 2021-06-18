import statistics
import pandas as pd
import numpy as np


def summary_statistics(data, title):
    print(title + " statistics:")
    if isinstance(data, list):
        print(
            "  #: {:-2} \n  Max: {:-2} \n  Min: {:-2} \n  Mean: {:-2} \n  Median: {:-2} \n  Variance: {:-2} \n  Std: {:-2} \n  25% quantile: {:-2} \n  75% quantile: {:-2}  \n  75% quantile + IQR*1.5: {:-2} \n  90% quantile: {:-2}  \n  91% quantile: {:-2}  \n  95% quantile: {:-2} \n  98% quantile: {:-2} \n  99% quantile: {:-2}".format(
                len(data),
                max(data),
                min(data),
                statistics.mean(data),
                statistics.median(data),
                statistics.variance(data),
                statistics.stdev(data),
                np.quantile(data, 0.25),
                np.quantile(data, 0.75),
                np.quantile(data, 0.75)
                + (np.quantile(data, 0.75) - np.quantile(data, 0.25)) * 1.5,
                np.quantile(data, 0.90),
                np.quantile(data, 0.91),
                np.quantile(data, 0.95),
                np.quantile(data, 0.98),
                np.quantile(data, 0.99),
            )
        )
    else:
        print("Not a list")


def Intersect_exonicHetSnps(
    parsed_mpileup,
    meta_file,
    prefix,
    out,
    min_totalcounts=1,
    min_singlecounts=0,
    DEBUG=False,
    percentile=None,
):
    altRatio = []
    tCount = []
    rCount = []
    aCount = []
    medAltRatio = []
    df = parsed_mpileup
    #   contig  position  variantID  ... lowBaseQDepth  rawDepth otherCount
    #       1    935222  rs2298214  ...             0         3          0
    #       1    949608     rs1921  ...             0      3322          9
    if DEBUG is True:
        print(df.head(n=2))
    # hetsMeta = pd.read_csv(hetSNP_file, index_col=False,sep='\t')
    # chr   | chrN  |     geneID     | genomicCoord_pos | transcriptCoord |    SNP_id |  genotype
    # chr1    1       ENSG00000227232.4       14930             -1         rs75454623      1|0
    meta_data = pd.read_csv(meta_file, index_col=False, sep="\t")
    meta_data["contig"] = meta_data["chr"].str.split("chr").str[1]
    # if DEBUG is True:
    #    print(meta_data.head(n=2))
    # "chr"   "pos"   "lag_pos"       "rsid"  "exon_start"    "exon_end"      "geneID"        "error" "distance"      "log10_distance"        "r2"    "d"     "EUR_MAF"       "lag_EUR_MAF"   "min_EUR_MAF"   "diff_EUR_MAF"
    # "chr1"  "   916549"     "   914940"     "rs6660139"     "   916516"     "   916553"     "ENSG00000187642"       "0"     "  1609"        "3.2065560"     "0.503" "1.000" "0.2416" "0.4105"        "0.2416"        "0.1689"
    # hetsMeta['pos']=hetsMeta['genomicCoord_pos']
    # fullinfo_hetsMeta = pd.merge(hetsMeta, sample_info,on=['chr','pos'], how="inner") #df6['gene']=df6.apply(lambda x: FindKey(str(x['contig']),str(x['position']),data,if_list=False),axis=1)
    # outfilename=out+prefix+"_logisticRegression_input_overlapped.tsv"
    # fullinfo_hetsMeta.to_csv(outfilename,index=False,sep='\t',header=True)
    meta_data["position"] = meta_data["pos"]
    df_hetsMeta = meta_data
    # df_hetsMeta= meta_data[['contig','position','geneID','rsid']]
    #  contig  position           geneID       rsid
    #   chr1    916549  ENSG00000187642  rs6660139
    if DEBUG is True:
        print("DEBUG USE - finish loading meta data file: %s" % (meta_data))
        print(
            "DEBUG USE - finish loading parsed mpileup file: %s\nStarting intersecting with biallelic-het-SNP file"
            % (parsed_mpileup)
        )
    if len(df):
        df2 = df[df["refCount"] >= int(min_singlecounts)]
        df3 = df2[df2["altCount"] >= int(min_singlecounts)]
        df4 = df3[df3["totalCount"] >= int(min_totalcounts)]
        if percentile is not None:
            pct = np.quantile(df4["totalCount"], percentile)
            df4 = df4[df4["totalCount"] <= pct]
        df5 = df4[df4["if_SNP"] == "Y"]
        df6 = df5[df5["if_biallelic"] == "Y"]
        freq = (
            df6["altCount"] / df6["totalCount"]
        )  # df6 = pd.concat([df_hetsMeta, df6], axis=1, join="inner")
        df6["contig"] = df6["contig"].astype(int)
        df6["position"] = df6["position"].astype(int)
        df_hetsMeta["contig"] = df_hetsMeta["contig"].astype(int)
        df_hetsMeta["position"] = df_hetsMeta["position"].astype(int)
        df_overlapped = pd.merge(
            df6, df_hetsMeta, on=["contig", "position"], how="inner"
        )  # df6['gene']=df6.apply(lambda x: FindKey(str(x['contig']),str(x['position']),data,if_list=False),axis=1)
        print(df_overlapped.head(n=2))
        df6_med = (
            df_overlapped.groupby(["geneID"])
            .agg({"altRatio": "median"})
            .reset_index()
            .rename(columns={"altRatio": "median altRatio"})
        )
        med_atr = np.ndarray.tolist(df6_med["median altRatio"].values)
        df6_nhets = (
            df_overlapped.groupby(by="geneID")
            .agg({"position": pd.Series.nunique})
            .rename(columns={"position": "number of hets"})
        )

        df6_totalref = (
            df_overlapped.groupby(by="geneID")
            .agg({"refCount": "sum"})
            .reset_index()
            .rename(columns={"totalRef": "total refAllele"})
        )
        df6_totalalt = (
            df_overlapped.groupby(by="geneID")
            .agg({"altCount": "sum"})
            .reset_index()
            .rename(columns={"totalAlt": "total altAllele"})
        )
        df_summary_1 = pd.merge(df6_med, df6_nhets, on=["geneID"], how="inner")
        df_summary_2 = pd.merge(df_summary_1, df6_totalref, on=["geneID"], how="inner")
        df_summary_3 = pd.merge(df_summary_2, df6_totalalt, on=["geneID"], how="inner")
        df_summary_3["totalCount"] = df_summary_3["refCount"] + df_summary_3["altCount"]
        medAltRatio = med_atr
        title = (
            "total alleles >= "
            + str(min_totalcounts)
            + ", single allele >= "
            + str(min_singlecounts)
        )
        med_atr = medAltRatio + np.ndarray.tolist(df6_med["median altRatio"].values)
        # summary_statistics(med_atr,'>>>>>> median alternative allele ratio per gene: '+title)
        nhets = np.ndarray.tolist(df6_nhets["number of hets"].values)
        # summary_statistics(nhets,'>>>>>> number of hets per gene: '+title)
        tCount = tCount + list(df_overlapped["totalCount"].values)
        # summary_statistics(tCount,'>>>>>> read depth per site: '+title)
    else:
        print("file no lines ")
    output1 = out + prefix + "_parsed_mpileup_overlapped_exonicHetSnps.tsv"
    df_overlapped.to_csv(output1, index=False, sep="\t", header=True)
    print(
        "DEBUG USE - Aggregated pasred mpile up data dim %s saved to %s"
        % (df_overlapped.size, output1)
    )
    output2 = out + prefix + "_parsed_mpileup_overlapped_exonicHetSnps_forLambda.tsv"
    df_summary_3.to_csv(output2, index=False, sep="\t", header=True)
    print(
        "DEBUG USE - Number of genes in this file for lambda prediction is %s saved to %s"
        % (df_summary_3.size, output2)
    )
    return (df_overlapped, df_summary_3)


def generate_modelCount(
    prefix, gene_df, out=None, DEBUG=False
):  # force the gene ID to be the target Gene
    out_modelinput = out + prefix + "_parsed_mpileup_exonicHetSnps_model_input.txt"
    unique_gene = gene_df["geneID"].unique()
    rst = ""
    for each in unique_gene:
        idx = each
        gene_lst = [idx, sum(gene_df["geneID"] == each)] + list(
            np.ravel(
                gene_df[gene_df["geneID"] == each][["altCount", "refCount"]].values
            )
        )
        rst += "\t".join([str(x) for x in gene_lst]) + "\n"
    if out is not None:
        file1 = open(out_modelinput, "w")
        file1.write(rst)
        file1.close()
    if DEBUG is True:
        print("DEBUG USE - model input saved to %s" % (out_modelinput))
    # return rst
