import pandas as pd
import os
from pathlib import Path
from statsmodels.stats.multitest import multipletests
import sys
import logging
import numpy as np
from scipy.stats import fisher_exact
"""
python novelgene_SNPcomposition.py ENSG00000164308

geneID  pvals_corrected combined_SNPs   sample  ASE_gene
ENSG00000164308 1.0469791929760155e-76  96211741_96231000_96232142_96232286_96232402_96235896_96237326_96245343_96245439_96245518_96245617_96245892_96249115_96254209_96254354_96254817   NA12829 ASE
"""

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_and_process_file(cutoff_totalcount: int, file_path: str, geneID: str, is_counts: bool) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep="\t", header=0)
    if is_counts:
        df['Gene.stable.ID'] = df['geneID'].str.split('.').str[0]
        df["total_count"] = df["totalCount"]
        df = df[df["Gene.stable.ID"] == geneID]
        df = df[df["total_count"] >= cutoff_totalcount]
    else:
        # Perform FDR correction
        _, pvals_corrected, _, _ = multipletests(df['mode_st_p_value'], method='fdr_bh')
        df['Gene.stable.ID'] = df['geneID']
        df['pvals_corrected'] = pvals_corrected
        df = df[["Gene.stable.ID", "n_hets","pvals_corrected"]][df["Gene.stable.ID"] == geneID]
    return df

def SNP_composition(cutoff_totalcount: int, data_dir: str, geneID: str, gene_name: str, SNP_composition_file: str) -> pd.DataFrame:
    if os.path.isfile(SNP_composition_file):
        df_combined = pd.read_csv(
            SNP_composition_file, sep="\t", header=0
        )
    else:
        logging.info("No pre-existed combined SNP composition information.")
        subfolders = [f.path for f in os.scandir(data_dir) if f.is_dir()]
        df_combined = pd.DataFrame()
        for sub_dir in subfolders:
            sample = os.path.basename(sub_dir)
            counts_file = os.path.join(sub_dir, "beastie/runModel_phased_even100/chr1-22_alignBiasp0.05_s0.7_a0.05_sinCov0_totCov1_W1000K1000/tmp", f"{sample}_real_alignBiasafterFilter.phasedByshapeit2.cleaned.tsv")
            qb_file = os.path.join(sub_dir, "qb", f"{sample}_qb_highestsite.tsv")

            if os.path.isfile(counts_file) and os.path.isfile(qb_file):
                df_counts = read_and_process_file(cutoff_totalcount, counts_file, geneID, True)
                df_qb = read_and_process_file(cutoff_totalcount, qb_file, geneID, False)

                if not df_counts.empty and not df_qb.empty:
                    chr = df_counts["chr"].iloc[0]
                    SNP_list = df_counts["pos"].tolist()
                    combined_SNPs = "_".join(map(str, SNP_list))
                else:
                    continue
                df_qb["Gene.name"] = gene_name
                df_qb["chr"] = chr
                df_qb["combined_SNPs"] = combined_SNPs
                df_qb["sample"] = sample
                df_qb["ASE_gene"] = "ASE" if df_qb["pvals_corrected"].iloc[0] < 0.05 else "no_ASE"
                df_combined = pd.concat([df_combined, df_qb]) if not df_combined.empty else df_qb

        df_combined.to_csv(SNP_composition_file, sep="\t", index=False)
    if df_combined.empty:
        logging.info("No data found for the specified gene.")
    else:
        logging.info(f"combined SNP composition processed successfully for gene {gene_candidate}.")
    return df_combined

def collect_info_one_gene(cutoff_totalcount: int, data_dir: str, geneID: str, gene_name: str, outfile: str):
    if os.path.isfile(common_SNP_file):
        df_combined = pd.read_csv(
            outfile, sep="\t", header=0
        )
    else:
        logging.info("No pre-existed common SNP information.")
        subfolders = [f.path for f in os.scandir(data_dir) if f.is_dir()]
        df_combined = pd.DataFrame()
        for sub_dir in subfolders:
            sample = os.path.basename(sub_dir)
            counts_file = os.path.join(sub_dir, "beastie/runModel_phased_even100/chr1-22_alignBiasp0.05_s0.7_a0.05_sinCov0_totCov1_W1000K1000/tmp", f"{sample}_real_alignBiasafterFilter.phasedByshapeit2.cleaned.tsv")
            qb_file = os.path.join(sub_dir, "qb", f"{sample}_qb_highestsite.tsv")

            if os.path.isfile(counts_file) and os.path.isfile(qb_file):
                df_counts = read_and_process_file(cutoff_totalcount, counts_file, geneID, True)
                df_qb = read_and_process_file(cutoff_totalcount, qb_file, geneID, False)

                if not df_counts.empty and not df_qb.empty:
                    df_counts["Gene.name"] = gene_name
                    df_counts = df_counts[["geneID","Gene.stable.ID", "Gene.name","chr","pos","rsid","ref","refCount","alt","altCount","total_count"]]
                else:
                    continue
                df_counts["sample"] = sample
                df_counts["ASE_gene"] = "ASE" if df_qb["pvals_corrected"].iloc[0] < 0.05 else "no_ASE"
                df_counts["compare_2allele"] = np.where(
                    df_counts["refCount"] > df_counts["altCount"],
                    "REF",
                    "ALT",
                )
                df_combined = pd.concat([df_combined, df_counts]) if not df_combined.empty else df_counts
        df_combined.to_csv(outfile, sep="\t", index=False)
    return df_combined

def collect_info_SNPs_combination(data, outfile):
    data_snp_list = data["combined_SNPs"].value_counts().index.tolist()
    data_snp = data["combined_SNPs"].value_counts(sort=True, ascending=False)

    #geneID = data["Gene.stable.ID"].iloc[0]
    data_snp.to_csv(
        outfile, sep="\t"
    )
    return data_snp, data_snp_list

def collect_info_common_SNPs(data, outfile):
    data_snp_list = data["pos"].value_counts().index.tolist()
    data_snp = data["pos"].value_counts(sort=True, ascending=False)
    #geneID = data["Gene.stable.ID"].iloc[0]
    data_snp.to_csv(outfile, sep="\t")
    return data_snp, data_snp_list

def get_variant_combination(data, data_snp, data_snp_list, cutoff_nind, n_total_ind):
    print(
        "geneID\tgeneName\tchr\tSNPs_combination\ttotal_n\tn_SNP\tn_noSNP\tn_individual\tpercentage_individual\tpval\tASE_SNP\tnoASE_SNP\tASE_noSNP\tnoASE_noSNP\n"
    )
    for pos in data_snp_list:
        if data_snp.to_dict()[pos] >= cutoff_nind:
            n_individual = data_snp.to_dict()[pos]
            data_w_SNP = data[data["combined_SNPs"] == pos]
            data_wo_SNP = data[data["combined_SNPs"] != pos]
            n_have_SNP = len(pd.unique(data_w_SNP["sample"]))
            n_nothave_SNP = len(pd.unique(data_wo_SNP["sample"]))
            chrN = data["chr"].iloc[0]
            (
                p,
                ASE_have,
                noASE_have,
                ASE_nothave,
                noASE_nothave,
            ) = calculate_2sidedfisherexact_variant(data_w_SNP, data_wo_SNP)
            gene_candidate = data["Gene.stable.ID"].iloc[0]
            gene_name = data["Gene.name"].iloc[0]
            n_hets = data["n_hets"].iloc[0]
            if p <= 0.1:
                info = [
                    [
                        gene_candidate,
                        gene_name,
                        n_hets,
                        chrN,
                        pos,
                        n_total_ind,
                        n_have_SNP,
                        n_nothave_SNP,
                        n_individual,
                        round(n_individual/n_total_ind,2),
                        round(p, 5),
                        ASE_have,
                        noASE_have,
                        ASE_nothave,
                        noASE_nothave,
                    ]
                ]
                for r in info:
                    print("\t".join(map(str, r)))

def get_variant(data, data_snp, data_snp_list, cutoff_nind, n_total):
    print(
        "geneID\tgeneName\tchr\tSNP\tvariant\ttotal_n\tn_SNP\tn_noSNP\tn_individual\tpercent_individual\tpval\tASE_SNP\tnoASE_SNP\tASE_noSNP\tnoASE_noSNP\n"
    )

    for pos in data_snp_list:
        if data_snp.to_dict()[pos] >= cutoff_nind:
            ind_w_SNP = data[data["pos"] == pos]["sample"].tolist()
            data_w_SNP = data[data["sample"].isin(ind_w_SNP)]
            data_wo_SNP = data[~data["sample"].isin(ind_w_SNP)]
            n_have_SNP = len(pd.unique(data_w_SNP["sample"]))
            n_nothave_SNP = len(pd.unique(data_wo_SNP["sample"]))
            count = data_snp.to_dict()[pos]
            rsid = data_w_SNP["rsid"].iloc[0]
            chrN = data_w_SNP["chr"].iloc[0]
            (
                p,
                ASE_have,
                ASE_nothave,
                noASE_have,
                noASE_nothave,
            ) = calculate_2sidedfisherexact_variant(data_w_SNP, data_wo_SNP)
            gene_candidate = data["Gene.stable.ID"].iloc[0]
            gene_name = data["Gene.name"].iloc[0]
            if p <= 0.1:
                info = [
                    [
                        gene_candidate,
                        gene_name,
                        chrN,
                        pos,
                        rsid,
                        n_total,
                        n_have_SNP,
                        n_nothave_SNP,
                        count,
                        round(count/n_total,2),
                        round(p, 5),
                        ASE_have,
                        ASE_nothave,
                        noASE_have,
                        noASE_nothave,
                    ]
                ]
                for r in info:
                    print("\t".join(map(str, r)))

def calculate_2sidedfisherexact_variant(data_w_SNP, data_wo_SNP):
    data_w_SNP = data_w_SNP[["sample", "ASE_gene"]]
    data_w_SNP = data_w_SNP.drop_duplicates()
    data_wo_SNP = data_wo_SNP[["sample", "ASE_gene"]]
    data_wo_SNP = data_wo_SNP.drop_duplicates()
    ASE_SNP = (data_w_SNP["ASE_gene"] == "ASE").sum()
    ASE_noSNP = (data_wo_SNP["ASE_gene"] == "ASE").sum()
    noASE_SNP = (data_w_SNP["ASE_gene"] == "no_ASE").sum()
    noASE_noSNP = (data_wo_SNP["ASE_gene"] == "no_ASE").sum()
    #                ASE   no ASE
    # have SNP       188     0
    # not have SNP    7     70
    table = np.array([[ASE_SNP, noASE_SNP], [ASE_noSNP, noASE_noSNP]])
    oddsr, p = fisher_exact(table, alternative="two-sided")
    return (
        p,
        ASE_SNP,
        noASE_SNP,
        ASE_noSNP,
        noASE_noSNP,
    )

if __name__ == "__main__":
    if len(sys.argv) < 3:
        logging.error("Gene candidate argument missing.")
        sys.exit(1)

    gene_candidate = sys.argv[1]
    gene_name = sys.argv[2]

    cutoff_totalcount = 10
    cutoff_nind = 1
    data_dir = "/data2/1000Genome"
    save_dir = "/data2/genes"

    output_dir = Path(os.path.join(save_dir, gene_name))
    output_dir.mkdir(parents=True, exist_ok=True)

    #######################################
    ####################################### step1: Combined SNPs
    #######################################
    SNP_composition_file = output_dir / "combined_SNP_composition.tsv"
    combined_SNP_frequency_file = output_dir / "combined_SNP_frequency_table.tsv"
        ## (1) SNP composition
    combined_SNP = SNP_composition(cutoff_totalcount, data_dir, gene_candidate, gene_name, SNP_composition_file)
    """
Gene.stable.ID  n_hets  pvals_corrected chr     combined_SNPs   sample  ASE_gene
ENSG00000164308 16      1.0469791929760155e-76  chr5    96231000_96232142_96237326_96245343_96245439_96249115_96254209_96254354_96254817        NA12829 ASE
ENSG00000164308 13      1.8654808436751065e-53  chr5    96231000_96232142_96237326_96245343_96245439_96249115_96254209_96254354_96254817        NA12272 ASE
    """
        ## (2) collect all SNP combination frequency information
    combined_SNP_data_snp, combined_SNP_data_snp_list = collect_info_SNPs_combination(combined_SNP, combined_SNP_frequency_file)
    n_total_ind_option1 = len(pd.unique(combined_SNP["sample"]))
    """"
combined_SNPs   count
96231000_96232142_96237326_96245343_96245439_96249115_96254209_96254354_96254817        99
96253448_96254301_96255017      7
    """

    #######################################
    ####################################### step2: common SNP identification
    #######################################
    common_SNP_file = output_dir / "common_SNP_composition.tsv"
    common_SNP_frequency_file = output_dir / "common_SNP_frequency_table.tsv"

        ## (1) common SNP
    """
geneID  Gene.stable.ID  chr     pos     rsid    ref     refCount        alt     altCount        total_count      sample  ASE_gene        compare_2allele
ENSG00000164308.12      ENSG00000164308 chr5    96211741        rs1230358       T       1       G       0       1       NA12829 ASE     REF
ENSG00000164308.12      ENSG00000164308 chr5    96231000        rs2549782       G       88      T       1       89      NA12829 ASE     REF
    """
    common_SNP = collect_info_one_gene(cutoff_totalcount, data_dir, gene_candidate, common_SNP_file)
        ## (2) collect all SNP combination frequency information
    """
    pos     count
96254817        219
96249115        216
    """
    common_SNP_data_snp, common_SNP_data_snp_list = collect_info_common_SNPs(
        common_SNP , common_SNP_frequency_file
    )
    n_total_ind_option2 = len(pd.unique(common_SNP["sample"]))

    if n_total_ind_option1 != n_total_ind_option2:
        logging.raiseExceptions
    n_total = n_total_ind_option1

    print(">>>>>>>>>>")
    print(
        f">>>>>>>>>> we only keep SNP(s) with >= {cutoff_totalcount} total counts, and shared in at least {cutoff_nind} individuals"
    )
    print(">>>>>>>>>> variant combination")
    #######################################
    ####################################### option1: variant_combination
    #######################################
    get_variant_combination(
        combined_SNP,
        combined_SNP_data_snp,
        combined_SNP_data_snp_list,
        cutoff_nind,
        n_total,
    )
    print(">>>>>>>>>> single variant")

    #######################################
    ####################################### option2: variant single
    #######################################
    get_variant(
        common_SNP,
        common_SNP_data_snp,
        common_SNP_data_snp_list,
        cutoff_nind,
        n_total,
    )
    #######################################
    ####################################### option3: haplotype
    #######################################