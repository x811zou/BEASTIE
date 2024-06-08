#!env python


import pandas as pd
import subprocess
import sys

import pandas as pd
import subprocess
import sys
import os
from pathlib import Path
from scipy.stats import hypergeom
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from novelgene_SNPcomposition import SNP_composition

# The probability that we would observe this or an even more imbalanced ratio by chance is about 3.5%. A commonly used significance level is 5%–if we adopt that, we can therefore conclude that our observed imbalance is statistically significant; whales prefer the Atlantic while sharks prefer the Indian ocean."


# main
data_dir = "/data2/1000Genome"
save_dir = "/data2/genes"

# genelist = pd.read_csv("/home/scarlett/data/genes/geneID_list.tsv", sep="\t", header=0)
# gene_list = genelist["geneID"].tolist()  # ["ENSG00000164308.12"]
gene_counter = 0
sig_counter = 0
snp_counter = 0

output = pd.DataFrame()
outputFilename = save_dir + "/significant_fisherexact.tsv"
# finished_genes = save_dir + "/finished_genelist.txt"
# if os.path.isfile(finished_genes):
#    finished_genelist = pd.read_csv(finished_genes, sep="\t", header=0)
#    finished_gene_list = finished_genelist["geneID"].tolist()
#    finished_genes = open(finished_genes, "a")
#    out_stream = open(outputFilename, "a")
# else:
#    finished_genes = open(finished_genes, "w")
#    finished_genes.write("geneID\n")
#    finished_genelist = []
#    out_stream = open(outputFilename, "w")

# python check_SNPs_composition.py ENSG00000177879
gene_candidate = sys.argv[1]

data = SNP_composition(data_dir, save_dir, gene_candidate)
if data.empty:
    print("empty")
sys.exit()
data_snp, data_snp_list = collect_info_SNPs(save_dir, data)
n_total_ind = len(pd.unique(data["sample"]))

print(
    "geneID\tchr\tpos\ttotal_n\tn_have\tn_nothave\tpval\tASE_have\tASE_nothave\tnoASE_have\tnoASE_nothave\n"
)

for pos in data_snp_list:
    if data_snp.to_dict()[pos] >= cutoff_nind:
        ind_w_SNP = data[data["pos"] == pos]["sample"].tolist()
        data_w_SNP = data[data["sample"].isin(ind_w_SNP)]
        data_wo_SNP = data[~data["sample"].isin(ind_w_SNP)]
        n_have_SNP = len(pd.unique(data_w_SNP["sample"]))
        n_nothave_SNP = len(pd.unique(data_wo_SNP["sample"]))
        count = data_snp.to_dict()[pos]
        chrN = data["chr"].iloc[0]
        (
            p,
            ASE_have,
            ASE_nothave,
            noASE_have,
            noASE_nothave,
        ) = calculate_2sidedfisherexact_variant(data_w_SNP, data_wo_SNP)
        info = [
            [
                gene_candidate,
                chrN,
                pos,
                n_total_ind,
                n_have_SNP,
                n_nothave_SNP,
                round(p, 2),
                ASE_have,
                ASE_nothave,
                noASE_have,
                noASE_nothave,
            ]
        ]
        for r in info:
            print("\t".join(map(str, r)))


def process(x):
    cmd = f"cat /hpc/group/allenlab/scarlett/output/RNAseq/1000Genome/{x['sample']}/beastie/beastie_SNPs_even_100/beastie_shapeit2/chr1-22_alignBiasp0.05_ase0.5_s0.5_a0.05_sinCov0_totCov1_W1000K1000/tmp/*_real_alignBiasafterFilter.phasedByshapeit2.tsv | grep \"ENSG00000089127.8\""
    output = subprocess.run(
        cmd, capture_output=True, shell=True, encoding="utf8"
    ).stdout
    num_line = len(output.split("\n"))
    SNPs = []
    print(output)
    for element in output.split("\n"):
        a = element.split("\t")
        SNPs = SNPs + a
    return SNPS.join("_")


def snp_before_genotypingError(x):
    # TP53
    # cmd = f"cat {x['sample']}/beastie/beastie_SNPs_even_100/filterGenotypingError/TEMP.*_hetSNP_intersected_filtered_underGenotypingErTesting.tsv | grep \"ENSG00000089127.8\" | grep \"7572154\""
    # OAS1
    cmd = f"cat {x['sample']}/beastie/beastie_SNPs_even_100/filterGenotypingError/TEMP.*_hetSNP_intersected_filtered_underGenotypingErTesting.tsv | grep \"ENSG00000089127.8\" | grep \"113357193\""
    output = subprocess.run(
        cmd, capture_output=True, shell=True, encoding="utf8"
    ).stdout
    num_line = len(output.split("\t"))
    # print(output)
    # print(num_line)
    if num_line != 1:
        a = output.split("\t")
        # print(a)
        if len(a) > 9:
            # print(a[8])
            return a[8]
        else:
            # print("No error")
            return "No error"
    else:
        # print("Not exit")
        return "Not exist"


def snp_before_alignbiasFilter(x):
    # TP53
    # cmd = f"cat {x['sample']}/beastie/beastie_SNPs_even_100/beastie_shapeit2/chr1-22_alignBiasp0.05_ase0.5_s0.5_a0.05_sinCov0_totCov1_W1000K1000/tmp/*_real_alignBiasbeforeFilter.tsv | grep \"ENSG00000089127.8\" | grep \"7572154\""
    # OAS1
    cmd = f"cat {x['sample']}/beastie/beastie_SNPs_even_100/beastie_shapeit2/chr1-22_alignBiasp0.05_ase0.5_s0.5_a0.05_sinCov0_totCov1_W1000K1000/tmp/*_real_alignBiasbeforeFilter.tsv | grep \"ENSG00000089127.8\" | grep \"113357193\""
    output = subprocess.run(
        cmd, capture_output=True, shell=True, encoding="utf8"
    ).stdout
    num_line = len(output.split("\t"))
    # print(">>")
    if num_line != 1:
        a = output.split("\t")
        # print(a)
        # print(len(a))
        # print(a[15])
        return a[15]
    else:
        # print("Not exit")
        return "Not exist"


############## het1 --> read the file generated from R studio "OAS1_het1.tsv"
# df = pd.read_csv('OAS1_het1.tsv', sep='\t', header=0)
# df['snp'] = df.apply(process, axis=1)
# df.to_csv('OAS1_het1_annotated.tsv', sep='\t')

############## het2
df = pd.read_csv("OAS1_het6.tsv", sep="\t", header=0)
df["113357193_genotypingError"] = df.apply(snp_before_genotypingError, axis=1)
df["113357193_alignmentBias"] = df.apply(snp_before_alignbiasFilter, axis=1)
df["snp"] = df.apply(process, axis=1)
df.to_csv("OAS1_het6_annotated.tsv", sep="\t")

print(df)
