#!/usr/bin/env python
# =========================================================================
# This version is written by Scarlett at 02/04/2022
# =========================================================================
import sys
import os
from BEASTIE import extractHets
import logging

# =========================================================================
# main()
# =========================================================================
logging.basicConfig(
    format="%(asctime)-15s [%(levelname)s] %(message)s",
    level=logging.DEBUG,
)

if len(sys.argv) != 8:
    exit(
        os.path.basename(sys.argv[0])
        + " <sample> <vcfgz-path-filename> <output-path-filename> <chr-start> <chr-end> <gencode_path> <debug_gene>\n"
    )
(
    sample,
    vcfFilename,
    output_path_filename,
    chr_start,
    chr_end,
    gencode_path,
    debug_gene,
) = sys.argv[1:]

sample = str(sample)
vcfgz_path_filename = str(vcfFilename)
chr_start = int(chr_start)
chr_end = int(chr_end)
output_path_filename = str(output_path_filename)
gencode_path = str(gencode_path)
debug_gene = str(debug_gene)
include_x_chromosome = False
include_x_chromosome = bool(include_x_chromosome)
if debug_gene == "None":
    debug_gene = None

extractHets.count_all_het_sites(
    vcfgz_path_filename,
    output_path_filename,
    chr_start,
    chr_end,
    gencode_path,
    debug_gene,
)

print("done!")
# [“ENG…”, “ENG…”]
# python run_extractHets.py HG00096 /hpc/group/allenlab/scarlett/output/RNAseq/1000Genome/HG00096/HG00096.no_chr.content.SNPs.hets.vcf.gz /hpc/group/allenlab/scarlett/output/RNAseq/1000Genome/HG00096 1 X /datacommons/allenlab/hg19/filter/gencode_chr/ None
# python run_extractHets.py NA18486 /home/scarlett/data/NA18486/NA18486.no_chr.content.SNPs.hets.filtered.vcf.gz /home/scarlett/data/NA18486/chr12_OAS1.tsv 12 12 /home/scarlett/data/reference/gencode_chr/ ENSG00000089127.8
