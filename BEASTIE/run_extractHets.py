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
# print("derp")
if len(sys.argv) != 8:
    exit(
        os.path.basename(sys.argv[0])
        + " <sample> <vcfgz-path-filename> <output-path-filename> <chr-start> <chr-start> <gencode_path> <debug_gene>\n"
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
# print("derp2")
sample = str(sample)
vcfgz_path_filename = str(vcfFilename)
chr_start = int(chr_start)
chr_end = int(chr_end)
output_path_filename = str(output_path_filename)
gencode_path = str(gencode_path)
debug_gene = str(debug_gene)
if debug_gene == "None":
    debug_gene = None

extractHets.count_all_het_sites(
    sample,
    vcfgz_path_filename,
    output_path_filename,
    chr_start,
    chr_end,
    gencode_path,
    debug_gene,
)

print("done!")

# python run_extractHets.py HG00099 /Users/scarlett/allenlab/BEASTIE_other_example/HG00099_50M/HG00099.remove_chr.content.SNPs.hets.header.vcf.gz /Users/scarlett/allenlab/BEASTIE_other_example/HG00099_50M/HG00099_hetsnp 10 10 /Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE/reference/gencode_chr ENSG00000237399.3
