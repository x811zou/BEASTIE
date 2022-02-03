#!/usr/bin/env python
# =========================================================================
# This version is written by Scarlett at 01/06/2022
# =========================================================================
import sys
import os
from BEASTIE.extractHets import count_all_het_sites
import logging

# =========================================================================
# main()
# =========================================================================
logging.basicConfig(
    format="%(asctime)-15s [%(levelname)s] %(message)s",
    level=logging.DEBUG,
)
# print("derp")
if len(sys.argv) != 7:
    exit(
        os.path.basename(sys.argv[0])
        + " <sample> <vcfgz-path-filename> <output-path-filename> <chr-start> <chr-start> <gencode_path>\n"
    )
(
    sample,
    vcfFilename,
    output_path_filename,
    chr_start,
    chr_end,
    gencode_path,
) = sys.argv[1:]
# print("derp2")
sample = str(sample)
vcfgz_path_filename = str(vcfFilename)
chr_start = int(chr_start)
chr_end = int(chr_end)
output_path_filename = str(output_path_filename)
gencode_path = str(gencode_path)

count_all_het_sites(
    sample,
    vcfgz_path_filename,
    output_path_filename,
    chr_start,
    chr_end,
    gencode_path,
)

print("done!")

# python run_extractHets.py HG00096 /Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr21/HG00096_chr21.remove_chr.content.SNPs.hets.header.vcf.gz /Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr21/test.txt 21 21 ./BEASTIE/reference/gencode.v19.annotation.filtered.gtf
