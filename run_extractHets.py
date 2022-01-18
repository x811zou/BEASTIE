#!/usr/bin/env python
# =========================================================================
# This version is written by Scarlett at 01/06/2022
# =========================================================================
import sys
from BEASTIE.misc_tools import ProgramName
from BEASTIE.extractHets import count_all_het_sites

# import count_all_het_sites

# =========================================================================
# main()
# =========================================================================
# if len(sys.argv) != 6:
#     exit(
#         ProgramName.get()
#         + " <sample> <vcfFilename> <full-output-path-filename> <chr-start> <chr-start> <gencode_path>\n"
#     )
# (sample, vcfFilename, outputFilename, chr_start, chr_end, gencode_path) = sys.argv[1:]
# chr_start = int(chr_start)
# chr_end = int(chr_end)
# outputFilename=full_path_output_filename
# chr_start=1
# chr_end=22
# gencode_path="/datacommons/allenlab/scarlett"

# count_all_het_sites(tmp=None, sample, vcfFilename, outputFilename, chr_start, chr_end,gencode_path)

count_all_het_sites(
    "HG00096",
    "/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr21/HG00096_chr21.remove_chr.content.SNPs.hets.header.vcf.gz",
    "/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr21/test.txt",
    21,
    21,
    "./BEASTIE/reference/gencode_chr",
)
