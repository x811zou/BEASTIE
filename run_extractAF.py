#!/usr/bin/env python
# =========================================================================
# This version is written by Scarlett at 02/04/2022
# =========================================================================
import sys
import os
from BEASTIE import extractAF
import logging

# =========================================================================
# main()
# =========================================================================
logging.basicConfig(
    format="%(asctime)-15s [%(levelname)s] %(message)s",
    level=logging.DEBUG,
)
# print("derp")
if len(sys.argv) != 6:
    exit(
        os.path.basename(sys.argv[0])
        + " <vcfgz-path> <output-path-filename> <chr-start> <chr-start> <gencode_path> \n"
    )
(
    vcfFile_path,
    output_path_filename,
    chr_start,
    chr_end,
    gencode_path
) = sys.argv[1:]
# print("derp2")
vcfgz_path = str(vcfFile_path)
chr_start = str(chr_start)
chr_end = str(chr_end)
output_path_filename = str(output_path_filename)
gencode_path = str(gencode_path)

extractAF.count_all_sites(
    vcfgz_path,
    output_path_filename,
    chr_start,
    chr_end,
    gencode_path
)

print("done!")
# [“ENG…”, “ENG…”]
# python /hpc/group/allenlab/scarlett/script/BEASTIE/run_extractAF.py /datacommons/allenlab/scarlett/data/VCF/1000_genome/20130502/bgzip /hpc/group/allenlab/scarlett/output/RNAseq/1000Genome/AF.csv X X /datacommons/allenlab/hg19/filter/gencode_chr/
