#!/usr/bin/env python
# =========================================================================
# This version is written by Scarlett at 01/06/2022
# =========================================================================
import sys
import ProgramName
import logging
from BEASTIE.extractHets import count_all_het_sites
import logging

# =========================================================================
# main()
# =========================================================================
logging.basicConfig(
    format="%(asctime)-15s [%(levelname)s] %(message)s",
    level=logging.DEBUG,
)

if len(sys.argv) != 7:
    exit(
        ProgramName.get()
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
