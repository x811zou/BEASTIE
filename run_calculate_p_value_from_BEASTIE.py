#!/usr/bin/env python
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
import sys
# example commands: python run_calculate_p_value_from_BEASTIE.py /data2/simulation/semi_empirical/CEU/g-1000 /data2/stan/iBEASTIE4/sigma0.7/semi_empirical/CEU/g-1000/output_pkl /home/scarlett/github/RNAseq-analysis/stan_models/iBEASTIE4 0.7
# example commands: python run_calculate_p_value_from_BEASTIE.py /data2/simulation/parametrized/ASE_0.05_error /data2/stan/BEASTIE3-pi0.05/sigma0.7/parametrized/ASE_0.05_error/output_pkl /home/scarlett/github/RNAseq-analysis/stan_models/BEASTIE3-pi0.05 0.7
# write a script to excute the calculate_p_value_from_qb.py on files in a directory and save the results in another directory with the same name
# pass an input directory and output directory
in_path = argv[1]
out_path = argv[2]
model = argv[3]
sigma = argv[4]

all_data = os.listdir(in_path)
all_data.sort()
# check whether output directory exists, if not, create one
if not os.path.exists(out_path):
    os.makedirs(out_path)

for files in all_data:
    outname = f"_s-{sigma}.tsv"
    outname = os.path.basename(files).replace(".txt", "")  + outname
    tsv_file = os.path.join(out_path, outname)
    if (("g-1000_h-1" in files in files) or ("g-1000_h-3" in files in files)) and (("t-0.5" in files) or ("t-1" in files)):
        print(">> start with %s !"%(files))
        cmd = (
            "python ./calculate_p_value_from_BEASTIE.py %s %s %s %s"
            % (f"{in_path}/{files}", f"{tsv_file}",model,sigma)
        )
        os.system(cmd)
        print(">> saved %s !"%(files))