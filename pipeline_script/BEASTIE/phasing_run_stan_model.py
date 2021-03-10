#!/usr/bin/env python
import pickle
import os
from sys import argv
#from sklearn import preprocessing

model="BEASTIE3"
folder="NA12878"
in_path="/data/allenlab/scarlett/data/ASE_simulation/BEASTIE_empirical_data/"
sigma=0.5
##################################################################### do not change
s=sigma
out_path="/data/allenlab/scarlett/software/cmdstan/examples/"+model+"/"+folder+"/output_pkl/"
if not os.path.exists(out_path):
    os.makedirs(out_path)

all_data = os.listdir(in_path)
all_data.sort()
#all_data = all_data[int(argv[1]):int(argv[1])+500]

all_data = os.listdir(in_path)[int(argv[1]):int(argv[1])+10000]
for files in all_data:
    if ("g-1000_" in files and "s14" in files) or ("g-1000_" in files and "s15" in files) or ("g-1000_" in files and "s16" in files) or ("g-1000_" in files and "s17" in files) or ("g-1000_" in files and "s12" in files) or ("g-1000_" in files and "s13" in files):
        #for s in sigma:
    	cmd = "python /data/allenlab/scarlett/python/stan_model/phasing_stan_wrapper_no_rep.py %s %s %s %s %s %s" % (files,s,argv[1],model,in_path,folder)
    	os.system(cmd)
#error=("03","05","10","20","30","40","50")
#in_path="/data/reddylab/scarlett/1000G/data/ASE_simulation/simulated_data_%s_error/"%(str(n))

print("Done")


