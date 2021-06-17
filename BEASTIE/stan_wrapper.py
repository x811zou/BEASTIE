#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
import os
import sys
import math
import pickle
import numpy as np
import statistics
import pandas as pd #(pip install pandas)
import math
import os.path
import time
import calendar
#from Pipe import Pipe
#from StanParser import StanParser
#========================================================================= cyvcf2 package
# The above imports should allow this program to run in both Python 3.8  You might need to update your version of module "future".
# This program can be installed via conda envs located in you lab directory /data/allenlab/scarlett/software/miniconda3/envs.  Anaconda cloud link https://anaconda.org/bioconda/cyvcf2
# -Your current working conda evn is /data/common/shared_conda_envs/miniconda3/bin/conda.
# -To see this output, run 'which conda`
# -Note that the following will change your ~/.bashrc file.  You can backup your .bashrc file prior to running conda init bash or run /data/common/shared_conda_envs/miniconda3/bin/conda init bash to go back to your original conda env.  To install cyvcf2 in your  /data/allenlab/scarlett/software/miniconda3/envs:
# -`cd /data/allenlab/scarlett/software/miniconda3/bin
# -`./conda init bash`
# -Close your working shell.
# -Log back into hardac-login 
# -Confirm your working conda env run `which conda`. S/B /data/allenlab/scarlett/software/miniconda3/bin/conda or ../miniconda3/condabin/conda
# -run `conda create -n cyvcf2`  Creates your cyvcf2 in the .../miniconda3/envs/ directory
# -run `conda activate cyvcf2`    
# -run `conda install -c bioconda cyvcf2`
# Conda installation will start and will display where the environment location for the cyvcf2
# -Once the installation completes `cd /data/allenlab/scarlett/software/miniconda3/envs/cyvcf2/bin`
# -Tests
# -`cyvcf2 --help`  your help options should be displayed.
#=========================================================================
from cyvcf2 import VCF #https://pypi.org/project/cyvcf2/ (pip install cyvcf2)
from parse_mpileup import Parse_mpileup_allChr
from intersect_hets import Intersect_exonicHetSnps,generate_modelCount
from prepare_model import update_model_input_lambda_phasing,significant_genes
#from summary import significant_genes


def main():
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> step1: parse mpileup input")
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   #parsed_mpileup = Parse_mpileup_allChr(prefix,vcfgz_file,pileup_file,min_total_cov,out,DEBUG=DEBUG)
   #input:  HG00096.pileup
   #output: HG00096_parsed_mpileup.tsv
   #local_dir="/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/"
   parsed_mpileup=pd.read_csv(local_dir+"output/"+prefix+"_parsed_mpileup.tsv",sep="\t",header=0,index_col=False)
   print(parsed_mpileup.head())

   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> step2: filter input - only keep exonic het SNPs, file for lambda prediction")
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   start_timestamp = calendar.timegm(time.gmtime())
   parsed_mpileup_exonicHetSNPs,_ = Intersect_exonicHetSnps(parsed_mpileup,meta_file,prefix,out,min_total_cov,min_single_cov,DEBUG=DEBUG)
   stop_timestamp = calendar.timegm(time.gmtime()) 
   print("\tTotal time to complete step2 : %d seconds"%(stop_timestamp-start_timestamp))
   #input:  NA12878_parsed_mpileup.tsv;NA12878_logisticRegression_input.tsv
   #output: NA12878_parsed_mpileup_exonicHetSnps.tsv; NA12878_parsed_mpileup_overlapped_exonicHetSnps_forLambda.tsv
   #contig  position        geneID        rsid      variantID       refAllele       refCount        altAllele       altCount        totalCount      altRatio        if_Indel        if_SV   if_SNP  if_biallelic    lowMAPQDepth    lowBaseQDepth   rawDepth        otherCount 
   # 1       935222  ENSG00000188290.6  rs2298214   rs2298214       C   0       A       3       3       1.0     N       N       Y       Y       0       0       3       0
   #  parsed_mpileup_exonicHetSNPs=pd.read_csv(local_dir+"output/"+str(prefix)+"_parsed_mpileup_exonicHetSnps.tsv",sep="\t",header=0,index_col=False)
   #  print(parsed_mpileup_exonicHetSNPs.head())
   #  init_start_timestamp = calendar.timegm(time.gmtime())
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> step3: save data in format for model input")
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   start_timestamp = calendar.timegm(time.gmtime())
   generate_modelCount(prefix,parsed_mpileup_exonicHetSNPs,out,DEBUG=DEBUG)
   stop_timestamp = calendar.timegm(time.gmtime()) 
   print("\tTotal time to complete step3 : %d seconds"%(stop_timestamp-start_timestamp))
   #  #input:  NA12878_parsed_mpileup_exonicHetSnps.tsv
   #output: NA12878_parsed_mpileup_exonicHetSnps_model_input.txt
   #ENSG00000177757.1	3	1	0	1	0	1	1
   print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> step4: predict phasing error - R")
   ##### 2.use converted format to predict lambda
   start_timestamp = calendar.timegm(time.gmtime())
   cmd="Rscript --vanilla predict_lambda_phasingError.R %s %s %s"%(alpha,local_dir,prefix)
   os.system(cmd)
   stop_timestamp = calendar.timegm(time.gmtime()) 
   print("\tTotal time to complete step4 : %d seconds"%(stop_timestamp-start_timestamp))
   #input: NA12878_parsed_mpileup_overlapped_exonicHetSnps_forLambda.tsv; NA12878_logisticRegression_input.tsv
   #output: NA12878_parsed_mpileup_overlapped_exonicHetSnps_Lambda_0.05_predicted.tsv
   #output: NA12878_logisticRegression_input_phasingErrorpredicted.tsv
   #      geneID      |tota_reads| num_hets | predicted lambda
   #"ENSG00000177757.1"     2	        1	   1.3978282039415
   print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> step5: update model input with phasing")
   start_timestamp = calendar.timegm(time.gmtime())
   update_model_input_lambda_phasing(local_dir,prefix,'pred_error_GIAB')
   stop_timestamp = calendar.timegm(time.gmtime()) 
   print("\tTotal time to complete step 5: %f seconds"%(stop_timestamp-start_timestamp))
   #input:  NA12878_parsed_mpileup_exonicHetSnps_model_input.txt; NA12878_logisticRegression_input_phasingErrorpredicted.tsv
   #output: NA12878_parsed_mpileup_exonicHetSnps_model_input_w_phasingError.txt
   #ENSG00000177757.1	3	1	0	1	0	1	1 prob1 prob2

   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> step6: run BEASTIE and summary model output")
   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   #print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> step4.1: run iBEASTIE")
   ##### 1. use data output from 1 and 2 to run BEASTIE model
   start_timestamp = calendar.timegm(time.gmtime())
   lambdas_file="./output/NA12878_parsed_mpileup_overlapped_exonicHetSnps_Lambda_0.05_predicted.tsv"
   #inFile="./output/NA12878_parsed_mpileup_exonicHetSnps_model_input_w_phasingError.txt"
   inFile="./output/"+prefix+"_parsed_mpileup_exonicHetSnps_model_input.txt"
   in_path=local_dir+"output/"
   out_path=local_dir+"output/output_pkl/"   
   if not os.path.exists( out_path):
      os.makedirs( out_path)
   default_lambda=True
   selfdefined_lambda=False 
   predicted_lambda_modeldata=True
   cmd = "python phasing_stan_wrapper_no_rep_iBEASTIE.py %s %s %s %s %s %s %s %s %s %s %s" % (inFile,sigma,99,model,in_path,out_path,default_lambda,selfdefined_lambda,predicted_lambda_modeldata,lambdas_file,alpha)
   os.system(cmd)
   stop_timestamp = calendar.timegm(time.gmtime()) 
   print("\tTotal time to complete step 6: %d seconds"%(stop_timestamp-start_timestamp))
    #input: NA12878_parsed_mpileup_exonicHetSnps_model_input_w_phasingError.txt
    #output: posterior theta

   #  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> step7: output gene list")
   #  ##### 2. summary
   start_timestamp = calendar.timegm(time.gmtime())
   siGenes=significant_genes(prefix,local_dir,cutoff)
   stop_timestamp = calendar.timegm(time.gmtime()) 
   print("\tTotal time to complete step 7: %d seconds"%(stop_timestamp-start_timestamp))
   #input : output/output_pkl/XX.pickle
   #output: NA12878_ASE.tsv

   init_stop_timestamp = calendar.timegm(time.gmtime()) 
   print("\t>>>>>Overall time needed from step3 to step7 costs %d seconds"%(init_stop_timestamp-init_start_timestamp))

#=========================================================================
# main()
#=========================================================================
if __name__ == '__main__':
    local_dir=sys.argv[1]
    vcfgz_file=sys.argv[2]
    pileup_file=sys.argv[3]
    hetSNP_file=sys.argv[4]
    meta_file=sys.argv[5]
    min_total_cov=sys.argv[6]
    min_single_cov=sys.argv[7]
    prefix=sys.argv[8]
    out=sys.argv[9]
    model=sys.argv[10]
    alpha=sys.argv[11]
    sigma=sys.argv[12]
    cutoff=sys.argv[13]
    DEBUG=False
    #cmd="conda activate cyvcf2"
    #os.system(cmd)
    main()
