#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
import logging
import os
import sys

import pandas as pd

import BEASTIE.ADM_for_real_data as ADM_for_real_data
import BEASTIE.binomial_for_real_data as binomial_for_real_data
import BEASTIE.run_model_stan_wrapper as run_model_stan_wrapper

from pkg_resources import resource_filename

from .beastie_step1 import create_output_directory
from .prepare_model import (generate_modelCount, significant_genes,
                            update_model_input_lambda_phasing)


def create_file_name(hetSNP_intersect_unique,meta,out,common):
    base = os.path.split(hetSNP_intersect_unique)
    base_modelin = os.path.join(out,'{0}.modelinput.tsv'.format(os.path.splitext(base[1])[0]))
    base_modelin_error = os.path.join(out,'{0}.modelinput_w_error.tsv'.format(os.path.splitext(base[1])[0]))
    base_meta = os.path.split(meta)
    meta_error = os.path.join(out,'{0}.w_prediction.tsv'.format(os.path.splitext(base_meta[1])[0]))
    return base_modelin, base_modelin_error,meta_error

def run(hetSNP_intersect_unique,meta,hetSNP_intersect_unique_forlambda_file,hetSNP_intersect_unique_lambdaPredicted_file,prefix,alpha,model,sigma,in_path,out,cutoff,SAVE_INT,WARMUP,KEEPER,total_cov,either_cov):
    out,common = create_output_directory(in_path,out,sigma,alpha,WARMUP,KEEPER,either_cov,total_cov)
    base_modelin, base_modelin_error,meta_error = create_file_name(hetSNP_intersect_unique,meta,out,common)
    ##########################
    logging.info('=================')
    logging.info('================= Starting step 2.1')
    logging.info('..... start converting data in format for model input')
    logging.info('..... hetSNP_intersect_unique file save to {0}'.format(hetSNP_intersect_unique))
    generate_modelCount(hetSNP_intersect_unique)

    ##########################
    logging.info('=================')
    logging.info('================= Starting step 2.2')
    logging.info('..... start predicting phasing error')
    if (os.path.isfile(meta_error)) and (os.path.isfile(hetSNP_intersect_unique_lambdaPredicted_file)):
        logging.info('..... output file save to {0}'.format(hetSNP_intersect_unique_lambdaPredicted_file))
        logging.info('..... output file save to {0}'.format(meta_error))
    else:
        logging.info('..... data exists, overwrites and saves at {0}'.format(hetSNP_intersect_unique_lambdaPredicted_file))
        logging.info('..... data exists, overwrites and saves at {0}'.format(meta_error))

    predict_lambda_phasing_error = resource_filename('BEASTIE', 'predict_lambda_phasingError.R')
    cmd = f"Rscript --vanilla {predict_lambda_phasing_error} {alpha} {common} {prefix} {model} {hetSNP_intersect_unique} {hetSNP_intersect_unique_forlambda_file} {hetSNP_intersect_unique_lambdaPredicted_file} {meta} {meta_error}"
    os.system(cmd)
    data22=pd.read_csv(hetSNP_intersect_unique_lambdaPredicted_file,sep="\t",header=0,index_col=False)
    logging.info('..... For type 1 error model, input alpha (family wise error rate) is {0}, adjusted after size of input {1} : {2}'.format(alpha,data22.shape[0],alpha/data22.shape[0]))
    #########################
    logging.info('=================')
    logging.info('================= Starting step 2.3')
    logging.info('..... start adding model input with phasing error information')
    if (os.path.isfile(base_modelin_error)):
        logging.info('..... output file save to {0}'.format(base_modelin_error))
    else:
        logging.info('..... data exists, overwrites and saves at : {0}'.format(base_modelin_error))
    update_model_input_lambda_phasing('pred_error_GIAB',base_modelin,base_modelin_error,meta_error)
    ##########################
    logging.info('=================')
    logging.info('================= Starting step 2.4')
    logging.info('..... start running BEASTIE model')
    df_beastie,picklename = run_model_stan_wrapper.run(prefix,base_modelin_error,sigma,alpha,model,out,hetSNP_intersect_unique_lambdaPredicted_file,WARMUP,KEEPER,either_cov,total_cov)
    if df_beastie.shape[0]>2:
        logging.info('..... model output is finished!')
    else:
        logging.error('..... model output is empty, please try again!')
        sys.exit(1)

    df_adm = ADM_for_real_data.run(prefix,base_modelin,out,picklename)
    logging.info('..... finishing running ADM method')

    df_binomial = binomial_for_real_data.run(prefix,base_modelin,out,picklename)
    logging.info('..... finishing running binomial')

    ##########################
    logging.info('=================')
    logging.info('================= Starting step 2.5')
    logging.info('..... start generating gene list')
    outfilename=out+"/"+prefix+"_ASE_all.tsv"
    outfilename_ase=out+"/"+prefix+"_ASE_cutoff_"+str(cutoff)+"_filtered.tsv"
    if (os.path.isfile(outfilename)) and (os.path.isfile(outfilename_ase)):
        logging.info('..... data exists, overwrites and saves at {0}'.format(outfilename))
        logging.info('..... data exists, overwrites and saves at {0}'.format(outfilename_ase))
    else:
        logging.info('..... output file save to {0}'.format(outfilename))
        logging.info('..... output file save to {0}'.format(outfilename_ase))
    significant_genes(df_beastie,df_binomial,df_adm,outfilename,outfilename_ase,cutoff,hetSNP_intersect_unique_lambdaPredicted_file)
    if str(SAVE_INT) == "True":
        for files in os.listdir(out):
            if "TEMP" in files:
                logging.info('..... remove TEMP files {0}'.format(files))
                os.remove(out+"/"+files)
    logging.info('=================')
    logging.info('>>  Finishing running BEASTIE!')
    data25_1=pd.read_csv(outfilename,sep="\t",header=0,index_col=False)
    data25_2=pd.read_csv(outfilename_ase,sep="\t",header=0,index_col=False)
    if data25_1.shape[0]>2 or data25_2.shape[0]>2:
        logging.info('..... Yep! You are done!')
    else:
        logging.error('..... output is empty, please try again!')
        sys.exit(1)
