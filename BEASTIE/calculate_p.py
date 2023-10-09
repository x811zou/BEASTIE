#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

import csv
import os
import subprocess
import numpy as np
from scipy.stats import t, norm
import matplotlib.pyplot as plt
import pandas as pd
import sys
from pathlib import Path
from scipy import stats
from scipy.stats import t,lognorm, gamma, weibull_min, expon, beta
import math
import pickle
import time
sys.path.append('/home/scarlett/github/BEASTIE')
from BEASTIE.run_model_stan_wrapper import runModel, summarize,computeBeastieScoreLog2
from BEASTIE.binomial_for_real_data import getBaseline,getBaseline_pooled
from BEASTIE.run_model_stan_wrapper import save_raw_theta_parallel
################################ command    
# python calculate_p.py /home/scarlett/github/ipython_notebook/BEASTIE_ipn/13_distribution_p/out /home/scarlett/github/BEASTIE/BEASTIE/iBEASTIE4 1 30 100 100
################################

# self-defined functions
def fit_data_to_dist(data, distribution):
    params = distribution.fit(data)
    return params

def fit_distribution(x,distribution):
    params = fit_data_to_dist(x, distribution)
    fitted_distribution = distribution.pdf(x, *params)
    return params,fitted_distribution

def generate_fields(geneID,M, D, theta):
    # calculate probability for binomial distribution
    p = theta / (1.0 + theta)
    # calculate alternative and reference read counts for each het
    alt_counts = np.random.binomial(D, p, M)
    ref_counts = D - alt_counts

    # construct output fields
    fields = [geneID, str(M)]
    for alt, ref in zip(alt_counts, ref_counts):
        fields.append(str(alt))
        fields.append(str(ref))

    # Adding number of missing pi and phasing error fields
    fields.append("0")  # Number of missing pi
    fields += ['0.000001'] * (M - 1)  # Phasing error (We assume all pairs to be -1)
    line="\t".join(fields)
    fields = line.rstrip().split()
    return fields

def get_thetas(out,fields,model,sigma,WARMUP=1000,KEEPER=1000):
    tmp_output_file=out+"/tmp.0.file"
    stan_output_file=out+"/stan.0.out"
    init_file=out+"/init.0.file"
    geneID, thetas = runModel(
        model,
        fields,
        tmp_output_file,
        stan_output_file,
        init_file,
        sigma,
        WARMUP,
        KEEPER,
        phasing_method="shapeit2",
    )
    return geneID, thetas

def calculate_right_tail_probability(values, X):
    # Sort the list in ascending order
    sorted_values = sorted(values)
    
    # Find the rank of X in the sorted list
    rank_x = sum(1 for value in sorted_values if value <= X)
    
    # Calculate the right tail probability
    tail_probability = (len(values) - rank_x + 1) / (len(values) + 1)
    
    return tail_probability

def get_stats(thetas,lambdas_choice_gam):
    mean,median,variance,CI_left,CI_right,mad,log2_mean,log2_median,log2_variance,abslog2_mean,abslog2_median,abslog2_variance = summarize(thetas,alpha=0.05)
    log2_thetas = np.log2(np.array(thetas))
    _, sum_prob_lambda_gam = computeBeastieScoreLog2(log2_thetas, lambdas_choice_gam)
    return mad,log2_mean,log2_median,abslog2_mean,abslog2_median,log2_variance,abslog2_variance,sum_prob_lambda_gam

def give_pkl_name(out_path,numgenes,numhets,readperhet,theta,sigma):
    names = f"g-{numgenes}_h-{numhets}_d-{readperhet}_t-{theta}_s-{sigma}.pickle"
    thetas_file = os.path.join(out_path, names)
    return thetas_file

def get_ALT(ALT_sample,model,sigma,lambdas_choice_gam,name):
    # ALT
    thetas = get_thetas(ALT_sample,model,sigma,WARMUP=1000,KEEPER=1000)
    mad, log2_mean, log2_median, abslog2_mean, abslog2_median, log2_variance, abslog2_variance, sum_prob_lambda_gam = get_stats(thetas,lambdas_choice_gam=lambdas_choice_gam)
    if name == "mad":
        ALT = mad
    elif name == "beastie_score":
        ALT = sum_prob_lambda_gam
    elif name == "mean-log2(theta)":
        ALT = log2_mean
    elif name =="median-log2(theta)":
        ALT = log2_median
    elif name =="median-abs(log2(theta))":
        ALT = abslog2_median
    elif name == "mean_over_std":
        ALT = log2_mean/math.sqrt(log2_variance)
    elif name == "abs_mean_over_std":
        ALT = abslog2_mean/math.sqrt(abslog2_variance)
    return ALT

def sim_nulls(out_path,numNulls,model,numhets,readperhet,theta,sigma,value_name="mean_over_std",lambdas_choice_gam=1.2):
    thetas_dict = check_nulls(out_path,numNulls,model,numhets,readperhet,theta=theta,sigma=sigma)
    # calculating stats
    list_mad = []
    list_posterior_mass_support_alt=[]
    list_log2_mean=[]
    list_log2_median=[]
    list_abslog2_median=[]
    list_mean_over_std=[]
    list_abs_mean_over_std=[]
    list_posterior_mass_support_alt=[]
    for geneID in thetas_dict:
        gene_thetas = thetas_dict[geneID][0]
        #NS_p = thetas_dict[geneID][1]
        #MS_p = thetas_dict[geneID][2]
        mad,log2_mean,log2_median,abslog2_mean,abslog2_median,log2_variance,abslog2_variance,sum_prob_lambda_gam = get_stats(gene_thetas,lambdas_choice_gam=lambdas_choice_gam)
        list_mad.append(mad)
        list_log2_mean.append(log2_mean)
        list_log2_median.append(log2_median)
        list_abslog2_median.append(abslog2_median)
        list_mean_over_std.append(log2_mean/math.sqrt(log2_variance))
        list_abs_mean_over_std.append(abslog2_mean/math.sqrt(abslog2_variance))
        list_posterior_mass_support_alt.append(sum_prob_lambda_gam)
    if value_name == "mad":
        LIST = list_mad
    elif value_name == "beastie_score":
        LIST = list_posterior_mass_support_alt
    elif value_name == "log2_mean":
        LIST = list_log2_mean
    elif value_name =="log2_median":
        LIST = list_log2_median
    elif value_name =="abslog2_median":
        LIST = list_abslog2_median
    elif value_name == "mean_over_std":
        LIST = list_mean_over_std
    elif value_name == "abs_mean_over_std":
        LIST = list_abs_mean_over_std
    return LIST


def check_nulls(out_path,numNulls,model,numhets,readperhet,theta=1,sigma=0.7):
    thetas_file = give_pkl_name(out_path,numNulls,numhets,readperhet,theta,sigma)
    if os.path.exists(thetas_file):
        print(f"... file {thetas_file} already exists.")
        thetas_dict = pickle.load(open(thetas_file, "rb"))
    else:
        print(f"... file {thetas_file} start generating. This may take a while.")
        thetas_dict = {}
        start_time = time.time()  # Get the start time
        for i in range(0, numNulls):
            # Check progress and print
            if i % 1000 == 0:  # Check every 1000 iterations
                elapsed_time = time.time() - start_time
                print(f"... {(i/numNulls)*100:.2f}% done. Time elapsed: {elapsed_time:.2f} seconds.")
            geneID='gene'+str(i+1)+'_t'+str(theta)
            field = generate_fields(geneID,M=numhets, D=readperhet, theta=theta)
            geneID,thetas = get_thetas(out_path,field,model,sigma,WARMUP=1000,KEEPER=1000)
            _,_,_,MS_p = getBaseline(field)
            _,NS_p,_,_=getBaseline_pooled(field)
            thetas_dict[geneID]=(thetas,NS_p,MS_p)
        # Save the thetas_dict
        pickle.dump(thetas_dict, open(thetas_file, "wb"))
        # Print the final time when done
        total_time = time.time() - start_time
        print(f"100% done. Total time: {total_time:.2f} seconds.")
    return thetas_dict

def sim_alts(out_path,numsamples,model,numhets,readperhet,theta,sigma,value_name,lambdas_choice_gam=1.2,NULL_list=None):
    thetas_dict = check_nulls(out_path,numsamples,model,numhets,readperhet,theta=theta,sigma=sigma)
    # calculating stats
    list_mad = []
    list_posterior_mass_support_alt=[]
    list_log2_mean=[]
    list_log2_median=[]
    list_abslog2_median=[]
    list_mean_over_std=[]
    list_abs_mean_over_std=[]
    em_p_list = []
    NS_p_list = []
    MS_p_list = []

    for geneID in thetas_dict:
        gene_thetas = thetas_dict[geneID][0]
        NS_p = thetas_dict[geneID][1]
        MS_p = thetas_dict[geneID][2]
        mad,log2_mean,log2_median,abslog2_mean,abslog2_median,log2_variance,abslog2_variance,sum_prob_lambda_gam = get_stats(gene_thetas,lambdas_choice_gam=lambdas_choice_gam)
        list_mad.append(mad)
        list_log2_mean.append(log2_mean)
        list_log2_median.append(log2_median)
        list_abslog2_median.append(abslog2_median)
        mean_over_std=log2_mean/math.sqrt(log2_variance)
        list_mean_over_std.append(mean_over_std)
        abs_mean_over_std=abslog2_mean/math.sqrt(abslog2_variance)
        list_abs_mean_over_std.append(abs_mean_over_std)
        NS_p_list.append(NS_p)
        MS_p_list.append(MS_p)
        if NULL_list is not None:
            em_p = calculate_right_tail_probability(NULL_list, mean_over_std)
            #print(f"find the tail prob of {abs_mean_over_std} among {NULL_list}")
            em_p_list.append(em_p)
        # list_posterior_mass_support_alt.append(sum_prob_lambda_gam)
    if value_name == "mad":
        LIST = list_mad
    elif value_name == "beastie_score":
        LIST = list_posterior_mass_support_alt
    elif value_name == "log2_mean":
        LIST = list_log2_mean
    elif value_name =="log2_median":
        LIST = list_log2_median
    elif value_name =="abslog2_median":
        LIST = list_abslog2_median
    elif value_name == "mean_over_std":
        LIST = list_mean_over_std
    elif value_name == "abs_mean_over_std":
        LIST = list_abs_mean_over_std 
    return LIST,NS_p_list,MS_p_list,em_p_list

def read_csv(outlist):
    if os.path.exists(outlist):
        # Open and read the CSV file
        dist_p_values = []
        em_p_values = []
        NS_p_list = []
        MS_p_list = []
        with open(outlist, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader)
            # Read and store each column's values
            for row in csvreader:
                dist_p_values.append(float(row[0]))
                em_p_values.append(float(row[1]))
                NS_p_list.append(float(row[2]))
                MS_p_list.append(float(row[3]))
        return dist_p_values,em_p_values,NS_p_list,MS_p_list
    else:
        print(f"... file {outlist} does not exist.")

def simulate_null_ans_calP(out_path,numNulls,num_samples,model,sigma,lambdas_choice_gam,numhets,readperhet,value_name="mean_over_std",save=False):
    sub_out_path = out_path + "/"+str(value_name)
    # Check if output directory exists
    if not os.path.exists(sub_out_path):
        os.makedirs(sub_out_path)
        print(f"... directory {sub_out_path} created.")
    else:
        print(f"... directory {sub_out_path} already exists.")
    # Check if output directory exists
    os.makedirs(out_path, exist_ok=True)
    NULL_list=sim_nulls(out_path,numNulls,model,numhets,readperhet,theta=1,sigma=sigma,value_name=value_name,lambdas_choice_gam=lambdas_choice_gam)
    distributions = [t,norm,lognorm, gamma, weibull_min, expon, beta]
    distribution_names = ["t","normal","Log-normal", "Gamma", "Weibull", "Exponential", "Beta"]

    #values = [list_mad,list_posterior_mass_support_alt,list_log2_mean,list_log2_median,list_abslog2_median,list_mean_over_var,list_abs_mean_over_var]
    #value_names = ["mad","beastie_score", "mean-log2(theta)", "median-log2(theta)", "median-abs(log2(theta))","mean_over_var","abs_mean_over_var"]
    # 
    cutoff=0.05/num_samples

    for distribution, dist_name in zip(distributions, distribution_names):
        params = fit_data_to_dist(NULL_list, distribution)
        print(f">>> {dist_name}")
        for theta in [1,0.75,0.5]:
            outlist = sub_out_path+"/null-"+str(numNulls)+"_sample-"+str(num_samples)+"_t"+str(theta)+"-"+str(dist_name)+"p_list.csv"
            if os.path.exists(outlist):
                dist_p_values,em_p_values,NS_p_list,MS_p_list = read_csv(outlist)
            else:
                value_list,NS_p_list,MS_p_list,em_p_values = sim_alts(out_path,num_samples,model,numhets,readperhet,theta=theta,sigma=sigma,value_name=value_name,lambdas_choice_gam=1.2,NULL_list=NULL_list)
                dist_p_values = [x for x in distribution.sf(value_list, *params)]
            count_b = sum(1 for value in dist_p_values if value <= cutoff)
            count_em = sum(1 for value in em_p_values if value <= cutoff)
            count_n = sum(1 for value in NS_p_list if value <= cutoff)
            count_m = sum(1 for value in MS_p_list if value <= cutoff)
            print(">>")
            print(f"theta {theta:.2f} -- {dist_name} distr p-values: {count_b/num_samples*100:.5f}% <= {cutoff}")
            print(f"theta {theta:.2f} -- empiri distr p-values: {count_em/num_samples*100:.5f}% <= {cutoff}")
            #print(f"theta {theta:.2f} -- empirical    p-values: {count_em/num_samples*100:.2f}% <= {cutoff}")
            print(f"theta {theta:.2f} -- NS           p-values: {count_n/num_samples*100:.5f}% <= {cutoff}")
            print(f"theta {theta:.2f} -- MS           p-values: {count_m/num_samples*100:.5f}% <= {cutoff}")
            if save is True and not os.path.exists(outlist):
                outlist = sub_out_path+"/null-"+str(numNulls)+"_sample-"+str(num_samples)+"_t"+str(theta)+"-"+str(dist_name)+"p_list.csv"
                write_p_list(outlist,dist_p_values, em_p_values,NS_p_list,MS_p_list)

def write_p_list(outlist,p_list,em_list,NS_list,MS_list):
    # write p_list to csv file
    with open(outlist, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['dist_p_value', 'em_p_val', 'NS_p_val', 'MS_p_val'])
        for i in range(len(p_list)):
            writer.writerow([p_list[i],em_list[i],NS_list[i],MS_list[i]])

################################ main function
if __name__ == "__main__":
    # read in parameters
    if len(sys.argv) < 9:
        print("Error: Not enough arguments provided.")
        print("Usage: python calculate_p.py <output directory> <model path> <number of hets> <reads per het> <number of samples> <number of nulls>")
        sys.exit(1)
    else:
        out_dir = sys.argv[1]
        model_path = sys.argv[2]
        numhets = int(sys.argv[3])
        readperhet = int(sys.argv[4])
        num_samples = int(sys.argv[5])
        numNulls = int(sys.argv[6])
        save = True if sys.argv[7].lower() == "true" else False
        value = sys.argv[8]
        print(f"Log: simulating {numNulls} null genes with {numhets} hets and {readperhet} read depth per het")
        print(f"Log: calculating p-values for {num_samples} simulated genes with {numhets} hets and {readperhet} read depth per het @ thetas 1, 0.75, 0.5")
    value_names = ["mean_over_std",'abs_mean_over_std']
    if value not in value_names:
        print(f"Error: {value} not in {value_names}")
        sys.exit(1)
    # Update the output directory path
    out_dir = out_dir + "/M" + str(numhets) + "_D" + str(readperhet)

    # Check if output directory exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"... directory {out_dir} created.")
    else:
        print(f"... directory {out_dir} already exists.")


    # pre-define variables 
    lambdas_choice_gam = 1.2
    sigma = 0.7

    # read in model
    simulate_null_ans_calP(out_dir,numNulls, num_samples, model_path, sigma,lambdas_choice_gam, numhets, readperhet, value_name=value,save=save)