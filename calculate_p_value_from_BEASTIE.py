import multiprocessing
import os
import sys
import subprocess
import tempfile
import shutil
import pandas as pd
from scipy.stats import norm
import numpy as np
from scipy import stats
from datetime import datetime, timedelta
import time
import logging
from multiprocessing import Pool
from numpy.random import uniform
from scipy.stats import t
from scipy.stats import rv_continuous
from scipy.optimize import minimize
from scipy.special import gamma
from scipy.stats import skew
from BEASTIE.run_model_stan_wrapper import save_raw_theta_parallel, summarize, runModel
import pickle
from scipy.stats import t, skewnorm
from scipy.optimize import minimize
# Configure logging
logging.basicConfig(level=logging.INFO)

"""
python calculate_p_value_from_BEASTIE.py /home/scarlett/github/RNAseq-analysis/run_quickBeast/test_data/test /home/scarlett/github/RNAseq-analysis/run_quickBeast/test_data/test_BEASTIE3.out /home/scarlett/github/RNAseq-analysis/stan_models/BEASTIE3-pi0.05 0.7
python calculate_p_value_from_BEASTIE.py /home/scarlett/github/RNAseq-analysis/run_quickBeast/test_data/test /data2/stan/BEASTIE3-pi0.05/sigma0.7/parametrized/ASE_0.05_error/tmp_output.out /home/scarlett/github/RNAseq-analysis/stan_models/BEASTIE3-pi0.05 0.7
""" 

def simulate_null_genes_helper(args):
    geneID, num_hets, average_count, model_path,sigma, WARMUP, KEEPER, phasing_method = args
    #print(f">>>>> geneID: {geneID}, num_hets: {num_hets}, total_count: {total_count}")
    mean, std, n_loc, n_scale, t_df, t_loc, t_scale, st_df, st_loc, st_scale = simulate_null_genes(num_hets, average_count, model_path,sigma, WARMUP, KEEPER, phasing_method)
    return geneID, mean, std, n_loc, n_scale, t_df, t_loc, t_scale, st_df, st_loc, st_scale

def calculate_time(start_t, end_t):
    elapsed_time_seconds = end_t - start_t
    elapsed_time_formatted = str(timedelta(seconds=elapsed_time_seconds))
    return elapsed_time_formatted

def generate_fields(geneID,M, D, theta,switching_error = 0.05):
    # calculate probability for binomial distribution
    p = theta / (1.0 + theta)
    # calculate alternative and reference read counts for each het
    alt_counts = np.random.binomial(D, p, M)
    ref_counts = D - alt_counts
    # construct output fields
    switched = False
    fields = [geneID, M]
    
    for counts in zip(alt_counts, ref_counts):
        if switched:
            fields.append(counts[1])
            fields.append(counts[0])
        else:
            fields.append(counts[0])
            fields.append(counts[1])
        if np.random.uniform() <= switching_error:
            switched = not switched
    
    # Adding number of missing pi and phasing error fields
    fields.append(0)  # Number of missing pi
    fields += [switching_error] * (M - 1)  # Phasing error (We assume all pairs to be -1)
    # Convert all items to strings
    # string_list = [str(item) for item in fields]
    # tab_separated_string = "\t".join(string_list)

    return fields

def simulate_null_genes(number_of_hets, average_read_depth_per_het, model_path, sigma, WARMUP=1000, KEEPER=1000, phasing_method="VCF", pi=0.05):
    #(f"het-count: {number_of_hets}, read-depth: {average_read_depth_per_het}")
    # SIMULATE GENE INPUTS
    simulated_genes = []
    for i in range(1000):
        gene_id = f"gene_{i}"
        field = generate_fields(gene_id, M=number_of_hets, D=average_read_depth_per_het, theta=1, switching_error = pi)
        #print(field)
        simulated_genes.append(field)
    
    # RUN qb simulated genes 
    model_thetas = save_raw_theta_parallel_list(
            simulated_genes,
            model_path,
            sigma,
            WARMUP,
            KEEPER,
            phasing_method,
        )
    
    simulated_gene_results = parse_stan_output_pval(simulated_genes, model_thetas)
    
    # compute summary statistics
    mean = simulated_gene_results["BEASTIE_zscore"].mean()
    std = simulated_gene_results["BEASTIE_zscore"].std()
    z_scores = simulated_gene_results["BEASTIE_zscore"].tolist()
    # df: degrees of freedom
    # loc: location
    # scale: scale of std
    n_loc, n_scale = stats.norm.fit(z_scores)
    t_df, t_loc, t_scale = stats.t.fit(z_scores)
    st_df, st_loc, st_scale = stats.skewnorm.fit(z_scores)

    return mean, std, n_loc, n_scale, t_df, t_loc, t_scale, st_df, st_loc, st_scale


def parse_stan_output_pval_initializer(thetas):
    global g_thetas
    g_thetas = thetas

def parse_stan_output_pval_worker(fields):
    global g_thetas
    #fields = line.rstrip().split("\t")
    gene_id = fields[0]
    gene_thetas = g_thetas.get(gene_id)

    if not gene_thetas:
        #print(fields)
        #print(gene_id, g_thetas.keys())
        return None
    
    (
        mean,
        median,
        variance,
        left_CI,
        right_CI,
        mad,
        log2_mean,
        log2_median,
        log2_variance,
        abslog2_mean,
        abslog2_median,
        abslog2_variance,
    ) = summarize(gene_thetas, 0.05)

    return (
        gene_id,
        round(median, 3),
        round(mean, 3),
        round(variance, 3),
        round(left_CI, 3),
        round(right_CI, 3),
        log2_median,
        log2_mean,
        log2_variance,
        abslog2_median,
        abslog2_mean,
        abslog2_variance,
    )

def parse_stan_output_pval(genes, thetas):
    # Assuming `thetas` is already defined somewhere
    initializer_args = (thetas,)

    with multiprocessing.Pool(initializer=parse_stan_output_pval_initializer, initargs=initializer_args) as pool:
        # The worker function now takes each list from `genes` as an argument
        rows = pool.map(parse_stan_output_pval_worker, genes)

    #print(rows)
    df = pd.DataFrame(
        rows,
        # columns=[
        #     "geneID",
        #     "posterior_median",
        #     "posterior_mean",
        #     "posterior_variance",
        #     "CI_left",
        #     "CI_right",
        #     "log2_posterior_median",
        #     "log2_posterior_mean",
        #     "log2_posterior_variance",
        #     "abslog2_posterior_median",
        #     "abslog2_posterior_mean",
        #     "abslog2_posterior_variance",
        # ],
    )
    column_names = [
            "geneID",
            "posterior_median",
            "posterior_mean",
            "posterior_variance",
            "CI_left",
            "CI_right",
            "log2_posterior_median",
            "log2_posterior_mean",
            "log2_posterior_variance",
            "abslog2_posterior_median",
            "abslog2_posterior_mean",
            "abslog2_posterior_variance",
        ]
    df.columns = column_names
    df[["log2_posterior_std"]] = np.sqrt(df[["log2_posterior_variance"]])
    df["BEASTIE_zscore"] = df["log2_posterior_mean"] / df["log2_posterior_std"]
    gene_df = df[["geneID", "log2_posterior_mean", "log2_posterior_variance","BEASTIE_zscore"]]
    return gene_df

def calculate_norm_2sided_pvalue(loc, scale, sample_X):
    # Get the cumulative probability of X (left tail)
    left_tail_p = stats.norm.cdf(sample_X, loc, scale)
    
    # Get the survival function value for X (right tail)
    right_tail_p = stats.norm.sf(sample_X, loc, scale)
    
    # Return the double-sided p-value
    return 2 * min(left_tail_p, right_tail_p)

def calculate_t_2sided_pvalue(df, loc, scale, sample_X):
    # Get the cumulative probability of X (left tail)
    left_tail_p = stats.t.cdf(sample_X, df, loc, scale)
    
    # Get the survival function value for X (right tail)
    right_tail_p = stats.t.sf(sample_X, df, loc, scale)
    
    # Return the double-sided p-value
    return 2 * min(left_tail_p, right_tail_p)

def calculate_st_2sided_pvalue(a, loc, scale, sample_X):
    # Get the cumulative probability of X (left tail)
    left_tail_p = stats.skewnorm.cdf(sample_X, a, loc, scale)
    
    # Get the survival function value for X (right tail)
    right_tail_p = stats.skewnorm.sf(sample_X, a, loc, scale)
    
    # Return the double-sided p-value
    return 2 * min(left_tail_p, right_tail_p)

def write_genes_to_BEASTIE_input_file(genes, file_path):
    with open(file_path, 'w') as file:
        for gene in genes:
            file.write(f'{gene}')
            file.write("\n")

def run_model_worker_initializer(
    tmp_dir,
    model,
    sigma,
    WARMUP,
    KEEPER,
    phasing_method,
):
    global g_model, g_sigma, g_WARMUP, g_KEEPER, g_phasing_method
    g_model = model
    g_sigma = sigma
    g_WARMUP = WARMUP
    g_KEEPER = KEEPER
    g_phasing_method = phasing_method

    global g_tmp_output_file, g_stan_output_file, g_init_file
    g_tmp_output_file = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False).name
    g_stan_output_file = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False).name
    g_init_file = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False).name

def run_model_worker(fields):
    global g_model, g_tmp_output_file, g_stan_output_file, g_init_file, g_sigma, g_WARMUP, g_KEEPER, g_phasing_method
    geneID, thetas = runModel(
        g_model,
        fields,
        g_tmp_output_file,
        g_stan_output_file,
        g_init_file,
        g_sigma,
        g_WARMUP,
        g_KEEPER,
        g_phasing_method,
    )
    return geneID, thetas


def save_raw_theta_parallel_list(
    genes,
    model,
    sigma,
    WARMUP,
    KEEPER,
    phasing_method,
):
    # run BEASTIE
    processes = os.cpu_count()
    with tempfile.TemporaryDirectory() as tmp_dir, multiprocessing.Pool(
        processes = processes,
        initializer = run_model_worker_initializer,
        initargs=(tmp_dir, model, sigma, WARMUP, KEEPER, phasing_method),
    ) as pool:
        items = []
        for item in pool.imap_unordered(run_model_worker, genes):
            items.append(item)
            if len(items) % 500 == 0:
                logging.info(f".... Processed thetas for {len(items)} genes")
        #logging.info(".... Processed thetas for all genes")

    model_theta = dict(items)
    return model_theta

def main():
    inFile = sys.argv[1]
    output_file_path = sys.argv[2]
    model_path = sys.argv[3]
    sigma = sys.argv[4]
    WARMUP = 1000
    KEEPER = 1000
    phasing_method = "nophasing"

    # define output file name
    out_path = os.path.dirname(output_file_path)
    out_path = os.path.join(out_path, "output_pkl")
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        # outname is the base file name of input_file_path without postix ".txt" + outname
    outname = f"_s-{sigma}.pickle"
    outname = os.path.basename(inFile).replace(".txt", "")  + outname
    thetas_file = os.path.join(out_path, outname)
    start_t = time.time()

    input_genes = []
    with open(inFile) as file:
        for line in file.readlines():
            fields = line.strip().split("\t")
            input_genes.append(fields)

    # step1: run BEASTIE on input genes, and then parsing input
    if os.path.isfile(thetas_file):
        logging.info("...... Already finshed running {0} and saved raw theta at : {1}".format(
            os.path.basename(model_path), thetas_file
        ))
    else:
        model_thetas = save_raw_theta_parallel_list(
            input_genes,
            model_path,
            sigma,
            WARMUP,
            KEEPER,
            phasing_method,
        )
        pickle.dump(model_thetas, open(thetas_file, "wb"))
        beastie_end_t = time.time()
        logging.info(
            "...... Just finshed running {0} and save raw theta at : {1}".format(
                os.path.basename(model_path), thetas_file
            )
        )
        logging.info(f"...... Running BEASTIE on input used ${calculate_time(start_t, beastie_end_t)}")

    # step2: parsing
    logging.info("...... Start parse_stan_output")
    thetas = pickle.load(open(thetas_file, "rb"))
    gene_df = parse_stan_output_pval(input_genes, thetas)

    parse_end_t = time.time()
    logging.info(f"...... Just finished parsing input in ${calculate_time(beastie_end_t,parse_end_t)}")

    # step3: NULL simulation
    null_simulation_cache = {}
    null_simulation_data = []
    runs_completed = 0

    all_runs = []
    for gene in input_genes:
        try:
            gene_id = gene[0]
            count = int(gene[1])
            sum_of_values = sum(map(int, gene[2:(count*2+1)]))  # Convert items to int before summing
            all_runs.append((gene_id, count, sum_of_values))
        except ValueError:
            # Handle the error or report it
            print(f"Skipping gene {gene[0]} due to invalid number format.")

    unique_runs = set([(f"h-{r[1]}_d-{int(r[2]/r[1])}", r[1], int(r[2]/r[1]), model_path, sigma, WARMUP, KEEPER, phasing_method) for r in all_runs])
    logging.info(f"...... {datetime.now()} Starting {len(unique_runs)} unique gene simulation runs {len(all_runs)} total runs")

    for params in unique_runs:
        # Call the simulation function directly with the parameters unpacked
        output = simulate_null_genes_helper(params)

        #print(output[0])
        #print(f"{output[0]} : {output[1]} {output[2]} {output[3]} {output[4]} {output[5]}")

        null_simulation_cache[output[0]] = output[1], output[2] , output[3] , output[4], output[5], output[6], output[7] , output[8] , output[9], output[10]

        runs_completed += 1

        if runs_completed % 10 == 0:
            genes_per_second = runs_completed / (time.time() - parse_end_t)
            seconds_remaining = (len(unique_runs) - runs_completed) / genes_per_second
            logging.info(f"{datetime.now()} Completed {runs_completed} gene simulations ({round(genes_per_second, 2)} genes/s) [ETA: {round(seconds_remaining, 1)}s]")

    for r in all_runs:
        #print(r)
        key = f"h-{r[1]}_d-{int(r[2]/r[1])}"
        mean, std, n_loc, n_scale, t_df, t_loc, t_scale, st_df, st_loc, st_scale = null_simulation_cache[key]
        null_simulation_data.append((r[0], mean, std, n_loc, n_scale, t_df, t_loc, t_scale, st_df, st_loc, st_scale))

    null_simulation_df = pd.DataFrame(null_simulation_data, columns=['geneID', 'null_zscore_mean', 'null_zscore_std' ,'n_loc' ,'n_scale','t_df' ,'t_loc' ,'t_scale','st_df' ,'st_loc' ,'st_scale'])
    
    #print(null_simulation_df)
    gene_df = pd.merge(gene_df, null_simulation_df, on="geneID")

    simulation_end_t = time.time()
    logging.info(f"...... Finished simulations in ${calculate_time(parse_end_t, simulation_end_t)}")

    gene_df['normal_p_value'] = gene_df.apply(lambda row: calculate_norm_2sided_pvalue(row['n_loc'], row['n_scale'], row['BEASTIE_zscore']), axis=1)
    gene_df['t_p_value'] = gene_df.apply(lambda row: calculate_t_2sided_pvalue(row['t_df'], row['t_loc'], row['t_scale'], row['BEASTIE_zscore']), axis=1)
    gene_df['st_p_value'] = gene_df.apply(lambda row: calculate_st_2sided_pvalue(row['st_df'], row['st_loc'], row['st_scale'], row['BEASTIE_zscore']), axis=1)
    columns_to_drop = ['null_zscore_mean', 'null_zscore_std','n_loc', 'n_scale','t_df','t_loc','t_scale','st_df','st_loc','st_scale']
    gene_df = gene_df.drop(columns=columns_to_drop)
    #(gene_df)
    # Save the DataFrame as a TSV file
    gene_df.to_csv(output_file_path, sep='\t', index=False)

    complete_time_t = time.time()
    logging.info(f"...... Completed BEASTIE on input {inFile} in ${calculate_time(start_t, complete_time_t)}")

# usage
if __name__ == "__main__":
    main()