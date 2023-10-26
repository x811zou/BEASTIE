#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import argparse
import logging
import multiprocessing
import pandas as pd
import numpy as np
from scipy.special import expit
import sys
def logit(p):
    return np.log(p) - np.log(1 - p)

def get_lambda_from_gam_pre_partition(
    model, hets, totalcount, expected_type1error, candidate_lambdas
):
    INITIAL_PREDICTION_COUNT = int(np.ceil(np.sqrt(len(candidate_lambdas))))
    log_hets=np.log(hets)
    log_totalcount=np.log(totalcount)

    initial_is = np.linspace(
        0, len(candidate_lambdas) - 1, INITIAL_PREDICTION_COUNT
    ).astype(int)
    initial_lambdas = candidate_lambdas[initial_is]
    initial_predictions = model.predict(
        [[log_hets, log_totalcount, np.log(lam-1+0.001)] for lam in initial_lambdas]
    )

    min_i = None
    for i in range(1, INITIAL_PREDICTION_COUNT):
        if initial_predictions[i] < logit(expected_type1error):
            min_i = initial_is[i - 1] if i > 0 else 0
            max_i = initial_is[i]
            break

    if min_i is not None:
        partitioned_candidate_lambdas = candidate_lambdas[min_i : max_i + 1]
        return get_lambda_from_gam(
            model,
            log_hets,
            log_totalcount,
            expected_type1error,
            partitioned_candidate_lambdas,
        )
    else:
        return 3


def get_lambda_from_gam(
    model, log_hets, log_totalcount, expected_type1error, candidate_lambdas
):
    # prepare input
    data = [[float(log_hets), float(log_totalcount), np.log(lam-1+0.001)] for lam in candidate_lambdas]

    # prediction
    prediction = model.predict(data)
    chosen_lambda = 3
    #print(f"hets: %s totalcount: %s" % (hets, totalcount))
    #print(prediction)
    #print(expit(prediction))
    if min(expit(prediction)) <= expected_type1error:
        chosen_lambda=np.exp(data[np.where(expit(prediction) <= expected_type1error)[0][0]][2])+1
    return round(chosen_lambda,10)


def predict_lambda_onrealdata(
    expected_type1error, in_filename, out_filename, modelByName
):
    in_data = pd.read_csv(in_filename, sep="\t")
    candidate_lambdas = np.linspace(1, 3, 2000)
    # for model_name, model in modelByName.items():
    #     in_data[model_name] = in_data.apply(lambda x: get_lambda_from_gam(model, x["number.of.hets"], x["totalCount"], expected_type1error, candidate_lambdas),axis=1,)
    with multiprocessing.Pool() as pool:
        for model_name, model in modelByName.items():
            column_data = pool.starmap(
                get_lambda_from_gam_pre_partition,
                [
                    (
                        model,
                        x["number.of.hets"],
                        x["totalCount"],
                        expected_type1error,
                        candidate_lambdas,
                    )
                    for _, x in in_data.iterrows()
                ],
            )
            in_data[model_name] = column_data

    in_data.to_csv(out_filename, sep="\t")
    logging.info(f"model input with GAM predicted lambda saved to {out_filename}")


def main():
    #################################################################
    # srun --mem=5000 --pty bash
    # python /hpc/group/allenlab/scarlett/python/stan_model/run_BEASTIE_compute_type1error_counts.py /hpc/group/allenlab/scarlett/ASE_simulation/iBEASTIE_empirical_data/g-1000/g-1000_h-8_d-5_t-0.5_semi1.txt iBEASTIE2 1 3 0.001 test.tsv 1
    ################################################################# read in arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--adjusted_alpha", type=float, required=True)
    parser.add_argument("--out_path", type=str, required=True)
    parser.add_argument("--prefix", type=str, required=True)
    parser.add_argument("--model", type=str, required=True)
    parser.add_argument("--in_path", type=str, required=True)
    parser.add_argument("--filename", type=str, required=True)
    args = parser.parse_args()

    ################################################################# parameters
    adjusted_alpha = args.adjusted_alpha
    tmp_path = args.out_path
    sample = args.prefix
    gam_model = args.model
    in_path = args.in_path
    file_for_lambda = args.filename
    ################################################################# parameters
    out_filename = tmp_path + "/" + file_for_lambda
    in_filename = in_path + "/" + file_for_lambda

    ################################################################# predicta lambda
    predict_lambda_onrealdata(adjusted_alpha, in_filename, out_filename, gam_model)


if __name__ == "__main__":
    main()
