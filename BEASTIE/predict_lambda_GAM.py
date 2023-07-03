#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import argparse
import logging
import multiprocessing
import pandas as pd
import numpy as np


def inv_logit(p):
    return np.exp(p) / (1 + np.exp(p))


def logit(p):
    return np.log(p) - np.log(1 - p)


def get_lambda_from_gam_pre_partition(
    model, hets, totalcount, expected_type1error, candidate_log_lambdas
):
    logit_expected_type1error = logit(expected_type1error)

    INITIAL_PREDICTION_COUNT = int(np.ceil(np.sqrt(len(candidate_log_lambdas))))

    initial_is = np.linspace(
        0, len(candidate_log_lambdas) - 1, INITIAL_PREDICTION_COUNT
    ).astype(int)
    initial_log_lambdas = candidate_log_lambdas[initial_is]
    initial_predictions = model.predict(
        [[hets, totalcount, lam] for lam in initial_log_lambdas]
    )

    min_i = None
    for i in range(1, INITIAL_PREDICTION_COUNT):
        if initial_predictions[i] < logit_expected_type1error:
            min_i = initial_is[i - 1] if i > 0 else 0
            max_i = initial_is[i]
            break

    if min_i is not None:
        partitioned_candidate_log_lambdas = candidate_log_lambdas[min_i : max_i + 1]
        return get_lambda_from_gam(
            model,
            hets,
            totalcount,
            expected_type1error,
            partitioned_candidate_log_lambdas,
        )
    else:
        return np.log(3)


def get_lambda_from_gam(
    model, hets, totalcount, expected_type1error, candidate_log_lambdas
):
    # prepare input
    data = [[hets, totalcount, lam] for lam in candidate_log_lambdas]

    # prediction
    #prediction = inv_logit(model.predict(data))
    prediction = model.predict(data)
    chosen_lambda = 3
    if min(10**(prediction)) <= expected_type1error:
        chosen_lambda = 10**(
            data[np.where(10**(prediction) <= expected_type1error)[0][0]][2]
        )
    return chosen_lambda


def predict_lambda_onrealdata(
    expected_type1error, in_filename, out_filename, modelByName
):
    in_data = pd.read_csv(in_filename, sep="\t")
    candidate_log_lambdas = np.log(np.linspace(1, 3, 3000))

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
                        candidate_log_lambdas,
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
