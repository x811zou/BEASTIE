#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import argparse
import pandas as pd
import numpy as np


def inv_logit(p):
    return np.exp(p) / (1 + np.exp(p))


def logit(p):
    return np.log(p) - np.log(1 - p)


def get_lambda_from_gam(model, hets, totalcount, expected_type1error):
    # preapre input
    data = [[hets, totalcount, np.log(lam)] for lam in np.linspace(1, 3, 3000)]
    # prediction
    prediction = model.predict(data)
    chosen_lambda = 3
    # print(f"hets: %s totalcount: %s" % (hets, totalcount))
    if min(inv_logit(prediction)) <= expected_type1error:
        chosen_lambda = np.exp(
            data[np.where(inv_logit(prediction) <= expected_type1error)[0][0]][2]
        )
    # print(f"hets: %s totalcount: %s - lambda: %s" % (hets, totalcount, chosen_lambda))
    return chosen_lambda


def predict_lambda_onrealdata(
    expected_type1error, in_filename, out_filename, model3, model4
):
    in_data = pd.read_csv(in_filename, sep="\t")
    # print(in_data)
    in_data["gam3_lambda"] = in_data.apply(
        lambda x: get_lambda_from_gam(
            model3, x["number.of.hets"], x["totalCount"], expected_type1error
        ),
        axis=1,
    )
    in_data["gam4_lambda"] = in_data.apply(
        lambda x: get_lambda_from_gam(
            model4, x["number.of.hets"], x["totalCount"], expected_type1error
        ),
        axis=1,
    )
    # print(in_data)
    in_data.to_csv(out_filename, sep="\t")
    print(f"model input with GAM predicted lambda saved to %s" % str(out_filename))


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
