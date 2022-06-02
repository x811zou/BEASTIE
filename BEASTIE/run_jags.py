#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

## import python packages
import pyjags
import numpy as np
from pyjags.model import *
from pyjags.modules import *
import sys
import argparse
import math
from pyrsistent import v

# from collections import Mapping

## replicate what jags is doing
np.random.seed(0)
np.set_printoptions(precision=1)


def genotype_bugs_model(
    smaller_sum,
    total_sum,
    nonzero_count,
    mean,
    variance,
    WARMUP=1000,
    KEEPER=1000,
    RE_ITERATE=False,
    pcutoff_low=None,
    pcutoff_high=None,
    N=None,
):
    # p is between 0 and 1, the probablity of having ALT allele
    # alpha: sum of lower counts (ALT), beta: sum of total counts
    # X is the predicted ALT allelel count under binomial distribution based on n and p
    # samples lots of values for n (the non-zero allele count)
    # what proportion of the predicted values are greater than observed n
    # for the SNP with 0 allele count, 0 is the ALT count, n is the REF count, n is also the total count
    # based on the p (ALT allele ratio) we calculated from other SNPs without 0 count, we use that to fit a binomial distribution and predict possible counts for ALT allele on the SNPs with one zero count.
    # Under binomial distribution, we have 1000 estimates of possible ALT counts. the observed ALT allele count is 0. p-val is the proportion of 1000 estimates that are smaller than or equal to observed ALT allele count?
    code = """
    model {
        p ~ dunif(0,1);                                                              
        alpha ~ dbinom(p,beta)T(1,beta-1)                                            
        x ~ dbin(p,n);                                                               
        n ~ dnegbin(1-(Var-mu)/Var,mu*mu/(Var-mu))I(1,1000000);   
    }
    """
    model = pyjags.Model(
        code,
        data=dict(
            alpha=smaller_sum,
            beta=total_sum,
            mu=mean,
            Var=variance,
            x=0,
        ),
        init=dict(
            p=0.5,
            n=5,
        ),
        chains=1,
        adapt=WARMUP,
        progress_bar=False,
    )
    samples = model.sample(KEEPER, vars=["n"])
    for varname in ["n"]:
        pval = summary(samples, varname, nonzero_count)
    if RE_ITERATE:
        if (pval <= pcutoff_high) and (pval >= pcutoff_low):
            samples = model.sample(N, vars=["n"])
            for varname in ["n"]:
                new_pval = summary(samples, varname, nonzero_count)
        else:
            new_pval = math.nan
    else:
        new_pval = math.nan
    return pval, new_pval


def summary(samples, varname, nonzero_count):
    values = samples[varname]
    N = values.shape[1]
    pval = np.sum(values[0, :, 0] >= nonzero_count) / values.shape[1]
    # pval = np.sum(values[0, :, 0] == 0) / values.shape[1]
    return pval


# ===========================================
# command line: ./run_jags.py 49 100 5 31 300 10000
# ===========================================


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "smaller_sum",
        help="prior count for smaller (come from the other (nonzero) sites, combined via pseudo-phasing)",
    )
    parser.add_argument(
        "total_sum",
        help="prior count for total (come from the other (nonzero) sites, combined via pseudo-phasing)",
    )
    parser.add_argument(
        "nonzero_count", help="the total count at the zero site (non-zero allele count)"
    )
    parser.add_argument("mu", help="mu (mean of the empirical distribution)")
    parser.add_argument("var", help="variance (variance of the empirical distribution)")
    parser.add_argument("--rep", help="number of replication", default=1)
    parser.add_argument(
        "--n_warmup", help="number of warm up samples (burn out)", default=1000
    )
    parser.add_argument("--n_keep", help="number of samples to keep", default=1000)
    args = parser.parse_args()
    alpha = args.smaller_sum
    beta = args.total_sum
    n = args.nonzero_count
    n = int(n)
    mu = args.mu
    var = args.var
    rep = int(args.rep)
    WARMUP = int(args.n_warmup)
    KEEPER = int(args.n_keep)
    # if rep > 1:
    #     pval_list=[]
    for i in range(rep):
        pval = genotype_bugs_model(alpha, beta, n, mu, var, WARMUP, KEEPER)
        print(pval)
    #         pval_list.append(pval)
    #     #print(pval_list)
    # else:
    #     pval = genotype_bugs_model(alpha, beta, n, mu, var, WARMUP, KEEPER)
    #     print(pval)
