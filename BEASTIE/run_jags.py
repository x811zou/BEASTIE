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

# from collections import Mapping

## replicate what jags is doing
np.random.seed(0)
np.set_printoptions(precision=1)


def genotype_bugs_model(
    smaller_sum, larger_sum, total_count, mean, variance, WARMUP=1000, KEEPER=1000
):
    code = """
    model {
        p ~ dbeta(alpha+1,beta+1);
        x ~ dbin(p,n);
        n ~ dnegbin(1-(Var-mu)/Var,mu*mu/(Var-mu))I(1,1000000);
    }
    """
    model = pyjags.Model(
        code,
        data=dict(
            alpha=smaller_sum,
            beta=larger_sum,
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
        pval = summary(samples, varname, total_count)
    return pval


def summary(samples, varname, total_count):
    values = samples[varname]
    N = values.shape[1]
    pval = np.sum(values[0, :, 0] >= total_count) / values.shape[1]
    return pval


# ===========================================
# command line: ./run_jags.py 49 51 5 31 300
# ===========================================
parser = argparse.ArgumentParser()
parser.add_argument(
    "smaller_sum",
    help="prior count for smaller (come from the other (nonzero) sites, combined via pseudo-phasing)",
)
parser.add_argument(
    "larger_sum",
    help="prior count for larger (come from the other (nonzero) sites, combined via pseudo-phasing)",
)
parser.add_argument(
    "nonzero_count", help="the total count at the zero site (non-zero allele count)"
)
parser.add_argument("mu", help="mu (mean of the empirical distribution)")
parser.add_argument("var", help="variance (variance of the empirical distribution)")
parser.add_argument(
    "--n_warmup", help="number of warm up samples (burn out)", default=1000
)
parser.add_argument("--n_keep", help="number of samples to keep", default=1000)
args = parser.parse_args()
alpha = args.smaller_sum
beta = args.larger_sum
n = args.nonzero_count
n = int(n)
mu = args.mu
var = args.var
WARMUP = args.n_warmup
KEEPER = args.n_keep

pval = genotype_bugs_model(alpha, beta, n, mu, var, WARMUP, KEEPER)
print(pval)
