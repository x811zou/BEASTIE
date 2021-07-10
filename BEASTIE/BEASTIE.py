#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
import os
import sys
import os.path
from numpy.core.numeric import outer
import logging
import argparse
import pandas as pd
from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe
import beastie_step1
import beastie_step2

def _build(args):
    logname=args.in_path+"log."+args.prefix
    if os.path.isfile(logname):
        os.remove(logname)
    logging.basicConfig(filename=logname,
                            filemode='a',
                            format='%(asctime)-15s [%(levelname)s] %(message)s',
                            level=logging.DEBUG)
    logging.info("Running BEASTIE")
    #logging.basicConfig(filename=logname,level=logging.INFO, format='%(asctime)-15s [%(levelname)s] %(message)s')

    logging.info('>>>>>>>>>>>>>>>>>>>> ')
    logging.info('>>>>>>>>>>>>>>>>>>>> step1: Processing raw data & annotating LD and AF information')
    logging.info('>>>>>>>>>>>>>>>>>>>> ')
    hetSNP_intersect_unique,meta,hetSNP_intersect_unique_forlambda_file,hetSNP_intersect_unique_lambdaPredicted_file = beastie_step1.run(args.prefix, args.vcf_sample_name, args.in_path,args.out,args.model,args.vcf,args.ref_dir, args.ancestry, args.chr_start, args.chr_end, args.min_total_cov,args.min_single_cov,args.read_length,args.LD_token,args.pileup,args.hetSNP,args.hetSNP_intersect_unique,args.parsed_pileup,args.meta,args.alpha)          
    logging.info('>>>>>>>>>>>>>>>>>>>> ')                
    logging.info('>>>>>>>>>>>>>>>>>>>> step2: Preparing input in a format required for BEASTIE model')
    logging.info('>>>>>>>>>>>>>>>>>>>> ')
    #print("alpha is %s"%(args.alpha))
    beastie_step2.run(hetSNP_intersect_unique,meta,hetSNP_intersect_unique_forlambda_file,hetSNP_intersect_unique_lambdaPredicted_file,args.prefix,args.read_length,args.min_total_cov,args.min_single_cov,args.alpha,args.model,args.sigma,args.in_path,args.out,args.cutoff,args.SAVE_INT)

def main():
    parser = argparse.ArgumentParser(description='Utilities for creating and working with BEASTIE.')
    subparsers = parser.add_subparsers(title='commands')
    build_parser = subparsers.add_parser('build', help='Process pileup data and VCF data.')
    # required arguments
    build_parser.add_argument('--prefix', required=True, help='sample name to use generate output files.')
    build_parser.add_argument('--vcf_sample_name', required=True, help='column name for sample in VCF file.')
    build_parser.add_argument('--ref_dir', required=True, help='reference directory containing AF file and gencode folder.')
    build_parser.add_argument('--vcf', required=True, help='VCF file to use.')
    build_parser.add_argument('--pileup', required=True, help='pileup file to use.')
    build_parser.add_argument('--in_path', required=True, help='path_to_input file to use.')
    build_parser.add_argument('--ancestry', required=True, help='example provided: EUR to HG00096.')
    build_parser.add_argument('--chr_start', required=True, help='Start chromosome')
    build_parser.add_argument('--chr_end', required=True, help='End chromosome')
    build_parser.add_argument('--read_length', required=True, help='length of RNA reads')
    build_parser.add_argument('--LD_token', required=True, help='Registered token to do LD annotation')
    # optional files that users can generate themselves
    build_parser.add_argument('--hetSNP', help='The file containing het SNPs to use. default name is ${prefix}_hetSNP.tsv')
    build_parser.add_argument('--hetSNP_intersect_unique', help='The file containing intersected between hetSNP and parsed_pileup data. One reads only count once. default name is ${prefix}_hetSNP_intersect_unique.tsv')
    build_parser.add_argument('--parsed_pileup', help='The parsed pileup file. default name is ${prefix}_parsed_pileup.tsv')
    build_parser.add_argument('--meta', help='The meta file containing information for the logistic model. default name is ${prefix}_logisticReg_input.tsv')
    # default setting
    build_parser.add_argument('--model', default='stan_path/iBEASTIE2/iBEASTIE2',help='model name.')
    build_parser.add_argument('--min_total_cov', type=int, default=1, help='Minimal total coverage required to filter varaints to be included in the graph. Defaults to 1.')
    build_parser.add_argument('--min_single_cov', type=int, default=0, help='Minimal coverage for one allele required to filter variants to be included in the graph. Defaults to 0.')
    build_parser.add_argument('--cutoff', type=float, default=0.5, help='Significance cutoff for ASE. Defaults to 0.5.')
    build_parser.add_argument('--alpha', type=float, default=0.5, help='Significance cutoff for ASE. Defaults to 0.5.')
    build_parser.add_argument('--sigma', type=float, default=0.5, help='Significance cutoff for ASE. Defaults to 0.5.')
    build_parser.add_argument('--out', default='output', help='Location to write output files. Defaults to \'output\' inside the sample folder')
    build_parser.add_argument('--SAVE_INT', default='False', help='Whether to save intermediate output. Defaults to remove \'TEMP\' folder inside \'output\' folder')

    build_parser.set_defaults(command=_build)

    args = parser.parse_args()
    args.command(args)

if __name__ == '__main__':
    main()
