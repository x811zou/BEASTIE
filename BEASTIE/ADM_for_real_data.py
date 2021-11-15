#!/usr/bin/env python
import logging
import os
import os.path
import pickle
import statistics

import numpy as np
import pandas as pd
from scipy.stats import percentileofscore


def AA_estimate(A,R):
    AA = abs(A-R)/(A+R)
    return AA

def fixed_simulator(D,M,N,estimate):
    esti = []
    for k in range(N):
        AR = []
        for i in range(M):
            A = np.random.binomial(D[i], 0.5)
            R = D[i] - A
            AR.append(AA_estimate(A,R))

        esti.append(statistics.mean(AR))
    # find the tail
    pval = 1- percentileofscore(esti,estimate)/100
    #print("pvalue: "+str(pval))
    return pval

def getAA(fields):
    if(len(fields)>=4):
        esti=[]
        Mreps=int(fields[1])
        #Mreps = hets
        #p = float(theta)/(float(theta)+1)
        depth=[]
        for rep in range(Mreps):
            A = float(fields[2+rep*2])
            R = float(fields[3+rep*2])
            estimate = abs(R-A)/(A+R)
            esti.append(estimate)
            total_counts=A+R
            depth.append(int(total_counts))
        avg_esti = statistics.mean(esti)
        hets = Mreps
        pval = fixed_simulator(depth,hets,1000,avg_esti)
    return round(avg_esti,3), round(pval,3)

def run(prefix,inFile,out,picklename):
    out_path=out+"/output_pkl/ADM/"
    if not os.path.isdir(out_path):os.mkdir(out_path)
    out_path_theta = out_path+"AA_esti/"
    out_path_pval = out_path+ "AA_pval/"
    if not os.path.isdir(out_path_theta):os.mkdir(out_path_theta)
    if not os.path.isdir(out_path_pval):os.mkdir(out_path_pval)
    esti_list=[]
    pval_list = []
    gene=[]
    filename=os.path.basename(inFile)
    out1 = out_path_theta + picklename
    out2 = out_path_pval + picklename

    with open(inFile,"rt") as IN:
        for line in IN:
            # logging.info('{0}'.format(line))
            fields=line.rstrip().split()
            ID=fields[0]
            gene.append(ID)
            # logging.info('gene {0}'.format(ID))
            avg_esti,pval = getAA(fields)
            esti_list.append(avg_esti)
            pval_list.append(pval)
    ADM_df = pd.DataFrame(np.column_stack([gene, esti_list, pval_list]),
                               columns=['geneID', 'ADM_esti', 'ADM_pval'])
    ADM_df.to_csv(out+"/"+prefix+"_ASE_ADM.tsv",sep="\t",header=True,index=False)
    logging.info('..... ADM estimates size {0}'.format(len(esti_list)))
    logging.info('..... ADM p val size {0}'.format(len(pval_list)))

    pickle.dump(esti_list,open(out1,'wb'))
    pickle.dump(pval_list,open(out2,'wb'))
    return ADM_df
