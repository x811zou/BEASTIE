#!/usr/bin/env python
from scipy.stats import binom
import os
from scipy import stats
import os.path
import pickle
import logging
import pandas as pd
import numpy as np
from pathlib import Path

# def create_binomial_library(depth):
#     count_p = {}
#     for i in range(int(depth)+1):
#         count_p[i] = binom.cdf(i, int(depth), 0.5)
#     return count_p

def getBaseline(fields,depth):
    if(len(fields)>=4):
	############
        #count_p = create_binomial_library(depth)	
        Mreps=int(fields[1])
        total_AR = 0
        for rep in range(Mreps):
            A = float(fields[2+rep*2])        #+1
            R = float(fields[3+rep*2])        #+1
            # logging.info('... {} site, A-{} R-{}'.format(rep,A,R))
            if A <= R:
                min_AR = A
                max_AR = R
            else:
                min_AR = R
                max_AR = A
            if total_AR < A+R:
                max_A=A
                max_R=R
                total_AR=A+R
            if str(rep) == "0":
                SS_AAR = A/(A+R)
                if SS_AAR == 1 :
                   SS_AAR = SS_AAR - 0.001 
                SS_esti = SS_AAR/(1-SS_AAR)
                # logging.info('... SS_esti - {}'.format(SS_esti))
                SS_prob = stats.binom_test(A, A+R, p=0.5, alternative='two-sided')
                # logging.info('... SS_prob - {}'.format(SS_prob))
            # esti3
            base3 = abs(0.5-max_AR/(A+R))
            # esti4
            diff_AR = float(abs(0.5-max_AR/(A+R)))
        MS_AAR = max_A/(max_A+max_R)
        if MS_AAR == 1 :
            MS_AAR = MS_AAR - 0.001 
        MS_esti = MS_AAR/(1-MS_AAR)
        # logging.info('... MS_esti - {}'.format(MS_esti)) 
        MS_prob = stats.binom_test(max_A, max_A+max_R, p=0.5, alternative='two-sided')
        # logging.info('... MS_prob - {}'.format(MS_prob))
        return round(SS_esti,3),round(SS_prob,3),round(MS_esti,3),round(MS_prob,3)
    else:
        return (None,None,None,None)
 
def getBaseline_pooled(fields,depth,hets):
    if(len(fields)>=4):
        base_thetas=[]
        Mreps=int(fields[1])
        pooled_A = 0
        pooled_R = 0
        pooled_min = 0
        for rep in range(Mreps):
            A = float(fields[2+rep*2])#+1
            R = float(fields[3+rep*2])#+1
            # logging.info('... {0} site: A-{1} R-{2}'.format(rep,A,R))
            pooled_A = pooled_A + A
            pooled_R = pooled_R + R
            pooled_min = pooled_min + min(A,R)
        sum_AR = pooled_A+pooled_R
        # logging.info('... {0} sites pooled: A-{1} R-{2}'.format(rep,pooled_A,pooled_R))
        # logging.info('... {0} sites pooled: minAR-{1} maxAR-{2}'.format(rep,pooled_min,sum_AR-pooled_min))
        # Naive Sum (NS): summing up the counts for A assuming phasing is correct
        NS_AAR = pooled_A/sum_AR
        if NS_AAR == 1 :            # if REF has 0 count
            NS_AAR = NS_AAR - 0.001
        NS_esti = NS_AAR/(1-NS_AAR)
        NS_prob=stats.binom_test(pooled_A, pooled_A+pooled_R, p=0.5, alternative='two-sided')
        # Pseudo phasing: lower counts assuming for ALT allele, high counts assuming for REF allele
        pseudo_AAR = pooled_min/sum_AR
        pseudo_esti = pseudo_AAR/(1-pseudo_AAR)
        pseudo_p = stats.binom_test(pooled_min, pooled_A+pooled_R, p=0.5, alternative='less')
        # logging.info('... NS_esti: {0}, NS_prob: {1}'.format(NS_esti,NS_prob))
        # logging.info('... pseudo_esti: {0}, pseudo_prob: {1}'.format(pseudo_esti,pseudo_p))
        return round(NS_esti,3), round(NS_prob,3),round(pseudo_esti,3),round(pseudo_p,3)
    else:
        return (None,None,None,None)

def run(prefix,inFile,out,picklename):
    # prefix="HG00096_chr21"
    # inFile="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr21/output/s-0.5_a-0.05_sinCov0_totCov1_W1000K1000/HG00096_chr21_hetSNP_intersected_filtered.TEMP.modelinput.tsv"
    # out="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr21/output/s-0.5_a-0.05_sinCov0_totCov1_W1000K1000"
    # picklename="HG00096_chr21_a-0.05_W1000K1000_s0t1.pickle"
    outfix=picklename
    out_path=out+"/output_pkl/binomial"
    # single site
    out_path_1 = out_path+"/SS_esti/"
    out_path_2 = out_path+"/SS_p/"
    # naive sum
    out_path_3 = out_path+"/NS_esti/"
    out_path_4 = out_path+"/NS_p/"
    # pseudo phasing
    out_path_5 = out_path+"/pseudo_esti/"
    out_path_6 = out_path+"/pseudo_p/"
    # major site
    out_path_7 = out_path+"/MS_esti/"
    out_path_8 = out_path+"/MS_p/"


    Path(out_path_1).mkdir(parents=True,exist_ok=True)
    Path(out_path_2).mkdir(parents=True,exist_ok=True)
    Path(out_path_3).mkdir(parents=True,exist_ok=True)
    Path(out_path_4).mkdir(parents=True,exist_ok=True)
    Path(out_path_5).mkdir(parents=True,exist_ok=True)
    Path(out_path_6).mkdir(parents=True,exist_ok=True)
    Path(out_path_7).mkdir(parents=True,exist_ok=True)
    Path(out_path_8).mkdir(parents=True,exist_ok=True)

    out1 = out_path_1 + str(outfix)
    out2 = out_path_2 + str(outfix)
    out3 = out_path_3 + str(outfix)
    out4 = out_path_4 + str(outfix)
    out5 = out_path_5 + str(outfix)
    out6 = out_path_6 + str(outfix)
    out7 = out_path_7 + str(outfix)
    out8 = out_path_8 + str(outfix)
    pseudo_esti_list=[]
    pseudo_p_list=[]
    SS_esti_list=[]
    SS_p_list=[]
    NS_esti_list=[]
    NS_p_list=[]
    MS_esti_list=[]
    MS_p_list=[]
    gene=[]
    if (os.path.isfile(out7)) and (os.path.isfile(out8)) and (os.path.isfile(out3)) and (os.path.isfile(out4)):
        logging.info('... binomial Already Processed {0}'.format(str(outfix)))
    else: 
        counter=0
        with open(inFile,"rt") as IN:
            #logging.info('... read}')
            #inFile="/Users/scarlett/Documents/Allen_lab/github/BEASTIE/BEASTIE_example/HG00096_chr21/output/s-0.5_a-0.05_sinCov0_totCov1_W1000K1000/HG00096_chr21_hetSNP_intersected_filtered.TEMP.modelinput.tsv"
            for line in IN:
                counter+=1
                # print(counter)
                # print(line)
                #logging.info('{0}'.format(line))
                fields=line.rstrip().split()   
                geneID=fields[0]
                gene.append(geneID)
                #logging.info('... gene {0}'.format(geneID))
                h = fields[1]
                d = int(fields[2])+int(fields[3])
                SS_esti,SS_prob,MS_esti,MS_prob = getBaseline(fields,int(d))
                NS_esti, NS_prob,pseudo_esti,pseudo_p = getBaseline_pooled(fields,int(d),int(h))
                #print("%s - %s - %s -%s"%(geneID,counter,SS_esti,NS_esti))
                #logging.info('... gene {0} - SS_esti {1} - SS_prob {2} - MS_esti {3} - MS_prob {4} -NS_esti {5}'.format(geneID,SS_esti,SS_prob,MS_esti,MS_prob,NS_esti))
                SS_esti_list.append(SS_esti)
                SS_p_list.append(SS_prob)
                NS_esti_list.append(NS_esti)
                NS_p_list.append(NS_prob)
                pseudo_esti_list.append(pseudo_esti)
                pseudo_p_list.append(pseudo_p)
                MS_esti_list.append(MS_esti)
                MS_p_list.append(MS_prob)

    binomial_df = pd.DataFrame(np.column_stack([gene, SS_esti_list, SS_p_list,NS_esti_list,NS_p_list,pseudo_esti_list,pseudo_p_list,MS_esti_list,MS_p_list]), 
                                columns=['geneID', 'FirstSite_esti', 'FirstSite_pval','NaiveSum_esti', 'NaiveSum_pval','Pseudo_esti', 'Pseudo_pval','MajorSite_esti', 'MajorSite_pval'])
    binomial_df.to_csv(out+"/"+prefix+"_ASE_binomial.tsv",sep="\t",header=True,index=False)
    logging.info('..... NS_p_list size {0}'.format(len(NS_p_list)))
    logging.info('..... SS_p_list size {0}'.format(len(SS_p_list)))
    logging.info('..... pseudo_p_list size {0}'.format(len(pseudo_p_list)))
    logging.info('..... MS_p_list size {0}'.format(len(MS_p_list)))
    pickle.dump(SS_esti_list,open(out1,'wb'))
    pickle.dump(SS_p_list,open(out2,'wb'))
    pickle.dump(NS_esti_list,open(out3,'wb'))
    pickle.dump(NS_p_list,open(out4,'wb'))
    pickle.dump(pseudo_esti_list,open(out5,'wb'))
    pickle.dump(pseudo_p_list,open(out6,'wb'))
    pickle.dump(MS_esti_list,open(out7,'wb'))
    pickle.dump(MS_p_list,open(out8,'wb'))
    return binomial_df
