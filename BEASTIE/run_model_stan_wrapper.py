#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
    generators,
    nested_scopes,
    with_statement,
)
from builtins import (
    bytes,
    dict,
    int,
    list,
    object,
    range,
    str,
    ascii,
    chr,
    hex,
    input,
    next,
    oct,
    open,
    pow,
    round,
    super,
    filter,
    map,
    zip,
)
import os
import sys
import math
import pickle
from StanParser import StanParser
import statistics
import pandas as pd
import math
import os.path
import TempFilename
import logging

def writeInitializationFile(filename):
    OUT=open(filename,"wt")
    print("theta <- 1",file=OUT)
    OUT.close()

def writeReadCounts(fields,start,numReps,varName,OUT):
    print(varName,"<- c(",file=OUT,end="")
    for rep in range(numReps):
        print(fields[start+rep*2],file=OUT,end="")
        if(rep+1<numReps): print(",",file=OUT,end="")
    print(")",file=OUT)

def writePi(fields,numReps,varName,OUT):
    print(varName,"<- c(",file=OUT,end="")
    start=numReps*2+3
    for rep in range(start,len(fields)):
        print(fields[rep],file=OUT,end="")
        if(rep+1<len(fields)): print(",",file=OUT,end="")
    print(")",file=OUT)

def writeInputsFile(fields,filename,sigma):
    Mreps=int(fields[1])
    OUT=open(filename,"wt")
    print("M <-",Mreps,file=OUT)
    writeReadCounts(fields,2,Mreps,"A",OUT) # alt
    writeReadCounts(fields,3,Mreps,"R",OUT) # ref
    print("sigma <-",sigma,file=OUT)
    OUT.close()

def writeInputsFile_i(fields,filename,sigma):
    Mreps=int(fields[1])
    OUT=open(filename,"wt")
    print("M <-",Mreps,file=OUT)
    writeReadCounts(fields,2,Mreps,"A",OUT) # alt
    writeReadCounts(fields,3,Mreps,"R",OUT) # ref
    print("sigma <-",sigma,file=OUT)
    writePi(fields,Mreps,"pi",OUT) # ref
    print("N_MISSING_PI <-",fields[Mreps*2+2],file=OUT)
    OUT.close()

def getBaseline(fields):
    if(len(fields)>=5):
        base_thetas=[]
        Mreps=int(fields[1])
        for rep in range(Mreps):
            A = float(fields[2+rep*2])
            R = float(fields[3+(rep)*2])
            base = (A+1)/(R+1)
            base_thetas.append(base)
        med_base_theta=statistics.median(base_thetas)
    return med_base_theta

def getFieldIndex(label,fields):
    numFields=len(fields)
    index=None
    for i in range(7,numFields):
        if(fields[i]==label): index=i
    return index

def getMaxProb_RMSE(thetas):
    p_less1 = len([i for i in thetas if i < 1])/len(thetas)
    p_more1 = 1-p_less1
    max_prob1 = max(p_less1,p_more1)
    # 2. transform thetas, and then calculate proportion
    thetas_log2= [math.log2(x) for x in thetas]
    p_less2 = len([i for i in thetas_log2 if i < 0])/len(thetas_log2)
    p_more2 = 1-p_less2
    max_prob2 = max(p_less2,p_more2)
    return max_prob1#, RMSE

def getMaxProb_lambda(thetas,Lambda):
    #print(thetas)
    # 1. no transformation 
    one_over_Lambda = float(1/float(Lambda))
    p_less1 = len([i for i in thetas if i < one_over_Lambda])/len(thetas)
    p_more1 = len([i for i in thetas if i > Lambda])/len(thetas)
    lambda_prob1 = max(p_less1,p_more1)
    # 2. transform thetas, and then calculate proportion
    #thetas_log2 = [math.log2(x) for x in thetas]
    #p_less2 = len([i for i in thetas_log2 if i < math.log2(one_over_Lambda)])/len(thetas)
    #p_more2 = len([i for i in thetas_log2 if i > math.log2(float(Lambda))])/len(thetas)
    #lambda_prob2 = max(p_less2,p_more2)
    # 3. sum tail
    lambda_sum1 = p_less1 + p_more1
    # 4. sum tail  transform thetas, and then calculate proportion
    #lambda_sum2 = p_less2 + p_more2
    return lambda_prob1,lambda_sum1

def runVariant(model,fields,tmp_output_file,stan_output_file,init_file,sigma,lambdas,KEEPER,WARMUP):
    if(len(fields)>=4):
        # if "iBEASTIE" in model:
        writeInputsFile_i(fields,tmp_output_file,sigma)
        # else:
        #     writeInputsFile(fields,tmp_output_file,sigma)
        writeInitializationFile(init_file)
        cmd = (
            "%s sample num_samples=%s num_warmup=%s data file=%s init=%s output file=%s refresh=0" 
            % (model,KEEPER,WARMUP,tmp_output_file,init_file,stan_output_file)
        )
        os.system(cmd)# Parse MCMC output
        parser=StanParser(stan_output_file)
        thetas=parser.getVariable("theta")
        med,_,_,_,_ = parser.getSummary("theta")
        max_prob = getMaxProb_RMSE(thetas)
        max_prob_lambda,sum_prob_lambda = getMaxProb_lambda(thetas,lambdas) 
        return (
            thetas,
            med,
            max_prob,
            max_prob_lambda,
            sum_prob_lambda,
        )
    else:
        logging.error('lines with no enough elements')

def runVariant_predictLambda_simulation(model,fields,tmp_output_file,stan_output_file,init_file,sigma,alphas,lambdas_file):
    if(len(fields)>=4):
        writeInputsFile(fields,tmp_output_file,sigma)
        #writeInputsFile(fields,tmp_output_file)
        writeInitializationFile(init_file)
        cmd = "%s sample data file=%s init=%s output file=%s" % (model,tmp_output_file,init_file,stan_output_file)
    #/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/ase sample data file=/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/tmp_output.txt init=/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/initialization_stan.txt output file=/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/output_theta.txt      
        #cmd = "%s sample data file=%s output file=%s" % (model,tmp_output_file,stan_output_file)
        os.system(cmd)# Parse MCMC output
        parser=StanParser(stan_output_file)
        thetas=parser.getVariable("theta")
        #med,_,_,_,_ = parser.getSummary("theta")
        #max_prob,_ = getMaxProb_RMSE(thetas,true_theta)
        prob_lambda1,sum_lambda1,prob_lambda2,sum_lambda2,prob_lambda3,sum_lambda3 = parse_lambda_validation(thetas,alphas,lambdas_file,1)
        prob_lambda11,sum_lambda11,prob_lambda22,sum_lambda22,prob_lambda33,sum_lambda33 = parse_lambda_validation(thetas,alphas,lambdas_file,2)
        return prob_lambda1,sum_lambda1,prob_lambda2,sum_lambda2,prob_lambda3,sum_lambda3,prob_lambda11,sum_lambda11,prob_lambda22,sum_lambda22,prob_lambda33,sum_lambda33

def parse_lambda_validation_simulation(thetas,alphas,lambdas_file,lm):
    with open(lambdas_file,"rt") as IN:
        for idx,line in enumerate(IN):
            #print(line)
            if idx == lm-1:
                fields=line.rstrip().split()
                print(fields)
                print("model %s - alpha %s - lambda: %s"%(lm,alphas[0],fields[0]))
                prob_lambda1,_,sum_lambda1,_ = getMaxProb_lambda(thetas,float(fields[0]))
                print("model %s - alpha %s - lambda: %s"%(lm,alphas[1],fields[1]))
                prob_lambda2,_,sum_lambda2,_ = getMaxProb_lambda(thetas,float(fields[1]))
                print("model %s - alpha %s - lambda: %s"%(lm,alphas[2],fields[2]))
                prob_lambda3,_,sum_lambda3,_ = getMaxProb_lambda(thetas,float(fields[2]))
                return prob_lambda1,sum_lambda1,prob_lambda2,sum_lambda2,prob_lambda3,sum_lambda3

def getMedian(thetas):
    # Precondition: thetas is already sorted
    thetas.sort()
    n=len(thetas)
    mid=int(n/2)
    if(n%2==0): return (thetas[mid-1]+thetas[mid])/2.0
    return thetas[mid]

def getCredibleInterval(thetas,alpha):
    halfAlpha=alpha/2.0
    n=len(thetas)
    leftIndex=int(halfAlpha*n)
    rightIndex=n-leftIndex
    left=thetas[leftIndex+1]
    right=thetas[rightIndex-1]
    return (left,right)

def summarize(thetas,alpha):
    thetas.sort()
    n=len(thetas)
    median=getMedian(thetas)
    (CI_left,CI_right)=getCredibleInterval(thetas,alpha)
    print(median,CI_left,CI_right,sep="\t")

def parse_stan_output(out1,out2,out3,out_max_lambda,out_sum_lambda,models,input_file,tmp_output_file,stan_output_file,init_file,sigma,out_path,outfix,alphas,lambdas_file,KEEPER=1000,WARMUP=1000):
    prob_tail_lambda = []
    prob_sum_lambda = []
    lambdas=pd.read_csv(lambdas_file, delimiter='\t', names = ['gene_ID','median_altratio','num_hets','totalRef','totalAlt','total_reads','predicted_lambda'])
    model_theta = []       # 150
    model_theta_med = []   # 150
    model_med_prob = []    # 150
    with open(input_file,"rt") as IN:
        i=0
        for line in IN:
            i+=1
            fields=line.rstrip().split()
            ID=fields[0]
            lambdas_choice=lambdas.loc[lambdas['gene_ID'] == ID].iloc[0,6]
            thetas,med,prob,tail_prob_lambda,sum_prob_lambda =runVariant(models,fields,tmp_output_file,stan_output_file,init_file,sigma,lambdas_choice,KEEPER,WARMUP)
            prob_tail_lambda.append(tail_prob_lambda)
            prob_sum_lambda.append(sum_prob_lambda)
            model_theta.extend(thetas)
            model_theta_med.append(med)
            model_med_prob.append(prob)        
    if not os.path.exists(out_path+"model_theta/"):os.makedirs(out_path+"model_theta/")
    if not os.path.exists(out_path+"model_theta_med/"):os.makedirs(out_path+"model_theta_med/")
    if not os.path.exists(out_path+"model_prob/"):os.makedirs(out_path+"model_prob/")

    pickle.dump(model_theta,open(out1,'wb'))
    pickle.dump(model_theta_med,open(out2,'wb'))
    pickle.dump(model_med_prob,open(out3,'wb'))
    pickle.dump(prob_sum_lambda,open(out_sum_lambda,'wb'))
    pickle.dump(prob_tail_lambda,open(out_max_lambda,'wb'))


def run(inFile,sigma,alpha,models,out_path,lambdas_file,WARMUP=1000,KEEPER=1000):
    if "txt" in inFile:
        outfix=os.path.split(inFile)[1].rsplit(".txt")[0]
    if "tsv" in inFile:
        outfix=os.path.split(inFile)[1].rsplit(".tsv")[0]
    tmpFile= "tmp_output.txt"
    initFile = "initialization_stan.txt"
    outFile= "stan_output.txt"
    out_path=out_path+"/output_pkl/"
    if not os.path.exists(out_path):os.makedirs(out_path)
    tmp_output_file=out_path+tmpFile
    init_file=out_path+initFile
    stan_output_file=out_path+outFile
    ###########################################################################################
    outname1=str(outfix)+"_s-"+str(sigma)+".pickle"
    out1 = out_path+"model_theta/"+outname1
    out2 = out_path+"model_theta_med/"+outname1
    out3 = out_path+"model_prob/"+outname1

    outname2=str(outfix)+"_s-"+str(sigma)+"_a-"+str(alpha)+".pickle"
    out_max = out_path+"model_prob_tail_lambda_predicted/"
    out_sum = out_path+"model_prob_sum_lambda_predicted/"
    if not os.path.exists(out_max):os.makedirs(out_max)
    if not os.path.exists(out_sum):os.makedirs(out_sum)
    out_max_lambda = out_max+outname2
    out_sum_lambda = out_sum+outname2

    if (os.path.isfile(out1)) and (os.path.isfile(out2)) and (os.path.isfile(out3)) and (os.path.isfile(out_max_lambda)) and (os.path.isfile(out_sum_lambda)):
        logging.info('.... Already Processed {0}'.format(out1))
        logging.info('.... Already Processed {0}'.format(out2))
        logging.info('.... Already Processed {0}'.format(out3))
        logging.info('.... Already Processed {0}'.format(out_max_lambda))
        logging.info('.... Already Processed {0}'.format(out_sum_lambda))
    else:
        parse_stan_output(out1,out2,out3,out_max_lambda,out_sum_lambda,models,inFile,tmp_output_file,stan_output_file,init_file,sigma,out_path,outfix,alpha,lambdas_file,WARMUP,KEEPER)
    return out_path,outname1,outname2