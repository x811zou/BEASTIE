#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
import logging
import math
import os
import os.path
import pickle
import statistics

import pandas as pd

from StanParser import StanParser


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
    # 1. no transformation
    one_over_Lambda = float(1/float(Lambda))
    # changes maded in 07/22
    min_l = min(Lambda,one_over_Lambda)
    max_l = max(Lambda,one_over_Lambda)
    #
    p_less1 = len([i for i in thetas if i < min_l])/len(thetas)
    p_more1 = len([i for i in thetas if i > max_l])/len(thetas)
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
    return round(lambda_prob1,3),round(lambda_sum1,3)

def runModel(model,fields,tmp_output_file,stan_output_file,init_file,sigma,WARMUP,KEEPER):
    if(len(fields)>=4):
        writeInputsFile_i(fields,tmp_output_file,sigma)
        writeInitializationFile(init_file)
        cmd = (
            "%s sample num_samples=%s num_warmup=%s data file=%s init=%s output file=%s refresh=0"
            % (model,KEEPER,WARMUP,tmp_output_file,init_file,stan_output_file)
        )
        os.system(cmd)# Parse MCMC output
        parser=StanParser(stan_output_file)
        thetas=parser.getVariable("theta")
        return thetas
    else:
        logging.error('lines with no enough elements')

def parse_lambda_validation_simulation(thetas,alphas,lambdas_file,lm):
    with open(lambdas_file,"rt") as IN:
        for idx,line in enumerate(IN):
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
    return left,right

def summarize(thetas,alpha):
    thetas.sort()
    median=getMedian(thetas)
    CI_left,CI_right=getCredibleInterval(thetas,alpha)
    return median,CI_left,CI_right

def parse_stan_output(input_file,out1,KEEPER,lambdas_file):
    thetas = pickle.load(
        open(
            out1,
            "rb",
        )
    )
    lambdas=pd.read_csv(lambdas_file, delimiter='\t', names = ['gene_ID','median_altratio','num_hets','totalRef','totalAlt','total_reads','predicted_lambda'])
    prob_sum_lambda = []
    model_theta_med = []   # 150
    CI_left=[]
    CI_right=[]
    geneID=[]
    with open(input_file,"rt") as IN:
        i=0
        for line in IN:
            fields=line.rstrip().split()
            ID=fields[0]
            geneID.append(ID)         # read the ith geneID
            j=i+int(KEEPER)-1
            gene_thetas=thetas[i:j]
            lambdas_choice=lambdas.loc[lambdas['gene_ID'] == ID].iloc[0,6]
            median,left_CI,right_CI = summarize(gene_thetas,0.05)
            max_prob = getMaxProb_RMSE(gene_thetas)
            max_prob_lambda,sum_prob_lambda = getMaxProb_lambda(gene_thetas,lambdas_choice)
            i=i+int(KEEPER)
            prob_sum_lambda.append(sum_prob_lambda)
            CI_left.append(round(left_CI,3))
            CI_right.append(round(right_CI,3))
            model_theta_med.append(round(median,3))
    logging.debug('size of thetas : {0}, size of output list :{1}'.format(len(thetas),len(prob_sum_lambda)))
    df={'gene_ID':geneID,'posterior_median':model_theta_med,'CI_left':CI_left,'CI_right':CI_right,'posterior_mass_support_ALT':prob_sum_lambda}
    df=pd.DataFrame(df)
    return df


def save_raw_theta(out0,models,input_file,tmp_output_file,stan_output_file,init_file,sigma,WARMUP,KEEPER):
    model_theta = []       # 150
    with open(input_file,"rt") as IN:
        i=0
        for line in IN:
            i+=1
            fields=line.rstrip().split()
            thetas =runModel(models,fields,tmp_output_file,stan_output_file,init_file,sigma,WARMUP,KEEPER)
            model_theta.extend(thetas)
    logging.debug('.... Number of genes : {0}, and length of thetas: {1}'.format(i,len(model_theta)))
    pickle.dump(model_theta,open(out0,'wb'))


def run(prefix,inFile,sigma,alpha,models,out,lambdas_file,WARMUP,KEEPER,either_cov,total_cov):
    logging.debug("Number of WARMUP samples is {0}, Number of posterior estimates is {1}".format(WARMUP,KEEPER))
    if "txt" in inFile:
        outfix=os.path.split(inFile)[1].rsplit("_hetSNP_intersected_filtered.TEMP.txt")[0]
    if "tsv" in inFile:
        outfix=os.path.split(inFile)[1].rsplit("_hetSNP_intersected_filtered.TEMP.tsv")[0]
    tmpFile= "tmp_output.txt"
    initFile = "initialization_stan.txt"
    outFile= "stan_output.txt"
    out=out+"/output_pkl"
    out_path=out+"/beastie/"
    if not os.path.exists(out):os.makedirs(out)
    if not os.path.exists(out_path):os.makedirs(out_path)
    tmp_output_file=out_path+tmpFile
    init_file=out_path+initFile
    stan_output_file=out_path+outFile
    ###########################################################################################
    outname1=str(prefix)+"_a-"+str(alpha)+"_W"+str(WARMUP)+"K"+str(KEEPER)+"_s"+str(either_cov)+"t"+str(total_cov)+"_s-"+str(sigma)+".pickle"
    out_path = out_path+"theta/"
    if not os.path.exists(out_path):os.makedirs(out_path)
    out1 = out_path+outname1
    if (os.path.isfile(out1)):
        logging.info('.... Already fiinshed running model and saved raw theta at : {0}'.format(out1))
    else:
        save_raw_theta(out1,models,inFile,tmp_output_file,stan_output_file,init_file,sigma,WARMUP,KEEPER)
        logging.info('.... Finish running model and save raw theta at : {0}'.format(out1))
    df=parse_stan_output(inFile,out1,KEEPER,lambdas_file)
    outname2=str(prefix)+"_a-"+str(alpha)+"_W"+str(WARMUP)+"K"+str(KEEPER)+"_s"+str(either_cov)+"t"+str(total_cov)+".pickle"
    return df,outname2
