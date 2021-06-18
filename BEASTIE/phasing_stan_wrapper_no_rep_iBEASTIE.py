#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
import os
import sys
import math
#from __future__ import print_function
import pickle
from StanParser import StanParser
#import numpy as np
import statistics
import pandas as pd
import math
import os.path
import TempFilename
# class RunStan:
#    def __init__(self, simulated_data_file):
#        self.simulated_data_file = simulated_data_file
#    def __enter__(self):
#    def __exit__(self ,type, value, traceback):


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
    #filename='output.txt'
    OUT=open(filename,"wt")
    print("M <-",Mreps,file=OUT)
    writeReadCounts(fields,2,Mreps,"A",OUT) # alt
    writeReadCounts(fields,3,Mreps,"R",OUT) # ref
    print("sigma <-",sigma,file=OUT)
    writePi(fields,Mreps,"pi",OUT) # ref
    #print()
    OUT.close()

def getBaseline(fields):
    if(len(fields)>=5):
        base_thetas=[]
        Mreps=int(fields[1])
        for rep in range(Mreps):
            A = float(fields[2+rep*2])
            R = float(fields[3+(rep)*2])
            base = (A+1)/(R+1)
            #abs_base = abs(base-1)
            base_thetas.append(base)
        med_base_theta=statistics.median(base_thetas)
    #true_theta = fields[-1]
    return med_base_theta

    #with tempfile.NamedTemporaryFile() as output_file:
    #    with tempfile.NamedTemporaryFile() as init_file:
    #        with tempfile.NamedTemporaryFile() as input_file:# Write inputs file for STAN

def getFieldIndex(label,fields):
    numFields=len(fields)
    index=None
    for i in range(7,numFields):
        if(fields[i]==label): index=i
    return index

def getMaxProb_RMSE(thetas):
    #print(thetas)
    # 1. no transformation 
    p_less1 = len([i for i in thetas if i < 1])/len(thetas)
    p_more1 = 1-p_less1
    max_prob1 = max(p_less1,p_more1)
    # 2. transform thetas, and then calculate proportion
    thetas_log2= [math.log2(x) for x in thetas]
    p_less2 = len([i for i in thetas_log2 if i < 0])/len(thetas_log2)
    p_more2 = 1-p_less2
    max_prob2 = max(p_less2,p_more2)
    #diff = [(x - 1)**2 for x in thetas]
    #rmse = np.sqrt(np.mean(diff))
    ############# theta
    #thetas_diff = [(x-true_theta) for x in thetas]
    #thetas_sq = list(map(lambda n:n**2,thetas_diff))
    #RMSE=math.sqrt(sum(thetas_sq)/len(thetas_sq))
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

def runVariant(model,fields,tmp_output_file,stan_output_file,init_file,sigma,predicted_lambda_modeldata,lambdas):
    if(len(fields)>=4):
        writeInputsFile(fields,tmp_output_file,sigma)
        #writeInputsFile(fields,tmp_output_file)
        writeInitializationFile(init_file)
        cmd = "%s sample num_samples=%s num_warmup=%s data file=%s init=%s output file=%s refresh=0" % (model,KEEPER,WARMUP,tmp_output_file,init_file,stan_output_file)
    #/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/ase sample data file=/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/tmp_output.txt init=/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/initialization_stan.txt output file=/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/output_theta.txt      
        #cmd = "%s sample data file=%s output file=%s" % (model,tmp_output_file,stan_output_file)
        os.system(cmd)# Parse MCMC output
        parser=StanParser(stan_output_file)
        thetas=parser.getVariable("theta")
        med,_,_,_,_ = parser.getSummary("theta")
        max_prob = getMaxProb_RMSE(thetas)
        max_prob_lambda12,sum_lambda12 = getMaxProb_lambda(thetas,1.2)
        max_prob_lambda13,sum_lambda13 = getMaxProb_lambda(thetas,1.3)
        max_prob_lambda14,sum_lambda14 = getMaxProb_lambda(thetas,1.4)
        max_prob_lambda15,sum_lambda15 = getMaxProb_lambda(thetas,1.5)
        max_prob_lambda16,sum_lambda16 = getMaxProb_lambda(thetas,1.6)
        max_prob_lambda18,sum_lambda18 = getMaxProb_lambda(thetas,1.8)
        if str(predicted_lambda_modeldata) == "True":
            max_prob_lambda,sum_prob_lambda = getMaxProb_lambda(thetas,lambdas)
        else:
            max_prob_lambda = None
            sum_prob_lambda = None    
        return thetas,med,max_prob,max_prob_lambda12,sum_lambda12,max_prob_lambda13,sum_lambda13,max_prob_lambda14,sum_lambda14,max_prob_lambda15,sum_lambda15,max_prob_lambda16,sum_lambda16,max_prob_lambda18,sum_lambda18,max_prob_lambda,sum_prob_lambda

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
            print(line)
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

def parse_stan_output(models,input_file,tmp_output_file,stan_output_file,init_file,sigma,out_path,outfix,alphas,alpha_list,lambdas_file=None,predicted_lambda_modeldata=False,default_lambda=False,selfdefined_lambda=False):
    ############################################################################################
    ##### 1. default
    #if default_lambda is True:
    if str(default_lambda) == "True":
        print(">>>>>>>>>>>>>>>>>>>>>>>>>> default_lambda")
        if str(predicted_lambda_modeldata) == "True":
            print(">>>>>>>>>>>>>>>>>>>>>>>>>> predicted_lambda_modeldata")
            prob_tail_lambda = []
            prob_sum_lambda = []
            lambdas=pd.read_csv(lambdas_file, delimiter='\t', names = ['gene_ID','median_altratio','num_hets','totalRef','totalAlt','total_reads','predicted_lambda'])
            #geneID	median altRatio	number of hets	refCount	altCount	totalCount
        model_theta = []       # 150
        model_theta_med = []   # 150
        model_med_prob = []    # 150
        model_prob_tail_lambda2 = []
        model_prob_tail_lambda3 = []
        model_prob_tail_lambda4 = []
        model_prob_tail_lambda5 = []
        model_prob_tail_lambda6 = []
        model_prob_tail_lambda8 = [] 
        model_prob_sum_lambda2 = []
        model_prob_sum_lambda3 = []
        model_prob_sum_lambda4 = []
        model_prob_sum_lambda5 = []
        model_prob_sum_lambda6 = []
        model_prob_sum_lambda8 = []
        with open(input_file,"rt") as IN:
            i=0
            for line in IN:
                print(line)
                i+=1
                fields=line.rstrip().split()
                ID=fields[0]
                nhets=fields[1]
                #true_theta=float(fields[2+int(nhets)*2])
                if str(predicted_lambda_modeldata) == "True":
                    lambdas_choice=lambdas.iloc[i-1,6]
                else:
                    lambdas_choice = None
                thetas,med,prob,lambda2,sum2,lambda3,sum3,lambda4,sum4,lambda5,sum5,lambda6,sum6,lambda8,sum8,tail_prob_lambda,sum_prob_lambda =runVariant(models,fields,tmp_output_file,stan_output_file,init_file,sigma,predicted_lambda_modeldata,lambdas_choice)
                if str(predicted_lambda_modeldata) == "True":
                    prob_tail_lambda.append(tail_prob_lambda)
                    prob_sum_lambda.append(sum_prob_lambda)
                model_theta.extend(thetas)
                model_theta_med.append(med)
                model_med_prob.append(prob)
                model_prob_tail_lambda2.append(lambda2)
                model_prob_tail_lambda3.append(lambda3)
                model_prob_tail_lambda4.append(lambda4)
                model_prob_tail_lambda5.append(lambda5)
                model_prob_tail_lambda6.append(lambda6)
                model_prob_tail_lambda8.append(lambda8)
                model_prob_sum_lambda2.append(sum2)
                model_prob_sum_lambda3.append(sum3)
                model_prob_sum_lambda4.append(sum4)
                model_prob_sum_lambda5.append(sum5)
                model_prob_sum_lambda6.append(sum6)
                model_prob_sum_lambda8.append(sum8)
        if str(predicted_lambda_modeldata) == "True":
            out_max = out_path+"model_prob_tail_lambda_predicted/"
            out_sum = out_path+"model_prob_sum_lambda_predicted/"
            if not os.path.exists(out_max):
                os.makedirs(out_max)
            if not os.path.exists(out_sum):
                os.makedirs(out_sum)
            out_max_lambda = out_max+str(outfix)+"_s-"+str(sigma)+"_a-"+str(alphas)+".pickle"
            out_sum_lambda = out_sum+str(outfix)+"_s-"+str(sigma)+"_a-"+str(alphas)+".pickle"
            pickle.dump(prob_sum_lambda,open(out_sum_lambda,'wb'))
            pickle.dump(prob_tail_lambda,open(out_max_lambda,'wb'))
            print("done with saving for predicted lambda")            
        if not os.path.exists(out_path+"model_theta/"):os.makedirs(out_path+"model_theta/")
        if not os.path.exists(out_path+"model_theta_med/"):os.makedirs(out_path+"model_theta_med/")
        if not os.path.exists(out_path+"model_prob/"):os.makedirs(out_path+"model_prob/")
        if not os.path.exists(out_path+"model_prob_tail_lambda1.2/"):os.makedirs(out_path+"model_prob_tail_lambda1.2/")
        if not os.path.exists(out_path+"model_prob_tail_lambda1.3/"):os.makedirs(out_path+"model_prob_tail_lambda1.3/")
        if not os.path.exists(out_path+"model_prob_tail_lambda1.4/"):os.makedirs(out_path+"model_prob_tail_lambda1.4/")
        if not os.path.exists(out_path+"model_prob_tail_lambda1.5/"):os.makedirs(out_path+"model_prob_tail_lambda1.5/")
        if not os.path.exists(out_path+"model_prob_tail_lambda1.6/"):os.makedirs(out_path+"model_prob_tail_lambda1.6/")
        if not os.path.exists(out_path+"model_prob_tail_lambda1.8/"):os.makedirs(out_path+"model_prob_tail_lambda1.8/")
        if not os.path.exists(out_path+"model_prob_sum_lambda1.2/"):os.makedirs(out_path+"model_prob_sum_lambda1.2/")
        if not os.path.exists(out_path+"model_prob_sum_lambda1.3/"):os.makedirs(out_path+"model_prob_sum_lambda1.3/")
        if not os.path.exists(out_path+"model_prob_sum_lambda1.4/"):os.makedirs(out_path+"model_prob_sum_lambda1.4/")
        if not os.path.exists(out_path+"model_prob_sum_lambda1.5/"):os.makedirs(out_path+"model_prob_sum_lambda1.5/")
        if not os.path.exists(out_path+"model_prob_sum_lambda1.6/"):os.makedirs(out_path+"model_prob_sum_lambda1.6/")
        if not os.path.exists(out_path+"model_prob_sum_lambda1.8/"):os.makedirs(out_path+"model_prob_sum_lambda1.8/")
        out1 = out_path+"model_theta/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out2 = out_path+"model_theta_med/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out3 = out_path+"model_prob/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out4 = out_path+"model_prob_tail_lambda1.2/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out5 = out_path+"model_prob_tail_lambda1.3/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out6 = out_path+"model_prob_tail_lambda1.4/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out7 = out_path+"model_prob_tail_lambda1.5/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out8 = out_path+"model_prob_tail_lambda1.6/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out9 = out_path+"model_prob_tail_lambda1.8/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out10 = out_path+"model_prob_sum_lambda1.2/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out11 = out_path+"model_prob_sum_lambda1.3/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out12 = out_path+"model_prob_sum_lambda1.4/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out13 = out_path+"model_prob_sum_lambda1.5/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out14 = out_path+"model_prob_sum_lambda1.6/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        out15 = out_path+"model_prob_sum_lambda1.8/"+str(outfix)+"_s-"+str(sigma)+".pickle"
        # if (os.path.isfile(out26)) and (os.path.isfile(out35)):# and (os.path.isfile(out3)):
        #     print ("Already Processed"+str(outfix)+"_s-"+str(sigma))
        #     os._exit(0)
        pickle.dump(model_theta,open(out1,'wb'))
        pickle.dump(model_theta_med,open(out2,'wb'))
        pickle.dump(model_med_prob,open(out3,'wb'))
        pickle.dump(model_prob_tail_lambda2,open(out4,'wb'))
        pickle.dump(model_prob_tail_lambda3,open(out5,'wb'))
        pickle.dump(model_prob_tail_lambda4,open(out6,'wb'))
        pickle.dump(model_prob_tail_lambda5,open(out7,'wb'))
        pickle.dump(model_prob_tail_lambda6,open(out8,'wb'))
        pickle.dump(model_prob_tail_lambda8,open(out9,'wb'))
        pickle.dump(model_prob_sum_lambda2,open(out10,'wb'))
        pickle.dump(model_prob_sum_lambda3,open(out11,'wb'))
        pickle.dump(model_prob_sum_lambda4,open(out12,'wb'))
        pickle.dump(model_prob_sum_lambda5,open(out13,'wb'))
        pickle.dump(model_prob_sum_lambda6,open(out14,'wb'))
        pickle.dump(model_prob_sum_lambda8,open(out15,'wb'))

    ############################################################################################
    ##### 2. predict lambda for simulation 
    #elif selfdefined_lambda is True:
    if str(selfdefined_lambda) == "True":
        print(">>>>>>>>>>>>>>>>>>>>>>>>>> selfdefined_lambda")
        prob_sum_lambda1 = []
        prob_max_lambda1 = []
        prob_sum_lambda2 = []
        prob_max_lambda2 = []
        prob_sum_lambda3 = []
        prob_max_lambda3 = []
        prob_sum_lambda11 = []
        prob_max_lambda11 = []
        prob_sum_lambda22 = []
        prob_max_lambda22 = []
        prob_sum_lambda33 = []
        prob_max_lambda33 = []
        with open(input_file,"rt") as IN:
            i=0
            for line in IN:
                print(line)
                i+=1
                fields=line.rstrip().split()
                ID=fields[0]
                nhets=fields[1]
                #true_theta=float(fields[2+int(nhets)*2])
                prob_lambda1,sum_lambda1,prob_lambda2,sum_lambda2,prob_lambda3,sum_lambda3,prob_lambda11,sum_lambda11,prob_lambda22,sum_lambda22,prob_lambda33,sum_lambda33 =runVariant_predictLambda_simulation(models,fields,tmp_output_file,stan_output_file,init_file,sigma,alpha_list,lambdas_file)
                prob_max_lambda1.append(prob_lambda1)
                prob_sum_lambda1.append(sum_lambda1)
                prob_max_lambda2.append(prob_lambda2)
                prob_sum_lambda2.append(sum_lambda2)
                prob_max_lambda3.append(prob_lambda3)
                prob_sum_lambda3.append(sum_lambda3)
                prob_max_lambda11.append(prob_lambda11)
                prob_sum_lambda11.append(sum_lambda11)
                prob_max_lambda22.append(prob_lambda22)
                prob_sum_lambda22.append(sum_lambda22)
                prob_max_lambda33.append(prob_lambda33)
                prob_sum_lambda33.append(sum_lambda33)
        out_max = out_path+"model_prob_tail_lambda_predicted/"
        out_sum = out_path+"model_prob_sum_lambda_predicted/"
        if not os.path.exists(out_max):
            os.makedirs(out_max)
        if not os.path.exists(out_sum):
            os.makedirs(out_sum)
        out_max_lambda1 = out_max+str(outfix)+"_s-"+str(sigma)+"_m-1_alpha-"+str(alpha_list[0])+".pickle"
        out_sum_lambda1 = out_sum+str(outfix)+"_s-"+str(sigma)+"_m-1_alpha-"+str(alpha_list[0])+".pickle"
        out_max_lambda2 = out_max+str(outfix)+"_s-"+str(sigma)+"_m-1_alpha-"+str(alpha_list[1])+".pickle"
        out_sum_lambda2 = out_sum+str(outfix)+"_s-"+str(sigma)+"_m-1_alpha-"+str(alpha_list[1])+".pickle"
        out_max_lambda3 = out_max+str(outfix)+"_s-"+str(sigma)+"_m-1_alpha-"+str(alpha_list[2])+".pickle"
        out_sum_lambda3 = out_sum+str(outfix)+"_s-"+str(sigma)+"_m-1_alpha-"+str(alpha_list[2])+".pickle"
        out_max_lambda11 = out_max+str(outfix)+"_s-"+str(sigma)+"_m-2_alpha-"+str(alpha_list[0])+".pickle"
        out_sum_lambda11 = out_sum+str(outfix)+"_s-"+str(sigma)+"_m-2_alpha-"+str(alpha_list[0])+".pickle"
        out_max_lambda22 = out_max+str(outfix)+"_s-"+str(sigma)+"_m-2_alpha-"+str(alpha_list[1])+".pickle"
        out_sum_lambda22 = out_sum+str(outfix)+"_s-"+str(sigma)+"_m-2_alpha-"+str(alpha_list[1])+".pickle"
        out_max_lambda33 = out_max+str(outfix)+"_s-"+str(sigma)+"_m-2_alpha-"+str(alpha_list[2])+".pickle"
        out_sum_lambda33 = out_sum+str(outfix)+"_s-"+str(sigma)+"_m-2_alpha-"+str(alpha_list[2])+".pickle"
        pickle.dump(prob_max_lambda1,open(out_max_lambda1,'wb'))
        pickle.dump(prob_sum_lambda1,open(out_sum_lambda1,'wb'))
        pickle.dump(prob_max_lambda2,open(out_max_lambda2,'wb'))
        pickle.dump(prob_sum_lambda2,open(out_sum_lambda2,'wb'))
        pickle.dump(prob_max_lambda3,open(out_max_lambda3,'wb'))
        pickle.dump(prob_sum_lambda3,open(out_sum_lambda3,'wb'))
        pickle.dump(prob_max_lambda11,open(out_max_lambda11,'wb'))
        pickle.dump(prob_sum_lambda11,open(out_sum_lambda11,'wb'))
        pickle.dump(prob_max_lambda22,open(out_max_lambda22,'wb'))
        pickle.dump(prob_sum_lambda22,open(out_sum_lambda22,'wb'))
        pickle.dump(prob_max_lambda33,open(out_max_lambda33,'wb'))
        pickle.dump(prob_sum_lambda33,open(out_sum_lambda33,'wb'))

def main():
    ####################### read in parameters from input
    inFile=sys.argv[1]
    sigma=sys.argv[2]
    model=sys.argv[4]
    in_path=sys.argv[5]
    out_path=sys.argv[6]
    default_lambda=sys.argv[7]
    #print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print("deafult lambda is %s"%(default_lambda))
    selfdefined_lambda=sys.argv[8]
    print("selfdefined_lambda is %s"%(selfdefined_lambda))
    predicted_lambda_modeldata=sys.argv[9]
    print("predicted_lambda_modeldata is %s"%(predicted_lambda_modeldata))
    lambdas_file=sys.argv[10]
    alphas=sys.argv[11]
    ####################### self-defined parameters
    outfix=inFile.rsplit(".txt")[0]
    selfdefined_alpha=[0.05,0.01,0.001]
    models =model
    tmpFile= "tmp_output."+sys.argv[3]+".txt"
    initFile = "initialization_stan."+sys.argv[3]+".txt"
    outFile= "stan_output."+sys.argv[3]+".txt"
    input_file=in_path+inFile
    tmp_output_file=out_path+tmpFile
    init_file=out_path+initFile
    stan_output_file=out_path+outFile
    ###########################################################################################

    parse_stan_output(models,input_file,tmp_output_file,stan_output_file,init_file,sigma,out_path,outfix,alphas,selfdefined_alpha,lambdas_file=lambdas_file,predicted_lambda_modeldata=predicted_lambda_modeldata,default_lambda=default_lambda,selfdefined_lambda=selfdefined_lambda)

if __name__ == '__main__':
    WARMUP=1000
    KEEPER=1000
    DEBUG=False
    STDERR=TempFilename.generate(".stderr")
    INPUT_FILE=TempFilename.generate(".staninputs")
    INIT_FILE=TempFilename.generate(".staninit")
    OUTPUT_TEMP=TempFilename.generate(".stanoutputs")
    main()
    print("finished running data file")
