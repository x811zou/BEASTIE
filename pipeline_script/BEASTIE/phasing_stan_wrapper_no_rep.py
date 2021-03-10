#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
import os
import tempfile
import sys
import math
import ProgramName
import TempFilename
import getopt
#from __future__ import print_function
from Pipe import Pipe
import pickle
from StanParser import StanParser
#import numpy as np
import statistics
import math
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
        print(fields[numReps*2+3],file=OUT,end="")
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
    #writePi(fields,Mreps,"pi",OUT) # ref
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
    true_theta = fields[-1]
    return med_base_theta, true_theta

    #with tempfile.NamedTemporaryFile() as output_file:
    #    with tempfile.NamedTemporaryFile() as init_file:
    #        with tempfile.NamedTemporaryFile() as input_file:# Write inputs file for STAN

def getFieldIndex(label,fields):
    numFields=len(fields)
    index=None
    for i in range(7,numFields):
        if(fields[i]==label): index=i
    return index


def getMaxProb(thetas):
    # 1. no transformation 
    p_less1 = len([i for i in thetas if i < 1])/len(thetas)
    p_more1 = 1-p_less1
    max_prob1 = max(p_less1,p_more1)
    # 2. transform thetas, and then calculate proportion
    thetas[:] = [math.log2(x) for x in thetas]
    p_less2 = len([i for i in thetas if i < 0])/len(thetas)
    p_more2 = 1-p_less2
    max_prob2 = max(p_less2,p_more2)
    #diff = [(x - 1)**2 for x in thetas]
    #rmse = np.sqrt(np.mean(diff))
    return max_prob1,max_prob2

def runVariant(model,fields,input_file,tmp_output_file,stan_output_file,init_file,sigma):
    if(len(fields)>=5):
        writeInputsFile(fields,tmp_output_file,sigma)
        #writeInputsFile(fields,tmp_output_file)
        writeInitializationFile(init_file)
        cmd = "%s sample data file=%s init=%s output file=%s" % (model,tmp_output_file,init_file,stan_output_file)
    #/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/ase sample data file=/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/tmp_output.txt init=/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/initialization_stan.txt output file=/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/output_theta.txt      
        #cmd = "%s sample data file=%s output file=%s" % (model,tmp_output_file,stan_output_file)
        #print (cmd)
        os.system(cmd)# Parse MCMC output
        #output=Pipe.run(cmd)
        parser=StanParser(stan_output_file)
        thetas=parser.getVariable("theta")
        med,_,_,_,_ = parser.getSummary("theta")
        max_prob1,max_prob2 = getMaxProb(thetas)
        # with open(stan_output_file,"rt") as IN:
        #     for line in IN:
        #         if(len(line)==0 or line[0]=="#"): continue
        #         fields=line.rstrip().split(",")
        #         numFields=len(fields)
        #         if(fields[0]=="lp__"):
        #             #printFields(fields,OUT)
        #             thetaIndex=getFieldIndex("theta",fields)
        #             continue
        #         theta=float(fields[thetaIndex])
        #         thetas.append(theta)
        #         #print("append")
        #         #thetas.sort()
        
        return med, max_prob1,max_prob2
    else:
        return (None,None,None)

def getMedian(thetas):
    # Precondition: thetas is already sorted
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

def summarize(thetas):
    thetas.sort()
    n=len(thetas)
    median=getMedian(thetas)
    (CI_left,CI_right)=getCredibleInterval(thetas,ALPHA)
    print(median,CI_left,CI_right,sep="\t")


#inFile = "simulate_data1.txt"
# outFile= "stan_output.txt"
# model = "ase"
# #folder=sys.argv[2]
# folder='theta_1'
# output_path="/data/reddylab/scarlett/1000G/software/cmdstan/examples/ase/"+folder+"/"
# THETA = output_path+"theta_stan.txt"

model=sys.argv[4]
in_path=sys.argv[5]+"/"
out_path="/data/allenlab/scarlett/software/cmdstan/examples/"+str(model)+"/"+sys.argv[6]+"/output_pkl/" #+sys.argv[3]+"/"
#out_path="./output_pkl/"
inFile=sys.argv[1]
#inFile="g-50_h-5_d-50_t-1.txt"
#sigma=0.1
sigma=sys.argv[2]
outfix=inFile.rsplit(".txt")[0]

import os.path

if not os.path.exists(out_path+"model_med/"):
    os.makedirs(out_path+"model_med/")
if not os.path.exists(out_path+"model_prob1/"):
    os.makedirs(out_path+"model_prob1/")
if not os.path.exists(out_path+"model_prob2/"):
    os.makedirs(out_path+"model_prob2/")

#out1 = out_path+"model_theta/"+str(outfix)+"_s-"+str(sigma)+".pickle"
out2 = out_path+"model_med/"+str(outfix)+"_s-"+str(sigma)+".pickle"
out3 = out_path+"model_prob1/"+str(outfix)+"_s-"+str(sigma)+".pickle"
out4 = out_path+"model_prob2/"+str(outfix)+"_s-"+str(sigma)+".pickle"
#out4 = out_path+"baseline_theta/"+str(outfix)+"_s-"+str(sigma)+".pickle"

import os
if (os.path.isfile(out2)) and (os.path.isfile(out4)) and (os.path.isfile(out3)):
   print ("Already Processed"+str(outfix)+"_s-"+str(sigma))
   os._exit(0)


model = "/data/allenlab/scarlett/software/cmdstan/examples/"+str(model)+"/"+str(model)
tmpFile= "tmp_output."+sys.argv[3]+".txt"
initFile = "initialization_stan."+sys.argv[3]+".txt"
outFile= "stan_output."+sys.argv[3]+".txt"
input_file=in_path+inFile
tmp_output_file=out_path+tmpFile
init_file=out_path+initFile
stan_output_file=out_path+outFile

ALPHA=0.05
#
#model_theta_list = []  # 150,000
model_theta_med = []   # 150
model_med_prob = []    # 150
model_med_prob2 = []
#   
baseline_theta_list=[] # 150

with open(input_file,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        ID=fields[0]
        #print(ID)
        #med_base_theta, true_theta = getBaseline(fields)
        #baseline_theta_list.append(med_base_theta)
        med,prob,prob2=runVariant(model,fields,input_file,tmp_output_file,stan_output_file,init_file,sigma)
        if med is not None:
            #model_theta_list.extend(thetas)
            model_theta_med.append(med)
            model_med_prob.append(prob)
            model_med_prob2.append(prob2)

#pickle.dump(model_theta_list,open(out_path+"model_theta/"+str(outfix)+"_s-"+str(sigma)+".pickle",'wb'))
pickle.dump(model_theta_med,open(out2,'wb'))
pickle.dump(model_med_prob,open(out3,'wb'))
pickle.dump(model_med_prob2,open(out4,'wb'))
#pickle.dump(baseline_theta_list,open(out_path+"baseline_theta/"+str(outfix)+"_s-"+str(sigma)+".pickle",'wb'))
