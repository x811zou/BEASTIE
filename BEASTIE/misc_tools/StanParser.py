#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .Rex import Rex
from .SummaryStats import SummaryStats

rex=Rex()

######################################################################
# Attributes:
#
# Methods:
#    parser=StanParser(filename)
#    samples=parser.getSamples()
#    varNames=parser.getVarNames()
#    var=parser.getVariable(name)
#    (median,mean,SD,min,max)=parser.getSummary(var)
#    (CI_left,CI_right)=parser.getCredibleInterval(0.95,variableName)
#    (median,CI_left,CI_right)=parser.getMedianAndCI(0.95,variableName)
#    P=parser.getLeftTail(variableName,value)
#    P=parser.getRightTail(variableName,value)
######################################################################

class StanParser:
    def __init__(self,filename):
        self.samples=[]
        self.varNames=[]
        self.varIndex={}
        self.parse(filename)

    def getCredibleInterval(self,percent,name):
        samples=self.getVariable(name)
        n=len(samples)
        half=(1.0-percent)/2.0
        samples.sort(key=lambda x: x)
        CI_left=samples[int(half*n)]
        CI_right=samples[n-int(half*n)]
        return (CI_left,CI_right)

    def getLeftTail(self,name,value):
        samples=self.getVariable(name)
        count=0
        for x in samples:
            if(x<value): count+=1
        return float(count)/float(len(samples))

    def getRightTail(self,name,value):
        samples=self.getVariable(name)
        count=0
        for x in samples:
            if(x>value): count+=1
        return float(count)/float(len(samples))

    def getMedianAndCI(self,percent,name):
        samples=self.getVariable(name)
        n=len(samples)
        half=(1.0-percent)/2.0
        samples.sort(key=lambda x: x)
        CI_left=samples[int(half*n)]
        CI_right=samples[n-int(half*n)]
        median=SummaryStats.median(samples)
        return (median,CI_left,CI_right)

    def parse(self,filename):
        with open(filename,"rt") as IN:
            return self.parseFile(IN)

    def parseFile(self,IN):
        firstIndex=None
        for line in IN:
            if(len(line)<1): continue
            if(line[0]=="#"): continue
            fields=line.rstrip().split(",")
            if(len(fields)<1): continue
            if(fields[0]=="lp__"):
                firstIndex=self.getFirstVariableIndex(fields)
                self.parseVarNames(fields,firstIndex)
            else: self.parseSample(fields,firstIndex)

    def getFirstVariableIndex(self,fields):
        for i in range(len(fields)):
            field=fields[i]
            L=len(field)
            lastChar=field[L-1]
            if(lastChar!="_"): return i
        return -1

    def parseVarNames(self,fields,firstIndex):
        self.varNames=fields[firstIndex:]
        #print("firstIndex=",firstIndex,"names=",self.varNames)
        for i in range(len(self.varNames)):
            self.varIndex[self.varNames[i]]=i

    def parseSample(self,fields,firstIndex):
        sample=fields[firstIndex:]
        for i in range(len(sample)):
            sample[i]=float(sample[i])
            #print("sample=",sample[i],"firstIndex=",firstIndex)
        self.samples.append(sample)

    def getSamples(self):
        return self.samples

    def getVarNames(self):
        return self.varNames

    def getVariable(self,name):
        i=self.varIndex.get(name,None)
        if(i is None): raise Exception("Cannot find variable: "+name)
        x=[]
        for sample in self.samples:
            x.append(sample[i])
        return x

    def getSummary(self,name):
        v=self.getVariable(name)
        med=SummaryStats.median(v)
        (mean,SD,min,max)=SummaryStats.summaryStats(v)
        return (med,mean,SD,min,max)

