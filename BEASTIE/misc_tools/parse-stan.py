#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
import sys

from . import ProgramName
from .StanParser import StanParser

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)<3):
    exit(ProgramName.get()+" <infile.txt> <var1> <var2> ...\n")
infile=sys.argv[1]
variables=sys.argv[2:]

parser=StanParser(infile)
#(median,mean,SD,min,max)=parser.getSummary(variable)
#print("# posterior median=",median,sep="")
samplesByVar=[]
n=0
for var in variables:
    samples=parser.getVariable(var)
    samplesByVar.append(samples)
    n=len(samples)
for i in range(n):
    line=[]
    for sample in samplesByVar: line.append(str(sample[i]))
    print("\t".join(line))



