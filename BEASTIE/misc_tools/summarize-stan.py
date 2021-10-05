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
print("variable\tmedian\tSD")
for var in variables:
    (median,mean,SD,min,max)=parser.getSummary(var)
    print(var,round(median,4),round(SD,5),sep="\t")


