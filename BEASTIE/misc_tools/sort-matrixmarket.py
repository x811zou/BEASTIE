#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
import gzip
import sys

from . import ProgramName, TempFilename
from .Pipe import Pipe

HEADERFILE=TempFilename.generate(".header")
SORTEDFILE=TempFilename.generate(".sorted")

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" in.mtx.gz column out.mtx.gz\n    where column = 1-based field index\n")
(infile,index,outfile)=sys.argv[1:]

OUT=open(HEADERFILE,"wt")
numHeader=1
with gzip.open(infile,"rt") as IN:
    for line in IN:
        if(len(line)==0): raise Exception("unexpected empty line")
        if(line[0]=="%"):
            numHeader+=1
            print(line,file=OUT,end="")
        else:
            print(line,file=OUT,end="")
            break
IN.close(); OUT.close()

cmd="cat "+infile+" | gunzip | tail -n +"+str(numHeader+1)+\
    " | sort -g -k "+index+" > "+SORTEDFILE
Pipe.run(cmd)
cmd="cat "+HEADERFILE+" "+SORTEDFILE+" | gzip > "+outfile
Pipe.run(cmd)
Pipe.run("rm "+SORTEDFILE)
Pipe.run("rm "+HEADERFILE)
