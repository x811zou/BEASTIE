#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys

from . import ProgramName

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <infile.txt>\n")
(infile,)=sys.argv[1:]

lines=[]
with open(infile,"rt") as IN:
    for line in IN:
        lines.append(line)
for i in range(len(lines)-1,-1,-1):
    print(lines[i],end="")


