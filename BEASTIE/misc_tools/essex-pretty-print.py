#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
import sys

from EssexParser import EssexParser

from . import ProgramName

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <in.essex>\n")
(infile,)=sys.argv[1:]

parser=EssexParser(infile)
while(True):
    tree=parser.nextElem()
    if(tree is None): break
    tree.print(sys.stdout)




