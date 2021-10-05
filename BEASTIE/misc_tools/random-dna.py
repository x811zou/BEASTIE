#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import random
import sys

from . import ProgramName
from .FastaWriter import FastaWriter

if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <length> <id>")
L=int(sys.argv[1])
id=sys.argv[2]
seq=""
alphabet=("A","C","G","T")
for i in range(L):
    index=int(random.random()*4)
    nuc=alphabet[index]
    seq+=nuc

writer=FastaWriter()
writer.addToFasta(">"+id,seq,sys.stdout)


