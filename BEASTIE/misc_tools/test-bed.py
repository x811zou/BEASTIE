#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .BedReader import BedReader

BASE="/Users/bmajoros/python/test/data"
filename=BASE+"/DEGs_downreg.FDR_0.1.TSS.protein_coding.bed"

reader=BedReader(filename)
while(True):
    rec=reader.nextRecord()
    if(rec is None): break
    begin=rec.getBegin()
    end=rec.getEnd()
    print(begin,"-",end,sep="",end="")
    print("\t"+rec.name+"\t"+str(rec.score)+"\t"+rec.strand)
reader.close()
