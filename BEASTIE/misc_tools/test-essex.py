#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import sys

from .EssexParser import EssexParser

BASE="/Users/bmajoros/python/test/data"
filename=BASE+"/HG00096-1-subset.essex"
parser=EssexParser(filename)
while(True):
    root=parser.nextElem()
    if(not root): break
    #root.printXML(sys.stdout); print("\n")
    elem=root.pathQuery("reference-transcript/variants")
    if(elem): elem.printXML(sys.stdout); print("\n")
