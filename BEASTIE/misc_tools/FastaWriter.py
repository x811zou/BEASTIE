#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import os
import re


#=========================================================================
# Attributes:
#   width
# Methods:
#   FastaWriter()
# Instance Methods:
#   writer=FastaWriter(optionalWidth)
#   writer.writeFasta(defline,sequence,filename)
#   writer.appendToFasta(defline,sequence,filename)
#   writer.addToFasta(defline,sequence,filehandle)
#=========================================================================
class FastaWriter:
    """FastaWriter"""
    def __init__(self,width=60):
        self.width=width

    def writeFasta(self,defline,seq,filename):
        with open(filename,"w") as fh:
            self.addToFasta(defline,seq,fh)

    def addToFasta(self,defline,seq,fh):
        defline=defline.rstrip()
        if(not re.search("^\s*>",defline)): defline=">"+defline
        fh.write(defline+"\n");
        length=len(seq)
        numLines=length//self.width
        if(length%self.width>0): numLines+=1
        start=0
        for i in range(0,numLines):
            line=seq[start:start+self.width]
            fh.write(line+"\n")
            start+=self.width
        if(length==0): fh.write("\n")

    def appendToFasta(self,defline,seq,filename):
        if(not os.path.exists(filename)):
            self.writeFasta(defline,seq,filename)
            return
        with open(filename,"a") as fh:
            self.addToFasta(defline,seq,fh)


