#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# 2018 William H. Majoros (bmajoros@allumni.duke.edu)
#=========================================================================
import gzip

from .Rex import Rex

rex=Rex()


#=========================================================================
# Attributes:
#   fh : file handle
# Instance Methods:
#   reader=FastqReader(filename) # can be gzipped!
#   (ID,seq,qual,qualSeq,pair)=reader.nextSequence() # returns None at EOF
#        * pair indicates which read of the pair: 1 or 2
#        * qual is an array of integer quality values
#        * qualSeq is the raw quality string
#   reader.close()
# Class Methods:
#=========================================================================
class FastqReader:
    """FastqReader"""
    def __init__(self,filename):
        if(filename is not None):
            if(rex.find("\.gz$",filename)): self.fh=gzip.open(filename,"rt")
            else: self.fh=open(filename,"r")

    def close(self):
        self.fh.close()

    def nextSequence(self):
        fh=self.fh
        line=fh.readline()
        if(line is None): return None
        if(len(line)==0): return None
        if(not rex.find("^(\S+)",line)):
            return None
            #raise Exception("Cannot parse fastq line: "+ID)
        ID=rex[1]
        pair=1
        if(rex.find("\s+(\d)",line)): pair=int(rex[1])
        seq=fh.readline().rstrip()
        junk=fh.readline()
        qualSeq=fh.readline().rstrip()
        qual=[ord(x)-33 for x in qualSeq]
        return (ID,seq,qual,qualSeq,pair)



