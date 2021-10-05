#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# 2018 William H. Majoros (bmajoros@allumni.duke.edu)
#=========================================================================
import gzip

from .CigarString import CigarString
from .Rex import Rex
from .SamRecord import SamRecord

rex=Rex()

#=========================================================================
# Attributes:
#   fh : file handle
#   headerLines : array of header lines
# Instance Methods:
#   reader=SamReader(filename)
#   samRecord=reader.nextSequence() # returns None at EOF
#   (record,line)=reader.nextSeqAndText() # returns None at EOF
#   reader.close()
# Class Methods:
#=========================================================================
class SamReader:
    """SamReader"""
    def __init__(self,filename):
        self.headerLines=[]
        if(filename is not None):
            if(rex.find("\.gz$",filename)): self.fh=gzip.open(filename,"rt")
            else: self.fh=open(filename,"r")

    def close(self):
        self.fh.close()

    def nextSequence(self):
        pair=self.nextSeqAndText()
        if(pair is None): return None
        (rec,line)=pair
        return rec

    def nextSeqAndText(self):
        headerChars=set(["@","["])
        fh=self.fh
        line=fh.readline()
        if(line is None): return None
        while(line is not None and len(line)>0 and line[0] in headerChars):
            if(line[0]=="@"): self.headerLines.append(line)
            line=fh.readline()
        if(line is None or len(line)==0): return None
        fields=line.rstrip().split()
        if(len(fields)<11): raise Exception("can't parse sam line: "+line)
        (ID,flags,refName,refPos,mapQual,cigar,rnext,pnext,templateLen,
         seq,qual)=fields[:11]
        refPos=int(refPos)-1 # convert 1-based to 0-based
        flags=int(flags)
        CIGAR=CigarString(cigar)
        tags=fields[11:]
        rec=SamRecord(ID,refName,refPos,CIGAR,seq,flags,tags)
        return (rec,line)

# M03884:303:000000000-C4RM6:1:1101:1776:15706    99      chrX:31786371-31797409  6687    44      150M    =       6813    271     ATACTATTGCTGCGGTAATAACTGTAACTGCAGTTACTATTTAGTGATTTGTATGTAGATGTAGATGTAGTCTATGTCAGACACTATGCTGAGCATTTTATGGTTGCTATGTACTGATACATACAGAAACAAGAGGTACGTTCTTTTACA  BBBBFFFFFFFGGGGGEFGGFGHFHFFFHHHFFHHHFHFHHHGFHEDGGHFHBGFHGBDHFHFFFHHHHFHHHHHGHGFFBGGGHFHFFHHFFFFHHHHGHGFHHGFHGHHHGFHFFHHFHHFFGFFFFGGEHFFEHHFGHHHGHHHHFB  AS:i:300        XN:i:0

