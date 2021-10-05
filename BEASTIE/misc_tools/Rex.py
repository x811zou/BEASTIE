#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import re


#=========================================================================
# Attributes:
#   match : returned from re.search()
# Instance Methods:
#   rex=Rex()
#   bool=rex.find("abc(\d+)def(\d+)ghi(\d+)",line)
#   rex.findOrDie("abc(\d+)def(\d+)ghi(\d+)",line)
#   x=rex[1]; y=rex[2]; z=rex[3]
#=========================================================================
class Rex:
    """Rex -- more compact regular expression matching similar to Perl"""

    def __init__(self):
        match=None

    def find(self,pattern,line):
        self.match=re.search(pattern,line)
        return self.match is not None

    def split(self,pattern,line):
        fields=re.split(pattern,line)
        nonEmpty=[]
        for x in fields:
            if(x!=""): nonEmpty.append(x)
        return nonEmpty

    def findOrDie(self,pattern,line):
        if(not self.find(pattern,line)): raise Exception("can't parse: "+line)

    def __getitem__(self,index):
        return self.match.group(index)

def test_regex():
    rex=Rex()
    line="chr1 HAVANA  initial-exon    34384   34457   .       -       0       transcript_id=ENST00000361813.5;gene_id=ENSG00000198952.7;\n"
    result=rex.find('transcript_id[:=]?\s*"?([^\s";]+)"?',line)
    #result=rex.find('transcript_id=([^\s";]+)',line)
    print(result)
    #x=y=z=None
    #if(rex.find("abc(\d+)abc(\d+)abc","ab123abc456abc789")):
    #    x=rex[1]; y=rex[2]
    #elif(rex.find("dog(\d+)cat(\d+)cow(\d+)chicken(\d+)",
    #              "dog1cat2cow8chicken100")):
    #    x=rex[1]; y=rex[4]
    #print(x,y)



