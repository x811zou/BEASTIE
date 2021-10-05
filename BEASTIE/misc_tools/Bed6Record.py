#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .Bed3Record import Bed3Record


#=========================================================================
# Inherited Attributes:
#   chr : string
#   interval : Interval
# Attributes:
#   name : string
#   score : float
#   strand : string
# Instance Methods:
#   record=Bed6Record(chr,begin,end,name,score,strand)
#   bool=isBed3()
#   bool=isBed6()
#   str=toString()
# Class Methods:
#
#=========================================================================
class Bed6Record(Bed3Record):
    """Bed6Record represents a record in a BED6 file"""
    def __init__(self,chr,begin,end,name,score,strand):
        Bed3Record.__init__(self,chr,begin,end)
        self.name=name
        self.score=score
        self.strand=strand

    def isBed3(self):
        return False

    def isBed6(self):
        return True

    def toString(self):
        s=self.chr+"\t"+str(self.interval.begin)+"\t"+str(self.interval.end)\
            +"\t"+self.name
        if(self.score is not None): s+="\t"+str(self.score)
        if(self.strand is not None): s+="\t"+self.strand
        return s
