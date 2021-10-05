#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .Interval import Interval


#=========================================================================
# Attributes:
#   chr : string
#   interval : Interval
# Instance Methods:
#   record=Bed3Record(chr,begin,end)
#   bool=record.isBed3()
#   bool=record.isBed6()
#   begin=record.getBegin()
#   end=record.getEnd()
#   line=record.toString()
# Class Methods:
#
#=========================================================================
class Bed3Record:
    """Bed3Record represents a record in a BED3 file"""
    def __init__(self,chr,begin,end):
        self.chr=chr
        self.interval=Interval(begin,end)

    def isBed3(self):
        return True

    def isBed6(self):
        return False

    def getBegin(self):
        return self.interval.begin

    def getEnd(self):
        return self.interval.end

    def toString(self):
        return self.chr+"\t"+str(self.interval.begin)+"\t"+\
            str(self.interval.end)
