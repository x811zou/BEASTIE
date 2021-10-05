#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================

#=========================================================================
# Attributes:
#
# Instance Methods:
#   SamMDtagParser()
# Class Methods:
#   N=SamMDtagParser.countMismatches(MDtags)
#=========================================================================
class SamMDtagParser:
    """SamMDtagParser"""
    def __init__(self):
        pass

    @classmethod
    def countMismatches(cls,MDtags):
        mismatches=0
        L=len(MDtags)
        for i in range(L):
            if(i%2==1):
                if(len(MDtags[i])==1): mismatches+=1
        return mismatches



