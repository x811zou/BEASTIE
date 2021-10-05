#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from enum import Enum


#=========================================================================
# Attributes:
#   FORWARD : int
#   REVERSE : int
# Instance Methods:
#   Strand()
# Class Methods:
#   toString(Strand)
#=========================================================================
class Strand(Enum):
    FORWARD=1
    REVERSE=0

    @classmethod
    def toString(cls,strand):
        return "+" if strand==strand.FORWARD else "-"




