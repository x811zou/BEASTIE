#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================

######################################################################
# bmajoros@duke.edu 10/15/2016
#
# Represents a single trinucleotide at a certain place within an exon.
# Stores its location both relative to the transcript and to the
# genomic axis.  All coordinates are zero-based and space-based.  On
# the forward strand, the absoluteCoord represents the position of the
# leftmost base in the codon, whereas on the reverse strand,
# the absoluteCoord points to the base immediately following the third
# base of the codon.  Thus, on the minus strand, a substring operation
# on the genomic axis should begin at absoluteCoord-3.  RelativeCoords
# are much simpler: a substring operation on the transcript always
# uses the relativeCoord as-is, regardless of strand.
#
# Attributes:
#   string triplet
#   int absoluteCoord : relative to genomic axis
#   int relativeCoord : relative to current exon
#   bool isInterrupted : exon ends before codon is complete
#   int basesInExon :  how many bases of this codon are in this exon (1-3)
#   Exon exon : which exon contains this codon
# Methods:
#   codon=Codon(exon,triplet,relative,absolute,isInterrupted)
######################################################################

class Codon:
    def __init__(self,exon,triplet,relative,absolute,isInterrupted):
        if(relative is None): raise Exception("relative is not set")
        self.triplet=triplet
        self.exon=exon
        self.relativeCoord=relative
        self.absoluteCoord=absolute
        self.isInterrupted=isInterrupted


