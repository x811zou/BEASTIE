#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .Codon import Codon

######################################################################
# bmajoros@duke.edu 10/15/2016
#
# Iterates through the codons of a transcript.  Each codon that it
# produces will have an absolute and relative coordinate, an
# indicator of which exon contains it, and whether the codon is
# interrupted by the end of the exon.
#
# Attributes:
#   transcript : Transcript
#   exon : which exon we are processing
#   relative : current position relative to current exon's 5-prime end
#   absolute : current position relative to genomic axis
#   stopCodons : pointer to hash containing stop codon sequences
# Methods:
#   iterator=CodonIterator(transcript,axisSequence,stopCodons)
#   codon=iterator.nextCodon(); # or None if no more
#   codons=iterator.getAllCodons()
######################################################################

class CodonIterator:
    def __init__(self,transcript,axisSequenceRef,stopCodons):
        self.transcript=transcript
        self.stopCodons=stopCodons

        # Advance to the exon containing the start codon
        exons=transcript.exons
        strand=transcript.strand
        startCodon=transcript.startCodonAbsolute
        if(startCodon):
            if(len(axisSequenceRef)==0): raise Exception("empty sequence")
            transcript.loadExonSequences(axisSequenceRef,transcript.exons)
            numExons=len(exons)
            exon=None
            for i in range(0,numExons):
                exon=exons[i]
                if(exon.containsCoordinate(startCodon) or
                   exon.containsCoordinate(startCodon-1)): break
            if(not exon.containsCoordinate(startCodon) and
               not exon.containsCoordinate(startCodon-1)):
                raise Exception(transcript.getID()+"start codon not found: "+str(startCodon))
            self.exon=exon
            self.relative=(startCodon-exon.begin) if strand=="+" else \
                exon.end-startCodon
            self.absolute=startCodon
        else:
            exon=exons[0]
            self.exon=exon
            frame=exon.getFrame()
            if(not frame): raise Exception("no frame")
            add=(3-frame)%3
            if(strand=="+"):
                self.absolute=exon.getBegin()+add
                self.relative=add
            else:
                self.absolute=exon.getEnd()-1-add
                self.relative=add

    def getAllCodons(self):
        codons=[]
        while(True):
            codon=self.nextCodon()
            if(not codon): break
            codons.append(codon)
        return codons

    def nextCodon(self):
        isStopCodon=self.stopCodons
        exon=self.exon
        if(not exon): return None
        exonSeq=exon.sequence
        relative=self.relative
        absolute=self.absolute
        exonBegin=exon.begin
        exonEnd=exon.end
        exonLen=exonEnd-exonBegin
        transcript=self.transcript
        transcriptId=transcript.getID()
        strand=transcript.strand
        isInterrupted=relative>exonLen-3
        triplet=None
        thisExonContrib=3;
        if(isInterrupted):
            thisExonContrib=exonLen-relative
            nextExonContrib=3-thisExonContrib
            triplet=exonSeq[relative:relative+thisExonContrib]
            transcript=self.transcript
            thisExonOrder=exon.order
            if(thisExonOrder is None): raise Exception("exon has no order")
            exons=transcript.exons
            numExons=len(exons)
            nextExonOrder=thisExonOrder+1
            if(nextExonOrder>=numExons): self.exon=None; return None
            else:
                nextExon=exons[nextExonOrder]
                exonSeq=nextExon.sequence
                triplet+=exonSeq[0:nextExonContrib]
                if(len(triplet)!=3): ### debugging
                    realSeqLen1=len(exon.sequence)
                    exonLen1=exon.getLength()
                    realSeqLen2=len(nextExon.sequence)
                    exonLen2=nextExon.getLength()
                    raise Exception("Error in transcript "+transcriptId+": nextContrib="+str(nextExonContrib)+" triplet=\""+triplet+"\" exonLen1="+str(exonLen1)+" realLen1="+str(realSeqLen1)+" exonLen2="+str(exonLen2)+" realLen2="+str(realSeqLen2))
                self.relative=nextExonContrib
                self.absolute=nextExon.begin+nextExonContrib if strand=="+" \
                               else nextExon.end-nextExonContrib
                self.exon=nextExon
        else: # codon was not interrupted by end of exon
            triplet=exonSeq[relative:relative+3]
            if(isStopCodon.get(triplet,False)): self.exon=None
            else:
                self.relative+=3
                self.absolute+=3 if strand=="+" else -3
                if(self.relative>=exon.end-exon.begin):
                    transcript=self.transcript
                    thisExonOrder=exon.order
                    exons=transcript.exons
                    numExons=len(exons)
                    nextExonOrder=thisExonOrder+1
                    if(nextExonOrder>=numExons):
                        self.exon=None # end of iteration
                    else:
                        nextExon=exons[nextExonOrder]
                        self.exon=nextExon
                        self.relative=0
                        self.absolute=nextExon.begin if strand=="+" \
                                       else nextExon.end
        codon=Codon(exon,triplet,relative,absolute,isInterrupted)
        codon.basesInExon=thisExonContrib
        return codon

