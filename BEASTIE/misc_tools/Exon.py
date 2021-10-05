#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .Interval import Interval

######################################################################
#
# Exon.pm
#
# bmajoros@duke.edu 10/15/2016
#
# A pair of coordinates representing the location of an exon on a
# genomic axis.  The first coordinate (begin) is always less than
# the second coordinate (end), regardless of strand.  The sequence
# is only loaded if loadExonSequences() or loadTranscriptSeq() is sent
# to the parent transcript.  The sequence is automatically reverse-
# complemented if on the reverse strand.  Has an order index telling
# you which exon in the transcript this is, from the first (0) to
# the last, in translation order.  Thus, the first exon (order=0
# is actually rightmost when on the reverse strand.  Note that frame
# (actually phase) denotes the phase of the first base of an exon, but
# on the reverse strand, the first base is actually the rightmost.
#
# Attributes:
#   begin : the left end of the exon (zero-based)
#   end   : one base past the right end of the exon (zero-based)
#   order : which exon this is (0,1,2,...), in translation order
#   sequence : NOT LOADED BY DEFAULT! (but automatically reverse-
#              complemented when it is loaded)
#   transcript : the parent transcript
#   frame : frame of first base in exon
#   type : type of object, usually "exon" or "internal-exon", etc...
#   score : stored as a string so that "." denotes "no score"
#   strand : string
#   substrate : string
# Methods:
#   exon=Exon(begin,end,transcript)
#   exon.containsCoordinate(x) : boolean
#   new=exon.copy()
#   length=exon.getLength()
#   exon.reverseComplement(seqLen)
#   exon.trimInitialPortion(numBases)
#   exon.trimFinalPortion(numBases)
#   strand=exon.getStrand() # "+" or "-"
#   bool=exon.overlaps(otherExon)
#   sequence=exon.getSequence()
#   transcript=exon.getTranscript()
#   gff=exon.toGff()
#   begin=exon.getBegin()
#   end=exon.getEnd()
#   interval=exon.asInterval()
#   frame=exon.getFrame()
#   exon.setFrame(frame)
#   type=exon.getType()
#   exon.setType(type)
#   exon.setScore(score)
#   score=exon.getScore()
#   exon.shiftCoords(delta)
#   exon.setStrand(strand)
#   substrate=exon.getSubstrate()
#   exon.setSubstrate(substrate)
#   exon.setBegin(begin)
#   exon.setEnd(end)
######################################################################

class Exon:
    def __init__(self,begin,end,transcript):
        self.begin=begin
        self.end=end
        self.transcript=transcript
        self.score="."
        self.type="exon"
        self.strand=transcript.getStrand() if transcript else None
        self.order=None
        self.frame=None

    def setBegin(self,begin):
        self.begin=begin

    def setEnd(self,end):
        self.end=end

    def containsCoordinate(self,x):
        return x>=self.begin and x<self.end

    def asInterval(self):
        return Interval(self.begin,self.end)

    def getLength(self):
        return self.end-self.begin

    def trimInitialPortion(self,numBases):
        """Trims a certain number of bases from the translationally early
        part of the exon.
        """
        if(self.getStrand()=="+"): self.begin+=numBases
        else: self.end-=numBases
        sequence=self.sequence
        if(sequence is not None):
            self.sequence=sequence[numBases:len(sequence)-numBases]

    def trimFinalPortion(self,numBases):
        """Trims a certain number of bases from the translationally late
        part of the exon
        """
        if(self.getStrand()=="+"): self.end-=numBases
        else: self.begin+=numBases
        sequence=self.sequence
        if(sequence):
            self.sequence=sequence[0:len(sequence)-numBases]

    def getStrand(self):
        return self.transcript.strand

    def overlaps(self,otherExon):
        return self.begin<otherExon.end and otherExon.begin<self.end

    def getSequence(self):
        return self.sequence

    def getTranscript(self):
        return self.transcript

    def getBegin(self):
        return self.begin

    def getEnd(self):
        return self.end

    def getFrame(self):
        frame=self.frame
        return frame if frame is not None else "."

    def getType(self):
        return self.type

    def toGff(self):
        transcript=self.getTranscript()
        substrate=self.getSubstrate()
        begin=self.getBegin()+1 # convert to 1-based coordinate system (1/B)
        end=self.getEnd()
        type=self.getType()
        if(not type): type="."
        source=transcript.getSource() if transcript else "."
        if(not source): source="."
        strand=self.getStrand()
        frame=self.getFrame()
        score=self.getScore()
        transcriptId=transcript.getID() if transcript else "."
        geneId=transcript.getGeneId() if transcript else "."
        return "\t".join([substrate,source,type,str(begin),str(end),score,
                          strand,str(frame),
                          "transcript_id \""+transcriptId+"\"; gene_id \""+
                          geneId+"\"\n"])

    def setScore(self,score):
        self.score=score

    def getScore(self):
        score=self.score
        return score if score else "."

    def setType(self,type):
        self.type=type

    def shiftCoords(self,delta):
        self.begin+=delta
        self.end+=delta

    def setFrame(self,frame):
        self.frame=frame

    def getStrand(self):
        strand=self.strand
        if(strand): return strand
        return self.transcript.getStrand()

    def setStrand(self,strand):
        self.strand=strand

    def getSubstrate(self):
        substrate=self.substrate if hasattr(self,'substrate') else None
        if(substrate): return substrate
        return self.transcript.getSubstrate()

    def setSubstrate(self,substrate):
        self.substrate=substrate

    def compStrand(self,strand):
        if(strand=="+"): return "-"
        if(strand=="-"): return "+"
        if(strand=="."): return "."
        raise Exception("Unknown strand \""+strand+"\"")

    def reverseComplement(self,seqLen):
        begin=self.getBegin()
        end=self.getEnd()
        self.begin=seqLen-end
        self.end=seqLen-begin
        self.strand=self.compStrand(self.strand)

    def copy(self):
        new=Exon(self.begin,self.end,self.transcript)
        new.score=self.score
        new.type=self.type
        new.strand=self.strand
        return new
