#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================

######################################################################
# Attributes:
#   transcripts : list of Transcript objects
#   ID
#   transcriptHash : transcripts hashed by their ID
# Methods:
#   gene=Gene()
#   gene.addTranscript(t)
#   n=gene.getNumTranscripts()
#   n=gene.numTranscripts()
#   t=gene.getIthTranscript(i)
#   t=gene.longestTranscript()
#   id=gene.getId()
#   id=gene.getID()
#   gene.setId(id)
#   begin=gene.getBegin() # leftmost edge
#   end=gene.getEnd()     # rightmost edge
#   strand=gene.getStrand()
#   substrate=gene.getSubstrate()
#   gff=gene.toGff()
#   exons=gene.getMergedExons()
#
######################################################################

class Gene:
    def __init__(self):
        self.transcripts=[]
        self.transcriptHash={}

    def getMergedExons(self):
        transcripts=self.transcripts
        exons=[]
        for transcript in transcripts:
            raw=transcript.getRawExons()
            exons.extend(raw)
            #print("RAW:",len(raw))
            #for i in range(len(raw)):
                #print("\t",raw[i].begin,raw[i].end)
            #print()
        exons.sort(key=lambda x: x.begin)
        n=len(exons)
        i=0
        while(i<n-1):
            if(exons[i].overlaps(exons[i+1])):
                exons[i].end=max(exons[i].end,exons[i+1].end)
                del exons[i+1]
                n-=1
            else: i+=1
        return exons

    def getStrand(self):
        transcripts=self.transcripts
        transcript=transcripts[0]
        return transcript.getStrand()

    def getSubstrate(self):
        transcripts=self.transcripts
        transcript=transcripts[0]
        return transcript.getSubstrate()

    def getBegin(self):
        transcripts=self.transcripts
        if(not transcripts):
            raise Exception("no transcripts in gene "+self.getID())
        begin=None
        for transcript in transcripts:
            b=transcript.getBegin()
            if(b is None):
                raise Exception("transcript has no begin: "+transcript.getID())
            if(begin is None or b<begin): begin=b
        return begin

    def getEnd(self):
        transcripts=self.transcripts
        end=None
        for transcript in transcripts:
            e=transcript.getEnd()
            if(end is None or e>end): end=e
        return end

    def addTranscript(self,transcript):
        id=transcript.getTranscriptId()
        hash=self.transcriptHash
        if(hash.get(id,None) is not None): return
        self.transcripts.append(transcript)
        hash[id]=transcript

    def getNumTranscripts(self):
        return len(self.transcripts)

    def numTranscripts(self):
        return len(self.transcripts)

    def getIthTranscript(self,i):
        return self.transcripts[i]

    def longestTranscript(self):
        transcripts=self.transcripts
        if(len(transcripts)==0): return None
        longest=transcripts[0]
        longestLength=longest.getExtent()
        for transcript in transcripts[1:]:
            length=transcript.getExtent()
            if(length>longestLength):
                longest=transcript
                longestLength=length
        return longest

    def getId(self):
        return self.ID

    def getID(self):
        return self.ID

    def setId(self,id):
        self.ID=id

    def getBeginAndEnd(self):
        transcripts=self.transcripts
        begin=None
        end=None
        for transcript in transcripts:
            b=transcript.getBegin()
            e=transcript.getEnd()
            if(begin is None): begin=b; end=e
            else:
                if(b<begin): begin=b
                if(e>end): end=e
        return (begin,end)

    def toGff(self):
        transcripts=self.transcripts
        gff=""
        for transcript in transcripts:
            gff+=transcript.toGff()
        return gff

    def __hash__(self):
        return hash(self.ID)

    def __eq__(self,other):
        return self.ID==other.ID


