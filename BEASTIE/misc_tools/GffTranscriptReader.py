#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import re

from .Exon import Exon
from .Gene import Gene
from .Integer import Integer
from .Rex import Rex
from .Transcript import Transcript

######################################################################
# Returns a list of Transcripts.  For each transcript, the Exons
# will be sorted according to order of translation, so that
# the exon containing the start codon will come before the exon
# containing the stop codon.  This means that for the minus strand,
# exon begin coordinates will be decreasing.  However, an individual
# exon's begin coordinate is always less than its end coordinate.
# The transcripts themselves are sorted along the chromosome left-to-
# right.  Note that although the GFF coordinates are 1-based/base-based
# (1/B), internally all coordinates are stored as 0-based/space-based
# (0/S).  This conversion is handled automatically.
#
# Attributes:
#   shouldSortTranscripts
#   exonsAreCDS : interpret "exon" features as "CDS" when reading GFF
# Methods:
#   reader=GffTranscriptReader()
#   reader.setStopCodons({"TAG":1,"TAA":1,"TGA":1})
#   transcriptArray=reader.loadGFF(filename)
#   geneList=reader.loadGenes(filename)
#   hashTable=reader.hashBySubstrate(filename)
#   reader.hashBySubstrateInto(filename,hash)
#   hashTable=reader.hashGenesBySubstrate(filename)
#   hashTable=reader.loadTranscriptIdHash(filename)
#   hashTable=reader.loadGeneIdHash(filename)
#   reader.doNotSortTranscripts()
######################################################################

class GffTranscriptReader:
    def __init__(self):
        self.shouldSortTranscripts=True
        self.exonsAreCDS=False
        self.stopCodons={"TAG":1,"TAA":1,"TGA":1}

    def loadGenes(self,filename):
        transcripts=self.loadGFF(filename)
        genes=set()
        for transcript in transcripts:
            gene=transcript.getGene()
            if(not gene):
                raise Exception("transcript "+transcript.getID()+
                                " has no gene")
            genes.add(gene)
        genes=list(genes)
        genes.sort(key=lambda gene: gene.getBegin())
        return genes

    def doNotSortTranscripts(self):
        self.shouldSortTranscripts=False

    def loadTranscriptIdHash(self,filename):
        transcriptArray=self.loadGFF(filename)
        hash={}
        for transcript in transcriptArray:
            id=transcript.getID()
            hash[id]=transcript
        return hash

    def loadGeneIdHash(self,filename):
        transcriptArray=self.loadGFF(filename)
        hash={}
        for transcript in transcriptArray:
            id=transcript.getGeneId()
            array=hash.get(id,None)
            if(array is None): array=hash[id]=[]
            array.append(transcript)
        return hash

    def hashBySubstrate(self,filename):
        hash={}
        self.hashBySubstrateInto(filename,hash)
        return hash

    def hashBySubstrateInto(self,filename,hash):
        transcriptArray=self.loadGFF(filename)
        for transcript in transcriptArray:
            id=transcript.getSubstrate()
            array=hash.get(id,None)
            if(array is None): array=hash[id]=[]
            array.append(transcript)

    def hashGenesBySubstrate(self,filename):
        geneArray=self.loadGenes(filename)
        hash={}
        for gene in geneArray:
            id=gene.getSubstrate()
            array=hash.get(id,None)
            if(array is None): array=hash[id]=[]
            array.append(gene)
        return hash

    def computeFrames(self,transcripts):
        for transcript in transcripts:
            if(not transcript.areExonTypesSet()): transcript.setExonTypes()
            strand=transcript.getStrand()
            exons=transcript.exons
            frame=0 # fine for both strands
            for exon in exons:
                exon.frame=frame
                length=exon.getLength()
                frame=(frame+length)%3 # this is fine, on both strands

    def exonContainsPoint(self,exon,point):
        return point>=exon.begin and point<=exon.end;
        # this '<=' is necessary for the minus strand!
        # do not change it back to '<' !!!

    def setStopCodons(self,stopCodons):
        self.stopCodons=stopCodons

    def adjustStartCodons_fw(self,transcript,totalIntronSize):
        startCodon=None
        exons=transcript.exons
        exons.sort(key=lambda exon: exon.begin)
        numExons=len(exons)
        #if(transcript.getID()=="ENST00000361390.2"):
        #    print("adjustStartCodons_fw",numExons,"exons")
        if(numExons==0): return None
        if(transcript.begin is None) :
            transcript.begin=exons[0].begin
        if(transcript.end is None):
            transcript.end=exons[numExons-1].end
        if(transcript.startCodon is not None):
            startCodon=transcript.startCodon-transcript.begin
        else:
            startCodon=0
            transcript.startCodon=transcript.begin
            transcript.startCodonAbsolute=transcript.begin
        for i in range(numExons): exons[i].order=i
        for i in range(numExons):
            exon=exons[i]
            if(i>0):
                prevExon=exons[i-1]
                intronSize=exon.begin-prevExon.end
                totalIntronSize+=intronSize
            if(transcript.startCodon is not None and
               self.exonContainsPoint(exon,transcript.startCodon)): break
        return startCodon

    def adjustStartCodons_bw(self,transcript,totalIntronSize):
        startCodon=None
        exons=transcript.exons
        exons.sort(key=lambda exon: exon.begin,reverse=True)
        numExons=len(exons)
        if(numExons==0): return None
        if(transcript.end is None): transcript.end=exons[0].end
        if(transcript.begin is None): transcript.begin=exons[numExons-1].begin
        if(transcript.startCodon is not None):
            startCodon=transcript.end-transcript.startCodon
        else:
            startCodon=0
            transcript.startCodon=transcript.end ###
            transcript.startCodonAbsolute=transcript.end ###
        for i in range(numExons):
            exon=exons[i]
            exon.order=i
            if(i>0):
                prevExon=exons[i-1]
                intronSize=prevExon.begin-exon.end
                totalIntronSize+=intronSize
            if(transcript.startCodon is not None and
               self.exonContainsPoint(exon,transcript.startCodon)): break
        return startCodon

    def adjustStartCodons(self,transcripts):
        for transcript in transcripts:
            transcript.sortExons()
            transcript.adjustOrders()
            strand=transcript.strand
            startCodon=None
            totalIntronSize=Integer(0)
            if(strand=="+"):
                startCodon=\
                    self.adjustStartCodons_fw(transcript,totalIntronSize)
            else:
                startCodon=\
                    self.adjustStartCodons_bw(transcript,totalIntronSize)
            if(startCodon is not None):
                startCodon-=int(totalIntronSize)
                transcript.startCodon=startCodon

    def loadGFF_transcript(self,fields,line,transcriptBeginEnd,GFF,
                           transcripts,readOrder,genes):
        begin=int(fields[3])-1
        end=int(fields[4])
        rex=Rex()
        if(rex.find('transcript_id[:=]?\s*"?([^\s";]+)"?',line)):
            transcriptId=rex[1]
            transcriptBeginEnd[transcriptId]=[begin,end]
            strand=fields[6]
            score=fields[5]
            transcriptExtraFields=""
            for i in range(8,len(fields)):
                transcriptExtraFields+=fields[i]+" "
            transcript=transcripts.get(transcriptId,None)
            if(transcript is None):
                transcripts[transcriptId]=transcript= \
	                                   Transcript(transcriptId,strand)
                transcript.setStopCodons(self.stopCodons)
                transcript.readOrder=readOrder;
                readOrder+=1
                transcript.substrate=fields[0]
                transcript.source=fields[1]
                transcript.setBegin(begin)
                transcript.setEnd(end)
            if(transcript.score is None and
               score!="."): transcript.score=float(score)
            geneId=None
            if(rex.find("genegrp=(\S+)",line)): geneId=rex[1]
            elif(rex.find('gene_id[:=]?\s*\"?([^\s\;"]+)\"?',line)):
                geneId=rex[1]
            if(not geneId): raise Exception("can't parse GTF: "+line)
            transcript.geneId=geneId
            gene=genes.get(geneId,None)
            if(not gene): genes[geneId]=gene=Gene(); gene.setId(geneId)
            transcript.setGene(gene)
            gene.addTranscript(transcript)
            transcript.extraFields=transcriptExtraFields

    def loadGFF_UTR(self,fields,line,transcriptBeginEnd,GFF,
                           transcripts,readOrder,genes):
        exonBegin=int(fields[3])-1
        exonEnd=int(fields[4])
        exonScore=fields[5]
        strand=fields[6]
        frame=fields[7]
        transcriptId=None
        rex=Rex()
        if(rex.find('transgrp[:=]\s*(\S+)',line)): transcriptId=rex[1]
        elif(rex.find('transcript_id[:=]?\s*"?([^\s";]+)"?',line)):
            transcriptId=rex[1]
        elif(rex.find('Parent=([^;,\s]+)',line)): transcriptId=rex[1]
        geneId=None
        if(rex.find('genegrp=(\S+)',line)): geneId=rex[1]
        elif(rex.find('gene_id[:=]?\s*"?([^\s\;"]+)"?',line)): geneId=rex[1]
        if(transcriptId is None): transcriptId=geneId
        if(geneId is None): geneId=transcriptId
        if(transcriptId is None):
            raise Exception(line+" : no transcript ID found")
        if(rex.find("(\S+);$",transcriptId)): transcriptId=rex[1]
        if(rex.find("(\S+);$",geneId)): geneId=rex[1]
        extra=""
        for i in range(8,len(fields)): extra+=fields[i]+" "
        if(exonBegin>exonEnd): (exonBegin,exonEnd)=(exonEnd,exonBegin)
        transcript=transcripts.get(transcriptId,None)
        if(not transcript):
            transcripts[transcriptId]=transcript= \
                Transcript(transcriptId,strand)
            transcript.setStopCodons(self.stopCodons)
            transcript.readOrder=readOrder
            readOrder+=1
            transcript.substrate=fields[0]
            transcript.source=fields[1]
            if(transcriptBeginEnd.get(transcriptId,None) is not None):
                (begin,end)=transcriptBeginEnd[transcriptId]
                transcript.setBegin(begin)
                transcript.setEnd(end)
            else:
                transcript.setBegin(exonBegin)
                transcript.setEnd(exonEnd)
        transcript.geneId=geneId
        gene=genes.get(geneId,None)
        if(gene is None):
            genes[geneId]=gene=Gene(); gene.setId(geneId)
        transcript.setGene(gene)
        exon=Exon(exonBegin,exonEnd,transcript)
        exon.extraFields=extra
        if(transcript.rawExons is not None):
            exon.frame=frame
            exon.score=exonScore
            exon.type=fields[2]
            transcript.rawExons.append(exon)
        elif(not transcript.exonOverlapsExon(exon)):
            exon.frame=frame
            exon.score=exonScore
            exon.type=fields[2]
            transcript.UTR.append(exon) # OK -- we sort later
        gene.addTranscript(transcript)

    def loadGFF_exon(self,fields,line,transcriptBeginEnd,GFF,
                           transcripts,readOrder,genes):
        exonBegin=int(fields[3])-1
        exonEnd=int(fields[4])
        exonScore=fields[5]
        strand=fields[6]
        frame=fields[7]
        transcriptId=None
        rex=Rex()
        if(rex.find("transgrp[:=]\s*(\S+)",line)): transcriptId=rex[1]
        elif(rex.find('transcript_id[:=]?\s*"?([^\s";]+)"?',line)):
            transcriptId=rex[1]
        elif(rex.find('Parent=([^;,\s]+)',line)): transcriptId=rex[1]
        geneId=None
        if(rex.find("genegrp=(\S+)",line)): geneId=rex[1]
        elif(rex.find('gene_id[:=]?\s*"?([^\s\;"]+)"?',line)): geneId=rex[1]
        if(transcriptId is None): transcriptId=geneId
        if(geneId is None): geneId=transcriptId
        if(rex.find("(\S+);$",transcriptId)): transcriptId=rex[1]
        if(rex.find("(\S+);$",geneId)): geneId=rex[1]
        extra=""
        for i in range(8,len(fields)): extra+=fields[i]+" "
        if(exonBegin>exonEnd): (exonBegin,exonEnd)=(exonEnd,exonBegin)
        transcript=transcripts.get(transcriptId,None)
        if(transcript is None):
            transcripts[transcriptId]=transcript= \
                                       Transcript(transcriptId,strand)
            transcript.setStopCodons(self.stopCodons)
            transcript.readOrder=readOrder
            readOrder+=1
            transcript.substrate=fields[0]
            transcript.source=fields[1]
            if(transcriptBeginEnd.get(transcriptId,None) is not None):
                (begin,end)=transcriptBeginEnd[transcriptId]
                transcript.setBegin(begin)
                transcript.setEnd(end)
            else:
                transcript.setBegin(exonBegin)
                transcript.setEnd(exonEnd)
        transcript.geneId=geneId
        gene=genes.get(geneId,None)
        if(gene is None):
           genes[geneId]=gene=Gene(); gene.setId(geneId)
        transcript.setGene(gene)
        exon=Exon(exonBegin,exonEnd,transcript)
        exon.extraFields=extra
        exon.score=exonScore
        exon.type=fields[2]
        if(transcript.rawExons is None): transcript.rawExons=[]
        transcript.rawExons.append(exon)
        gene.addTranscript(transcript)

    def loadGFF_CDS(self,fields,line,transcriptBeginEnd,GFF,
                    transcripts,readOrder,genes):
        exonBegin=int(fields[3])-1
        exonEnd=int(fields[4])
        exonScore=fields[5]
        strand=fields[6]
        frame=fields[7]
        transcriptId=None
        rex=Rex()
        if(rex.find('transgrp[:=]\s*(\S+)',line)): transcriptId=rex[1]
        elif(rex.find('transcript_id[:=]?\s*"?([^\s";]+)"?',line)):
            transcriptId=rex[1]
        elif(rex.find('Parent=([^;,\s]+)',line)): transcriptId=rex[1]
        geneId=None
        if(rex.find('genegrp=(\S+)',line)): geneId=rex[1]
        elif(rex.find('gene_id[:=]?\s*"?([^\s\;"]+)"?',line)): geneId=rex[1]
        if(transcriptId is None): transcriptId=geneId
        if(geneId is None): geneId=transcriptId
        if(transcriptId is None):
            raise Exception(line+" : no transcript ID found")
        if(rex.find('(\S+);$',transcriptId)): transcriptId=rex[1]
        if(rex.find('(\S+);$',geneId)): geneId=rex[1]
        extra=""
        for i in range(8,len(fields)): extra+=fields[i]+" "
        if(exonBegin>exonEnd): (exonBegin,exonEnd)=(exonEnd,exonBegin)
        transcript=transcripts.get(transcriptId,None)
        if(transcript is None):
            transcripts[transcriptId]=transcript= \
	                               Transcript(transcriptId,strand)
            transcript.setStopCodons(self.stopCodons)
            transcript.readOrder=readOrder
            readOrder+=1
            transcript.substrate=fields[0]
            transcript.source=fields[1]
            if(transcriptBeginEnd.get(transcriptId,None) is not None):
                (begin,end)=transcriptBeginEnd[transcriptId]
                transcript.setBegin(begin)
                transcript.setEnd(end)
            else:
                transcript.setBegin(exonBegin)
                transcript.setEnd(exonEnd)
        transcript.geneId=geneId
        gene=genes.get(geneId,None)
        if(gene is None):
            genes[geneId]=gene=Gene(); gene.setId(geneId)
        transcript.setGene(gene)
        exon=Exon(exonBegin,exonEnd,transcript)
        exon.extraFields=extra
        if(not transcript.exonOverlapsExon(exon)):
            exon.frame=frame
            exon.score=exonScore
            exon.type=fields[2]
            transcript.exons.append(exon) # OK -- we sort later
        gene.addTranscript(transcript)

    def loadGFF(self,gffFilename):
        transcripts={}
        genes={}
        readOrder=Integer(1)
        GFF=open(gffFilename,"r")
        transcriptBeginEnd={}
        while(True):
            line=GFF.readline()
            if(not line): break
            if(not re.search("\S+",line)): continue
            if(re.search("^\s*\#",line)): continue
            fields=line.split("\t") ### \t added 3/24/2017
            if(len(fields)<8): raise Exception("can't parse GTF:"+line)
            if(fields[2]=="transcript" or fields[2]=="mRNA"):
                self.loadGFF_transcript(fields,line,transcriptBeginEnd,GFF,
                                   transcripts,readOrder,genes)
            elif("UTR" in fields[2] or "utr" in fields[2]):
                self.loadGFF_UTR(fields,line,transcriptBeginEnd,GFF,
                            transcripts,readOrder,genes)
            elif(fields[2]=="exon"):
                if(self.exonsAreCDS):
                    self.loadGFF_CDS(fields,line,transcriptBeginEnd,GFF,
                                     transcripts,readOrder,genes)
                else:
                    self.loadGFF_exon(fields,line,transcriptBeginEnd,GFF,
                                      transcripts,readOrder,genes)
            elif("CDS" in fields[2] or "-exon" in fields[2]):
                self.loadGFF_CDS(fields,line,transcriptBeginEnd,GFF,
                            transcripts,readOrder,genes)
        GFF.close()
        transcripts=list(transcripts.values())
        for transcript in  transcripts: transcript.parseRawExons()
        self.adjustStartCodons(transcripts)
        self.computeFrames(transcripts);
        if(self.shouldSortTranscripts):
            transcripts.sort(key=lambda t: t.substrate+" "+str(t.begin)+" "+
                             str(t.end))
        else:
            transcripts.sort(key=lambda t: t.readOrder)
        return transcripts

