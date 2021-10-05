#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import re

from .Exon import Exon
from .Gene import Gene
from .Rex import Rex
from .Transcript import Transcript


#=========================================================================
# Attributes:
#
# Instance Methods:
#   reader=GFF3Parser() =TESTED
#   transcriptArray=reader.loadGFF(filename) =TESTED
#   geneList=reader.loadGenes(filename) =TESTED
#   hashTable=reader.hashBySubstrate(filename) =TESTED
#   hashTable=reader.hashGenesBySubstrate(filename) =TESTED
#   hashTable=reader.loadTranscriptIdHash(filename) =TESTED
#   hashTable=reader.loadGeneIdHash(filename) =TESTED
#
# Private Methods:
#   structure=parser.loadStructure(filename)
#   records=parser.loadRecord(filename)
#   record=self.parseRecord(fields)
#=========================================================================
class GFF3Parser:
    """GFF3Parser"""
    def __init__(self):
        pass

    def loadTranscriptIdHash(self,filename):
        transcripts=self.loadGFF(filename)
        hash={}
        for transcript in transcripts:
            hash[transcript.getID()]=transcript
        return hash

    def loadGeneIdHash(self,filename):
        genes=self.loadGenes(filename)
        hash={}
        for gene in genes:
            hash[gene.getID()]=gene
        return hash

    def hashGenesBySubstrate(self,filename):
        genes=self.loadGenes(filename)
        hash={}
        for gene in genes:
            substrate=gene.getSubstrate()
            array=hash.get(substrate,None)
            if(array is None): array=hash[substrate]=[]
            array.append(gene)
        return hash

    def hashBySubstrate(self,filename):
        transcripts=self.loadGFF(filename)
        hash={}
        for transcript in transcripts:
            substrate=transcript.getSubstrate()
            array=hash.get(substrate,None)
            if(array is None): array=hash[substrate]=[]
            array.append(transcript)
        return hash

    def loadGFF(self,filename):
        genes=self.loadGenes(filename)
        transcripts=[]
        for gene in genes:
            n=gene.getNumTranscripts()
            for i in range(n):
                transcript=gene.getIthTranscript(i)
                transcripts.append(transcript)
        return transcripts

    def makeGene(self,root):
        gene=Gene()
        root["object"]=gene
        gene.ID=root["extra"]["ID"]
        children=root.get("children",None)
        if(children is None): return gene
        for child in children:
            obj=self.labelStructure(child)
            if(obj is None): continue
            if(type(obj)==Transcript):
                gene.addTranscript(obj)
                obj.gene=gene
                obj.geneId=gene.getId()
        extra=root["extra"]
        gene.extraFields=""
        for key in extra:
            gene.extraFields+=key+"="+extra[key]+";"
        return gene

    def makeTranscript(self,root):
        id=root["extra"]["ID"]
        strand=root["strand"]
        transcript=Transcript(id,strand)
        root["object"]=transcript
        transcript.substrate=root["substrate"]
        transcript.source=root["source"]
        transcript.begin=int(root["begin"])
        transcript.end=int(root["end"])
        children=root.get("children",None)
        if(children is None): return transcript
        for child in children:
            obj=self.labelStructure(child)
            if(obj is None): continue
            if(type(obj)==Exon):
                childType=child["type"]
                if(childType=="CDS" or
                   re.search("-exon",childType)): transcript.addExon(obj)
                elif(childType=="exon"): transcript.addRawExon(obj)
                elif(re.search("UTR",childType)): transcript.addUTR(obj)
                obj.transcript=transcript
        transcript.parseRawExons()
        transcript.setExonTypes()
        transcript.setUTRtypes()
        transcript.sortExons()
        transcript.adjustOrders()
        extra=root["extra"]
        transcript.extraFields=""
        for key in extra:
            transcript.extraFields+=key+"="+extra[key]+";"
        return transcript

    def makeExon(self,root):
        begin=int(root["begin"])
        end=int(root["end"])
        exon=Exon(begin,end,None)
        exon.strand=root["strand"]
        exon.frame=root["frame"]
        exon.type=root["type"]
        exon.score=root["score"]
        exon.substrate=root["substrate"]
        extra=root["extra"]
        exon.extraFields=""
        for key in extra:
            exon.extraFields+=key+"="+extra[key]+";"
        return exon

    def labelStructure(self,root):
        obj=None
        t=root["type"]
        if(t=="gene"):
            obj=self.makeGene(root)
        elif(t=="transcript" or t=="mRNA"):
            obj=self.makeTranscript(root)
        elif(t=="exon" or t=="CDS"):
            obj=self.makeExon(root)
        return obj

    def loadGenes(self,filename):
        roots=self.loadStructure(filename)
        genes=[]
        for root in roots:
            obj=self.labelStructure(root)
            if(type(obj)==Gene): genes.append(obj)
        return genes

    def loadStructure(self,filename):
        records=self.loadRecords(filename)
        idHash={}
        self.hashRecordsByID(records,idHash)
        self.connectParentsChildren(records,idHash)
        roots=self.findRoots(records)
        #for root in roots: self.printStructure(root)
        return roots

    def printStructure(self,rec,depth=0):
        print("\t"*depth+rec["type"]+" "+rec["extra"]["ID"])
        children=rec.get("children",None)
        if(not children): return
        for child in children:
            self.printStructure(child,depth+1)

    def findRoots(self,records):
        roots=[]
        for record in records:
            if(record.get("parent",None) is None):
                roots.append(record)
        return roots

    def addChild(self,parent,child):
        if(parent.get("children",None) is None): parent["children"]=[]
        parent["children"].append(child)
        child["parent"]=parent

    def connectParentsChildren(self,records,idHash):
        for record in records:
            parent=record["extra"].get("Parent",None)
            if(parent is not None):
                parents=parent.split(",")
                for parent in parents:
                    parentRec=idHash.get(parent,None)
                    if(parentRec is None):
                        raise Exception("Cannot find GFF3 parent "+parent)
                    self.addChild(parentRec,record)

    def hashRecordsByID(self,records,hash):
        for record in records:
            extraHash=record["extra"]
            ID=extraHash.get("ID",None)
            if(ID is not None): hash[ID]=record

    def loadRecords(self,filename):
        fh=open(filename,"rt")
        records=[]
        while(True):
            line=fh.readline()
            if(line==""): break
            fields=line.split(sep="\t")
            numFields=len(fields)
            if(numFields<9): continue
            rec=self.parseRecord(fields)
            records.append(rec)
        fh.close()
        return records

    def parseRecord(self,fields):
        if(len(fields)>9):
            raise Exception("too many fields in GFF3 record"+"\t".join(fields))
        (substrate,source,type,begin,end,score,strand,frame,extra)=fields
        extra=extra.rstrip()
        extraFields=extra.split(";")
        extraHash={}
        rex=Rex()
        for field in extraFields:
            if(not rex.find("(.+)=(.+)",field)):
                raise Exception("Can't parse GFF3 field: "+field)
            key=rex[1]; value=rex[2]
            extraHash[key]=value
        rec={"substrate":substrate,
             "source":source,
             "type":type,
             "begin":int(begin)-1,
             "end":int(end),
             "score":score,
             "strand":strand,
             "frame":frame,
             "extra":extraHash}
        return rec



# =========================================================================
def test_parser1(filename):
    parser=GFF3Parser()
    genes=parser.loadGenes(filename)
    for gene in genes:
        id=gene.getId()
        numTrans=gene.getNumTranscripts()
        print(gene.getId(),gene.getStrand(),numTrans)
        for i in range(numTrans):
            trans=gene.getIthTranscript(i)
            trans.shiftCoords(1000)
            numExons=trans.numExons()
            print("\t",trans.getTranscriptId(),trans.getStrand(),numExons,
                  trans.getSubstrate(),trans.getSource(),trans.getGeneId(),
                  trans.getBegin(),trans.getEnd(),trans.getScore())
            for i in range(numExons):
                exon=trans.getIthExon(i)
                print("\t\tEXON ",exon.getBegin(),exon.getEnd(),
                      exon.getStrand(),
                      exon.order,exon.getTranscript().getID(),exon.getFrame(),
                      exon.getScore(),exon.getSubstrate(),exon.getType())
            numUTR=trans.numUTR()
            for i in range(numUTR):
                exon=trans.getIthUTR(i)
                print("\t\tUTR ",exon.getBegin(),exon.getEnd(),exon.getStrand(),
                      exon.order,exon.getTranscript().getID(),exon.getFrame(),
                      exon.getScore(),exon.getSubstrate(),exon.getType())

def test_parser2(filename):
    parser=GFF3Parser()
    transcripts=parser.loadGFF(filename)
    numTrans=len(transcripts)
    for i in range(numTrans):
        trans=transcripts[i]
        id=trans.getID()
        numExons=trans.numExons()
        print("\t",trans.getId(),trans.getStrand(),numExons)
        for i in range(numExons):
            exon=trans.getIthExon(i)
            print("\t\tEXON ",exon.getBegin(),exon.getEnd(),
                  exon.getStrand(),
                  exon.order,exon.getTranscript().getID(),exon.getFrame(),
                  exon.getScore(),exon.getSubstrate(),exon.getType())
        numUTR=trans.numUTR()
        for i in range(numUTR):
            exon=trans.getIthUTR(i)
            print("\t\tUTR ",exon.getBegin(),exon.getEnd(),exon.getStrand(),
                  exon.order,exon.getTranscript().getID(),exon.getFrame(),
                  exon.getScore(),exon.getSubstrate(),exon.getType())

def test_parser3(filename):
    reader=GFF3Parser()
    hashTable=reader.hashBySubstrate(filename)
    for substrate in hashTable:
        transcripts=hashTable[substrate]
        for transcript in transcripts:
            print(transcript.getID())

def test_parser4(filename):
    reader=GFF3Parser()
    hashTable=reader.hashGenesBySubstrate(filename)
    for substrate in hashTable:
        genes=hashTable[substrate]
        for gene in genes:
            print(gene.getID())

def test_parser5(filename):
    reader=GFF3Parser()
    hash=reader.loadTranscriptIdHash(filename)
    for id in hash:
        transcript=hash[id]
        print(transcript.getId())

def test_parser6(filename):
    reader=GFF3Parser()
    hash=reader.loadGeneIdHash(filename)
    for id in hash:
        gene=hash[id]
        print(gene.getId())

def test_parser7(filename):
    reader=GFF3Parser()
    genes=reader.loadGenes(filename)
    for gene in genes:
        print(gene.toGff())

#test_parser7("/Users/bmajoros/python/test/data/subset.gff3")
