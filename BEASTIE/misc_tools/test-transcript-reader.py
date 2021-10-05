#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .GffTranscriptReader import GffTranscriptReader

#filename="/home/bmajoros/1000G/assembly/local-genes.gff"
#filename="/home/bmajoros/1000G/assembly/tmp.gff"
#filename="test/data/tmp.gff"
#filename="test/data/local-genes.gff"
filename="/home/bmajoros/ensembl/protein-coding.gff"

reader=GffTranscriptReader()
genes=reader.loadGenes(filename)
for gene in genes:
    exons=gene.getMergedExons()
    unmerged=0
    for transcript in gene.transcripts:
        unmerged+=len(transcript.getRawExons())
    print(unmerged,"exons merged to",len(exons))
    #for i in range(len(exons)):
    #    print("MERGED TO:",exons[i].begin,exons[i].end)
    #    print()

#transcripts=reader.loadGFF(filename)
#for transcript in transcripts:
    #print(transcript.getID())
    #gff=transcript.toGff()
    #print(gff)

#genes=reader.loadGenes(filename)
#for gene in genes:
#    print("gene",gene.getID())
#    n=gene.getNumTranscripts()
#    for i in range(n):
#        transcript=gene.getIthTranscript(i)
#        transID=transcript.getID()
#        print("\t"+transID+"\t"+str(transcript.getBegin())+"\t"
#              +str(transcript.getEnd()))

#hashTable=reader.hashBySubstrate(filename)
#keys=hashTable.keys()
#for key in keys:
#    print(key)

#hashTable=reader.hashGenesBySubstrate(filename)
#keys=hashTable.keys()
#for key in keys:
#    print(key)

#hashTable=reader.loadTranscriptIdHash(filename)
#keys=hashTable.keys()
#for key in keys:
#    print(key)

#hashTable=reader.loadGeneIdHash(filename)
#keys=hashTable.keys()
#for key in keys:
#    print(key)

