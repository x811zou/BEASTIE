#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .FastaReader import FastaReader
from .FastaWriter import FastaWriter
from .Translation import Translation

#filename="/home/bmajoros/1000G/assembly/combined/HG00096/1.fasta"

#reader=FastaReader("/home/bmajoros/1000G/assembly/combined/HG00096/1.fasta")
#while(True):
#    [defline,seq]=reader.nextSequence()
#    if(not defline): break
#    print("defline="+defline);
#    L=len(seq)
#    print("length="+str(L))

#filename="/home/bmajoros/1000G/assembly/BRCA1-NA19782.fasta";
filename="/Users/bmajoros/python/test/data/subset.fasta"
print(FastaReader.getSize(filename))

[defline,seq]=FastaReader.firstSequence(filename)
print(len(seq))

#filename="/home/bmajoros/1000G/assembly/test.fasta"
filename="/Users/bmajoros/python/test/data/subset.fasta"
hash=FastaReader.readAllAndKeepDefs(filename)
for key in hash.keys():
    [defline,seq]=hash[key]
    print(defline)
    [id,attrs]=FastaReader.parseDefline(defline)
    print("id="+id)
    for key,value in attrs.items():
        print(key+"="+value)

writer=FastaWriter()
writer.writeFasta(">ABCD","ATCGATCGTAGCTAGTCTGCGCGTATCGTCAGTCTCTATCGATCGTACTGCGATCTAGCTAGCTGATCGTAGCTTCTATGACTGCTAGTCATCTAGCTAGCTGATCGTAGCTGCGCGCGATATATTGCATCTATGCTATCATTGCATGCTAGCTCTAGCTAGTCGATGCTATCTTAGCTAC","test1.fasta")

writer.appendToFasta(">XYZ","GATTACA","test1.fasta")

print(Translation.translate(seq))
print("forward:",seq)
print("revcomp: ",Translation.reverseComplement(seq))




