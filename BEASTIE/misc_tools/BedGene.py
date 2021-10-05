#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .Interval import Interval


#=========================================================================
# Attributes:
#   ID : name of gene
#   chr : string (chromosome name)
#   strand : string ("+" or "-", or "." if unknown)
#   CDS : array of Interval
#   UTR : array of Interval
#   exons : array of Interval : includes both CDS and UTR
# Instance Methods:
#   gene=BedGene(chr,strand)
#   interval=gene.getInterval()
#   gene.addCDS(Interval(begin,end))
#   gene.addUTR(Interval(begin,end))
#   gene.addExon(Interval(begin,end))
#   gene.coalesce() # combines UTR and CDS elements into exons; sorts by coord
# Class Methods:
#=========================================================================
class BedGene:
    """BedGene"""
    def __init__(self,ID,chr,strand):
        self.ID=ID
        self.exons=[]
        self.chr=chr
        self.strand=strand
        self.CDS=[]
        self.UTR=[]
        self.exons=[]

    def getInterval(self):
        cdsInterval=self.getInterval_array(self.CDS)
        utrInterval=self.getInterval_array(self.UTR)
        exonInterval=self.getInterval_array(self.exons)
        begin=cdsInterval.begin
        end=cdsInterval.begin
        if(utrInterval.begin):
            if(not begin or utrInterval.begin<begin):
                begin=utrInterval.begin
            if(not end or utrInterval.end>end):
                end=utrInterval.end
        if(exonInterval.begin):
            if(not begin or exonInterval.begin<begin):
                begin=exonInterval.begin
            if(not end or exonInterval.end>end):
                end=exonInterval.end
        return Interval(begin,end)

    def getInterval_array(self,array):
        begin=end=None
        for interval in array:
            if(not begin or interval.begin<begin): begin=interval.begin
            if(not end or interval.end>end): end=interval.end
        return Interval(begin,end)

    def addCDS(self,interval):
        self.CDS.append(interval)

    def addUTR(self,interval):
        self.UTR.append(interval)

    def addExon(self,interval):
        self.exons.append(interval)

    def coalesce(self):
        exons=self.exons=[]
        for cds in self.CDS: exons.append(cds.clone())
        for utr in self.UTR:
            added=False
            for exon in exons:
                if(utr.begin<exon.begin and utr.end>=exon.begin):
                    exon.begin=utr.begin
                    added=True
                    break
                elif(utr.end>exon.end and utr.begin<=exon.end):
                    exon.end=utr.end
                    added=True
                    break
            if(not added): exons.append(utr)
        exons.sort(key=lambda exon: exon.begin)
