#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .FastaWriter import FastaWriter
from .Interval import Interval

######################################################################
#
# FastbTrack.py bmajoros@duke.edu 10/14/2016
#
# A track in a Fastb file.
#
# Attributes:
#   type : "discrete" or "continuous"
#   id : name of track
#   data : string (for discrete) or array of float (for continuous)
#   deflineExtra : extra info for defline
# Methods:
#   track=FastbTrack(type,id,data,deflineExtra="") # type="discrete" or "continuous"
#   type=track.getType()
#   data=track.getData()
#   track.setSequence(string) # discrete
#   track.setData(values) # continuous
#   id=track.getID()
#   L=track.getLength()
#   track.rename(newID)
#   bool=track.isDiscrete()
#   bool=track.isContinuous()
#   track.save(FILEHANDLE)
#
# Methods for discrete tracks:
#   getDiscreteContiguousRegions(self) # returns array of Interval with
#      "value" attribute added
#
# Methods for continuous tracks:
#   array=track.getNonzeroRegions() # returns array of Interval
#   array=track.getZeroRegions() # returns array of Interval
#   array=track.getRegionsAbove(cutoff) # returns array of Interval
#   bool=track.anyZeroValues() # only for continuous tracks
#   array=track.getContiguousRegions() # returns an array of Interval with
#      "value" attribute added
#   newTrack=track.slice(begin,end) # [begin,end) => end not inclusive
#   meanValue=track.getMean(interval=None) # only for continuous data
#   (maxValue,maxPos)=track.getMax(interval=None) # only for continuous data
######################################################################

class FastbTrack:
    def __init__(self,type,id,data,deflineExtra=""):
        self.type=type
        self.id=id
        self.data=data
        self.deflineExtra=deflineExtra

    def getType(self):
        return self.type

    def getData(self):
        return self.data

    def isDiscrete(self):
        return self.type=="discrete"

    def isContinuous(self):
        return self.type=="continuous"

    def getID(self):
        return self.id

    def save(self,fh):
        writer=FastaWriter()
        id=self.id
        data=self.data
        deflineExtra=self.deflineExtra
        if(self.isDiscrete()):
            writer.addToFasta(">"+id+" "+deflineExtra,data,fh)
        else:
            fh.write("%"+id+" "+deflineExtra+"\n")
            n=len(data)
            for i in range(0,n): fh.write(str(data[i])+"\n")

    def getContiguousRegions(self):
        """getContiguousRegions() returns an array of Interval objects
        with a "value" attribute added
        """
        data=self.data
        if(self.isDiscrete()): raise Exception("track is not continuous")
        L=len(data)
        intervals=[]
        if(L==0): return intervals
        begin=0
        for i in range(1,L):
            x=data[i]
            prev=data[i-1]
            if(x==prev): continue
            interval=Interval(begin,i)
            interval.value=prev
            intervals.append(interval)
            begin=i
        interval=Interval(begin,L)
        interval.value=data[L-1]
        intervals.append(interval)
        return intervals

    def getDiscreteContiguousRegions(self):
        """getDiscreteContiguousRegions() returns an array of Interval objects
        with a "value" attribute added
        """
        data=self.data
        if(not self.isDiscrete()): raise Exception("track is not discrete")
        L=len(data)
        intervals=[]
        if(L==0): return intervals
        begin=0
        for i in range(1,L):
            x=data[i]
            prev=data[i-1]
            if(x==prev): continue
            interval=Interval(begin,i)
            interval.value=prev
            intervals.append(interval)
            begin=i
        interval=Interval(begin,L)
        interval.value=data[L-1]
        intervals.append(interval)
        return intervals

    def getZeroRegions(self):
        """getNonzeroRegions() returns array of Intervals"""
        data=self.data
        if(self.isDiscrete()): raise Exception("track is not continuous")
        L=len(data)
        intervals=[]
        begin=None
        for i in range(L):
            x=data[i]
            if(x==0 and (i==0 or data[i-1]!=0)): begin=i
            elif(x!=0 and i>0 and data[i-1]==0):
                intervals.append(Interval(begin,i))
        if(L>0 and data[L-1]==0):
            intervals.append(Interval(begin,L))
        return intervals

    def getNonzeroRegions(self):
        """getNoneroRegions() returns array of Intervals"""
        data=self.data
        if(self.isDiscrete()): raise Exception("track is not continuous")
        L=len(data)
        intervals=[]
        begin=None
        for i in range(0,L):
            x=data[i]
            if(x!=0 and (i==0 or data[i-1]==0)): begin=i
            elif(x==0 and i>0 and data[i-1]!=0):
                intervals.append(Interval(begin,i))
        if(L>0 and data[L-1]!=0):
            intervals.append(Interval(begin,L))
        return intervals

    def getRegionsAbove(self,cutoff):
        """getRegionsAbove() returns array of Intervals"""
        data=self.data
        if(self.isDiscrete()): raise Exception("track is not continuous")
        L=len(data)
        intervals=[]
        begin=None
        for i in range(0,L):
            x=data[i]
            if(x>cutoff and (i==0 or data[i-1]<=cutoff)): begin=i
            elif(x<=cutoff and i>0 and data[i-1]>cutoff):
                intervals.append(Interval(begin,i))
        if(L>0 and data[L-1]>cutoff):
            intervals.append(Interval(begin,L))
        return intervals

    def rename(self,newID):
        self.id=newID

    def getLength(self):
        return len(self.data)

    def slice(self,begin,end):
        return FastbTrack(self.type,self.id,self.data[begin:end],
                          self.deflineExtra)

    def setSequence(self,string):
        self.data=string

    def setData(self,values):
        self.data=values

    def getMean(self,interval=None):
        data=self.data
        L=len(data)
        begin=0
        end=L
        if(interval is not None):
            begin=interval.begin
            end=interval.end
            L=end-begin
        sum=0.0
        for i in range(begin,end): sum+=data[i]
        mean=sum/L
        return mean

    def getMax(self,interval=None):
        data=self.data
        L=len(data)
        begin=0
        end=L
        if(interval is not None):
            begin=interval.begin
            end=interval.end
            L=end-begin
        max=data[begin]; maxPos=begin
        for i in range(begin+1,end):
            if(data[i]>max):
                max=data[i]
                maxPos=i
        return (max,maxPos)

    def anyZeroValues(self):
        data=self.data
        for x in data:
            if(x==0.0): return True
        return False

