#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import re

from .FastbTrack import FastbTrack

######################################################################
#
# Fastb.py bmajoros@duke.edu 10/14/2016
#
# Loads a Fastb file into memory and provides access to tracks.
#
# Attributes:
#   trackHash : hash mapping name to FastbTrack
#   trackArray : array of FastbTrack
# Methods:
#   fastb=Fastb(filename=None)
#   n=fastb.numTracks()
#   L=fastb.getLength()
#   track=fastb.getIthTrack(i)
#   track=fastb.getTrackByName(id)
#   fastb.renameTrack(oldName,newName)
#   fastb.addTrack(fastbTrack)
#   fastb.save(filename)
#   newFastb=fastb.slice(begin,end)
#   newFastb=fastb.sliceInterval(Interval)
#   fastb.dropTrack(trackName)
# Private:
#   fastb=Fastb()
#   load(filename)
######################################################################

class Fastb:
    def __init__(self,filename=None):
        self.trackArray=[]
        self.trackHash={}
        if(filename and len(filename)>0): self.load(filename)

    def numTracks(self):
        return len(self.trackArray)

    def getIthTrack(self,i):
        return self.trackArray[i]

    def getTrackByName(self,id):
        return self.trackHash[id]

    def addTrack(self,track):
        self.trackArray.append(track)
        self.trackHash[track.getID()]=track

    def save(self,filename):
        with open(filename,"w") as OUT:
            N=self.numTracks()
            for i in range(0,N): self.getIthTrack(i).save(OUT)

    def renameTrack(self,oldName,newName):
        hash=self.trackHash
        oldTrack=hash[oldName]
        if(not oldTrack): return
        oldTrack.rename(newName)
        hash[newName]=oldTrack
        del hash[oldName]

    def getLength(self):
        if(self.numTracks()==0): return 0
        return self.getIthTrack(0).getLength()

    def slice(self,begin,end):
        newFastb=Fastb()
        n=self.numTracks()
        for i in range(0,n):
            track=self.getIthTrack(i)
            newFastb.addTrack(track.slice(begin,end))
        return newFastb

    def sliceInterval(self,interval):
        return self.slice(interval.begin,interval.end)

    def dropTrack(self,name):
        n=self.numTracks()
        index=None
        for i in range(n):
            if(self.getIthTrack(i).getID()==name):
                index=i
                break
        if(index is None): raise Exception("can't find track "+name)
        del self.trackArray[index]
        del self.trackHash[name]

    def load(self,filename):
        lines=[]
        with open(filename,"r") as IN:
            while(True):
                line=IN.readline()
                if(not line): break
                line=line.rstrip()
                line=line.lstrip()
                if(re.search("\S",line)): lines.append(line)
        numLines=len(lines)
        i=0
        while(i<numLines):
            line=lines[i]
            match=re.search("^\s*([%>])\s*(\S+)(.*)",line)
            if(match):
                (op,id,rest)=(match.group(1),match.group(2),match.group(3))
                if(op==">"):
                    seq=""
                    i+=1
                    while(i<numLines):
                        line=lines[i]
                        if(re.search("^\s*([%>])\s*\S+.*",line)):
                            break
                        seq+=line
                        i+=1
                    track=FastbTrack("discrete",id,seq,rest)
                    self.addTrack(track)
                else:
                    data=[]
                    i+=1
                    while(i<numLines):
                        line=lines[i]
                        if(re.search("^\s*[%>]\s*\S+.*",line)):
                            break
                        match=re.search("(\S+)",line)
                        data.append(float(match.group(1)))
                        i+=1
                    track=FastbTrack("continuous",id,data,rest)
                    self.addTrack(track)
            else: raise Exception("can't parse line: "+line+"\n")
