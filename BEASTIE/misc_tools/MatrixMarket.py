#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
import gzip

from .Rex import Rex

rex=Rex()

#=========================================================================
# Attributes:
#    FH : file handle
#    header : array of int
#    nextLine : string
# Instance Methods:
#    MatrixMarket(filename)
#    nextGroup(self,colIndex)
#    getHeader()
# Class Methods:
#    allGroups=loadFile(filename,colIndex)
#=========================================================================
class MatrixMarket:
    def __init__(self,filename):
        if(rex.find("\.gz$",filename)):
            self.FH=gzip.open(filename,"rt")
        else:
            self.FH=open(filename,"rt")
        self.header=None
        self.nextLine=None

    def getHeader(self):
        return self.header

    def nextGroup(self,colIndex):
        line=None
        # First, see if the header needs to be parsed
        while(True):
            line=self.nextLine
            if(line is None): line=self.FH.readline()
            if(line is None): return None
            L=len(line)
            if(L>0 and line[0]=="%"): continue
            break
        # The first non-comment line contains the totals
        if(self.header is None):
            self.header=line.rstrip().split()
            self.header=[int(x) for x in self.header]
            line=self.FH.readline()
        # Now we can read in the next group of lines
        prevID=None
        group=[]
        if(self.nextLine is not None): # buffered from previous call
            prevID=int(self.nextLine.rstrip().split()[colIndex])
        while(True):
            fields=line.rstrip().split()
            if(len(fields)==0): return None
            if(prevID is None): prevID=int(fields[colIndex])
            if(colIndex>len(fields)-1):
                raise Exception("colIndex=",colIndex,"len(fields)=",len(fields))
            thisID=int(fields[colIndex])
            if(thisID!=prevID):
                self.nextLine=line
                return group
            group.append(fields)
            line=self.FH.readline()
            if(line is None): return None

    @classmethod
    def loadFile(self,filename,colIndex):
        reader=MatrixMarket(filename)
        groups=[]
        while(True):
            group=reader.nextGroup(colIndex)
            if(group is None): break
            groups.append(group)
        return groups
