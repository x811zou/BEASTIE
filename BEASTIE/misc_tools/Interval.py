#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import sys

from .Rex import Rex

rex=Rex()

#=========================================================================
# Attributes:
#   begin
#   end
# Methods:
#   i=Interval(begin,end)
#   s=interval.toString()
#   print(file=STDOUT)
#   bool=interval.overlaps(other)
#   bool=interval.contains(position)
#   bool=interval.containsInterval(other)
#   distance=interval.distance(other)
#   distance=interval.distanceFromPoint(x)
#   intersection=interval.intersect(other)
#   union=interval.union(other) # returns an array of intervals
#   diff=interval.minus(other)  # returns an array of intervals
#   length=interval.getLength()
#   begin=interval.getBegin()
#   end=interval.getEnd()
#   length=interval.length()
#   bool=interval.equals(other)
#   other=interval.clone()
#   bool=interval.isEmpty()
#   d=interval.relativeDistanceFromBegin(pos)
#   d=interval.relativeDistanceFromEnd(pos)
#   interval.shift(delta)
#   center=interval.floatCenter()
#   center=interval.intCenter()
#   center=interval.center() # same as floatCenter()
# Class methods:
#   interval=Interval.parseInt("(10,15)")
#   interval=Interval.parseFloat("(3.5,7.2)")
#=========================================================================

class Interval:
   """"Interval represents an interval [a,b) in which a is inclusive and
       b is not
   """
   def __init__(self,begin=0,end=0):
      self.begin=begin
      self.end=end

   @classmethod
   def parseInt(cls,interval):
      if(rex.find("\(([^,]+),([^\)]+)\)",interval)):
         return Interval(int(rex[1]),int(rex[2]))
      if(rex.find("([^,]+):([^\)]+)",interval)):
         return Interval(int(rex[1]),int(rex[2]))
      if(rex.find("([^,]+)-([^\)]+)",interval)):
         return Interval(int(rex[1]),int(rex[2]))
      return None

   @classmethod
   def parseFloat(cls,interval):
      if(rex.find("\(([^,]+),([^\)]+)\)",interval)):
         return Interval(float(rex[1]),float(rex[2]))
      return None

   def print(self,file=sys.stdout):
      print("(",self.begin,",",self.end,")",sep="",end="",file=file)

   def toString(self):
      return "("+str(self.begin)+","+str(self.end)+")"

   def overlaps(self,other):
      return self.begin<other.end and other.begin<self.end

   def distanceFromPoint(self,x):
      if(self.contains(x)): return 0
      d1=abs(self.begin-x)
      d2=abs(self.end-x)
      return d1 if d1<d2 else d2

   def distance(self,other):
      if(self.overlaps(other)): return 0
      d=self.begin-other.end
      if(d>0): return d
      return other.begin-self.end

   def contains(self,index):
      return index>=self.begin and index<self.end

   def containsInterval(self,other):
      return self.begin<=other.begin and self.end>=other.end

   def isEmpty(self):
      return self.begin>=self.end

   def clone(self):
      return Interval(self.begin,self.end)

   def intersect(self,other):
      if(self.isEmpty()): return self.clone()
      if(other.isEmpty()): return other.clone()
      begin=max(self.begin,other.begin)
      end=min(self.end,other.end)
      return Interval(begin,end)

   def length(self):
      L=self.end-self.begin
      return 0 if L<0 else L

   def getLength(self):
      return self.length()

   def getBegin(self):
      return self.begin

   def getEnd(self):
      return self.end

   def equals(self,other):
      return self.begin==other.begin and self.end==other.end

   def relativeDistanceFromBegin(self,pos):
      if(not self.contains(pos)): raise IndexError(pos)
      return (pos-self.begin)/self.length()

   def relativeDistanceFromEnd(self,pos):
      if(not self.contains(pos)): raise IndexError(pos)
      return (self.end-pos)/self.length()

   def union(self,other):
      s=[]
      if(self.overlaps(other)):
         s.append(Interval(self.begin,other.end))
      else:
         s.append(self,other)
      return s

   def minus(self,other):
      if(not self.overlaps(other)): return [self]
      s=[]
      if(self.begin<other.begin):
         s.append(Interval(self.begin,other.begin))
      if(self.end>other.end):
         s.append(Interval(other.end,self.end))
      return s

   def __str__(self):
      return "("+str(self.begin)+","+str(self.end)+")"

   def __repr__(self):
      return "("+str(self.begin)+","+str(self.end)+")"

   def shift(self,delta):
      self.begin+=delta
      self.end+=delta

   def floatCenter(self):
      return float(self.begin+self.end)/2.0

   def intCenter(self):
      return int((self.begin+self.end)/2)

   def center(self):
      return self.floatCenter()


