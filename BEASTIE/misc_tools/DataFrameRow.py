#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# 2018 William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================

#=========================================================================
# Attributes:
#   label : string
#   values : array of values
# Methods:
#   row=DataFrameRow()
#   elem=row[i] # first element is at 0 (the label is not counted)
#   label=row.getLabel()
#   raw=raw.getRaw()
#   row.rename(label)
#   n=row.length()
#   row.toInt()
#   row.toFloat()
#   row.append(value)
#   row.print(handle)
#   newRow=row.clone()
#   row.log()
#   row.log2()
#   row.log10()
#=========================================================================

class DataFrameRow:
   def __init__(self):
      self.label=""
      self.values=[]

   def log(self):
      values=self.values
      for i in range(len(values)):
         values[i]=log(values[i])

   def log2(self):
      values=self.values
      for i in range(len(values)):
         values[i]=log2(values[i])

   def log10(self):
      values=self.values
      for i in range(len(values)):
         values[i]=log10(values[i])

   def getRaw(self):
      return self.values

   def clone(self):
      r=DataFrameRow()
      r.label=self.label
      for x in self.values:
         r.values.append(x)
      return r

   def __getitem__(self,i):
      return self.values[i]

   def __setitem__(self,i,value):
      self.values[i]=value

   def print(self,handle):
      if(self.label!=""): print(self.label+"\t",end="",file=handle)
      print("\t".join([str(x) for x in self.values]),sep="",file=handle)

   def append(self,value):
      self.values.append(value)

   def length(self):
      return len(self.values)

   def getLabel(self):
      return self.label

   def rename(self,x):
      self.label=x

   def toInt(self):
      self.values=[int(x) for x in self.values]

   def toFloat(self):
      self.values=[float(x) for x in self.values]

