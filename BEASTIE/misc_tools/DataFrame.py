#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# 2018 William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from .DataFrameRow import DataFrameRow
from .Rex import Rex

rex=Rex()

#=========================================================================
# Attributes:
#   header
#   matrix : array of rows, each of which is a DataFrameRow
#   rowHash : dictionary mapping row names to row indices
#   colHash : dictionary mapping column names to column indices
# Methods:
#   df=DataFrame()
#   df.save(filename)
#   rowNames=df.getRowNames()
#   colNames=df.getColumnNames()
#   df.addRow(DataFrameRow)
#   n=df.nrow()
#   n=df.ncol()
#   row=df[index]
#   rows=df.getRows()
#   elem=df[i][j]
#   df.toInt()
#   df.toFloat()
#   df.colToFloat(colIndex)
#   header=df.getHeader()
#   df.removeQuotes()
#   df.hashRowNames()
#   df.hashColNames()
#   row=df.getRowI(i)
#   col=df.getColI(i)
#   row=df.getRow(rowName) # call hashRowNames() first!
#   col=df.getColumn(columnName) # call hashColNames() first!
#   bool=df.rowExists(rowName) # call hashRowNames() first!
#   bool=df.columnExists(colName) # call hashColNames() first!
#   index=df.getColumnIndex(colName) # call hashColNames() first!
#   newDataFrame=df.subsetColumns(colIndices)
#   newDataFrame=df.subsetRows(rowIndices)
#   idx=df.addColumn(colName,defaultValue) # returns index of new column
#   df.print(handle)
#   array=df.toDataArray()
#   df.appendDF(otherDF) # does NOT do a deep copy!
# Class methods:
#   df=DataFrame.readTable(filename,header=False,rowNames=False)
#=========================================================================

class DataFrame:
   def __init__(self):
      self.header=[]
      self.matrix=[]
      self.rowHash=None
      self.colHash=None

   def save(self,filename):
      with open(filename,"wt") as OUT:
         self.print(OUT)

   def appendDF(self,other):
      self.matrix.extend(other.matrix)

   def addRow(self,row):
      self.matrix.append(row)

   def getRows(self):
      return self.matrix

   def removeQuotes(self):
      for row in self.matrix:
         raw=row.getRaw()
         for i in range(len(raw)):
            if(rex.find("\"\s*(\S+)\"",raw[i])):
               raw[i]=rex[1]
      self.unquoteHeader()

   def unquoteHeader(self):
      raw=self.header
      for i in range(len(raw)):
         if(rex.find("\"\s*(\S+)\"",raw[i])):
            raw[i]=rex[1]

   def toDataArray(self):
      array=[]
      for row in self.matrix:
         array.append(row.values)
      return array

   def print(self,handle):
      print("\t".join(self.header),file=handle)
      for row in self.matrix: row.print(handle)

   def addColumn(self,name,defaultValue):
      colIndex=len(self.header)
      self.header.append(name)
      for row in self.matrix:
         row.append(defaultValue)
      return colIndex

   def subsetColumns(self,colIndices):
      newDF=DataFrame()
      header=self.header
      newHeader=newDF.header
      for i in colIndices: newHeader.append(header[i])
      for i in range(self.nrow()):
         row=self[i]
         newRow=DataFrameRow()
         newRow.rename(row.getLabel())
         for j in colIndices: newRow.values.append(row[j])
         newDF.matrix.append(newRow)
      return newDF

   def subsetRows(self,rowIndices):
      newDF=DataFrame()
      newDF.header=self.header
      for i in rowIndices:
         newDF.addRow(self[i].clone())
      return newDF

   def rowExists(self,rowName):
      if(self.rowHash is None): raise Exception("call hashRowNames() first")
      return self.rowHash.get(rowName,None) is not None

   def getColumnIndex(self,colName):
      return self.colHash.get(colName)

   def columnExists(self,colName):
      if(self.colHash is None): raise Exception("call hashColNames() first")
      return self.colHash.get(colName,None) is not None

   def getRowNames(self):
      names=[]
      for row in self.matrix:
         names.append(row.label)
      return names

   def getColumnNames(self):
      return self.header

   def getRowI(self,rowIndex):
      return self.matrix[rowIndex]

   def getColI(self,colIndex):
      column=DataFrameRow()
      for row in self.matrix:
         column.values.append(row[colIndex])
      return column

   def getRow(self,rowName):
      if(self.rowHash is None): raise Exception("call hashRowNames() first")
      rowIndex=self.rowHash.get(rowName,None)
      if(rowIndex is None): raise Exception("row not found: "+rowName)
      return self.matrix[rowIndex]

   def getColumn(self,colName):
      if(self.colHash is None): raise Exception("call hashColNames() first")
      colIndex=self.colHash.get(colName,None)
      if(colIndex is None): raise Exception("column not found: "+colName)
      column=DataFrameRow()
      column.label=colName
      for row in self.matrix:
         column.values.append(row[colIndex])
      return column

   def hashRowNames(self):
      h=self.rowHash={}
      numRows=self.nrow()
      for i in range(numRows):
         row=self.matrix[i]
         h[row.label]=i

   def hashColNames(self):
      h=self.colHash={}
      numCols=self.ncol()
      for i in range(numCols):
         h[self.header[i]]=i

   def getHeader(self):
      return self.header

   def nrow(self):
      return len(self.matrix)

   def ncol(self):
      if(len(self.header)!=0): return len(self.header)
      if(len(self.matrix)==0): return 0
      return self.matrix[0].length()

   def __getitem__(self,i):
      return self.matrix[i]

   def toInt(self):
      for row in self.matrix: row.toInt()

   def colToFloat(self,colIndex):
      for row in self.matrix: row[colIndex]=float(row[colIndex])

   def toFloat(self):
      for row in self.matrix: row.toFloat()

   @classmethod
   def readTable(cls,filename,header=False,rowNames=False):
      df=DataFrame()
      with open(filename,"rt") as IN:
         if(header):
            df.header=IN.readline()
            df.header=df.header.rstrip().split() #("\t")
         for line in IN:
            fields=line.rstrip().split() #("\t")
            if(len(fields)<1): continue
            label=""
            if(rowNames):
               label=fields[0]
               fields=fields[1:]
            row=DataFrameRow()
            row.label=label
            row.values=fields
            df.matrix.append(row)
      if(len(df.matrix)>0 and df.matrix[0].length()<len(df.header)):
         df.header=df.header[1:]
      return df

