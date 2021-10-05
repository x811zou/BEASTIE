#!/usr/bin/env python
import operator
import sys

if(len(sys.argv)!=2):
   print sys.argv[0]+" <in.txt>"
   sys.exit(0)
[infile]=sys.argv[1:]

f=open(infile)
header=f.readline()
variants=header.split()
counts={}
for line in f:
   line.rstrip("\n")
   fields=line.split()
   n=len(fields)
   for i in range(1,n):
      id=variants[i-1]
      if fields[i]!="1|1": continue
      if not id in counts: counts[id]=1
      else: counts[id]+=1
ranked=sorted(counts.items(), key=operator.itemgetter(1),
       reverse=True)
for variant in ranked: print variant[0],"\t",variant[1]
