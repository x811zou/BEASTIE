#!/usr/bin/env python
import re
import sys

name=sys.argv[0];
if(len(sys.argv)!=2):
   print(name+" <*.vcf>")
   sys.exit(0)
[name,vcfFile]=sys.argv;

IN=open(vcfFile,"r")
while(True):
  line=IN.readline()
  #if(line is ""):
  if not line: break
  line.rstrip("\n")
  if(re.search("#",line)): continue
  fields=line.split()
  if(len(fields)<10): continue
  [chr,pos,id,ref,alt]=fields[:5]
  if(len(ref)!=1 or len(alt)!=1): continue
  if(not re.search("rs",id)):
     id=chr+"at"+pos;
  print(id+"\t"+chr+"\t"+pos+"\t"+ref+"\t"+alt)
IN.close();

