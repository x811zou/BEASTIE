#!/usr/bin/env python
import gzip
import sys

if(len(sys.argv)!=4):
   print(sys.argv[0]+" <in.vcf.gz> <begin> <end>")
   sys.exit(0)
[infile,begin,end]=sys.argv[1:]
begin=int(begin)
end=int(end)

f=gzip.open(infile)
for line in f:
   line.rstrip("\n")
   fields=line.split()
   if len(fields)<7: continue
   if fields[0]=="#CHROM":
      header=fields
      print(line)
   elif fields[6]=="PASS":
      [chr,pos,id,ref,alt]=fields[:5]
      pos=int(pos)
      if pos>=begin and pos<=end: print(line)


