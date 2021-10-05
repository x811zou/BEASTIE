#!/usr/bin/env python
import os
import sys

name=sys.argv[0];

if(len(sys.argv)!=1):
   print name+" <dna.sam>";
   sys.exit(0)
[name,dnaRep]=sys.argv;

alignmentDir="aligned";

def System(cmd):
   print "Excecuting: "+cmd
   os.system(cmd)

System("get-haplo-read-counts.pl "+alignmentDir+" "+dnaRep+" > haplo-read-counts.txt");

System("get-SNP-allele-counts.pl > SNP-read-counts.txt");

System("filter-SNP-read-counts.pl SNP-read-counts.txt > SNP-read-counts-filtered.txt");

System("SNP-fisher.R SNP-read-counts-filtered.txt > SNP-fdr.txt");

System("analyze-haplotypes.pl "+alignmentDir+" "+dnaRep+" > analyze-haplotypes.txt");

System("filter-zeros.pl > nonzero.txt");

System("fdr.R > haplotype-fdr.txt");

System("join2.pl > join.txt");



