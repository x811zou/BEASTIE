#!/usr/bin/env python
import distutils
import os
import sys

name=sys.argv[0];

if(len(sys.argv)!=5):
   print "\n"+name+""" <SNPs.vcf> <chrom.fasta> <amplicons.bed> <fastq-dir>
      SNPs.vcf      = VCF file containing variants
      chrom.fasta   = FASTA file containing a single chromosome DNA sequence
      amplicons.bed = BED file containing chromosome coordinates of amplicons
      fastq-dir     = path to directory containing FASTQ files
   """
   sys.exit(0)
[name,vcfFile,chrFile,ampliconsBed,fastqDir]=sys.argv
workingDir=os.getcwd()
alignmentsDir="aligned"
if(not os.path.exists(alignmentsDir)):
    os.makedirs(alignmentsDir)

# First, do some sanity checks
found=distutils.spawn.find_executable("bowtie2");
if(found is None):
    print "please install bowtie2\n"
    sys.exit(0)

# Make sure scripts are executable
if(not os.path.exists("get-SNP-IDs.pl") or not os.path.exists("fdr.R")):
   print "please copy *.pl and *.R scripts into the current directory";
   sys.exit(1)
os.system("chmod +x *.pl *.R");

def System(cmd):
   print "Excecuting: "+cmd
   os.system(cmd)

System("get-SNP-IDs.pl "+vcfFile+" > SNPs.txt");

System("get-amplicons.pl "+chrFile+" "+ampliconsBed+" > amplicons.txt");

System("make-haplotypes.pl "+vcfFile+" haplotypes.fasta > haplotypes.txt");

System("cat haplotypes.txt | awk '{print $2 \"\t\" $1}' > amp-hap.txt");

System("get-SNP-haplotypes.pl haplotypes.txt > SNP-haplotypes.txt");

System("rm -f hap*bt2");

System("bowtie2-build haplotypes.fasta hap");

System("make-bowtie-slurms.pl "+fastqDir+" "+alignmentsDir+" "+workingDir);

print "Please run bowtie2 now -- see commands in bowtie-slurms directory or bowtie-commands.sh\n";

