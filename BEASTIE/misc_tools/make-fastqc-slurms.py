#!/usr/bin/env python
import glob
import os
import re

# Global variables
samplesDir="/data/chilab/RNAseq_2015-07"
slurmDir="/data/chilab/bill/slurm-fastqc"
outputDir="/data/chilab/bill/fastqc"
fastqc="/data/chilab/bill/software/FastQC/fastqc"

# Make output directory
if(not os.path.exists(slurmDir)):
    os.makedirs(slurmDir)

# Get list of sample directories
samples=glob.glob(samplesDir+"/Sample_*")

# Process each sample
jobID=1
for sample in samples:
  match=re.search("(Sample_\S+)",sample); id=match.group(0)
  outfile=slurmDir+"/"+id+".slurm"
  OUT=open(outfile,"w")
  header="\n".join(["#!/bin/bash",
  "#",
  "#SBATCH -J FASTQC%(jobID)i" % locals(),
  "#SBATCH -o FASTQC%(jobID)i.output" % locals(),
  "#SBATCH -e FASTQC%(jobID)i.output" % locals(),
  "#SBATCH -A FASTQC%(jobID)i" % locals(),
  "#\n"])
  print >>OUT, header
  print >>OUT, "cd "+outputDir

  # Process each file
  files=glob.glob(sample+"/*.fastq.gz")
  for file in files:
    command=fastqc+" -o "+outputDir+" "+file
    print >>OUT, command

  OUT.close()
  jobID=jobID+1
