#!/usr/bin/env python
import glob
import os
import re

# Global variables
jobName="bam"
slurmDir="/data/chilab/bill/slurm-bam"
samDir="/data/chilab/bill/sam"
memory=10000
samtools="/usr/local/bin/samtools"

# Make output directories
if(not os.path.exists(slurmDir)):
    os.makedirs(slurmDir)

# Get list of sam files
samFiles=glob.glob(samDir+"/*.sam")

# Process each sam file
jobID=1
for samfile in samFiles:
  match=re.search("([^/]+)\.sam",samfile);
  bamfile=samDir+"/"+match.group(1)+".bam"
  bamStem=samDir+"/"+match.group(1)
  slurmFile=slurmDir+"/"+str(jobID)+".slurm"
  OUT=open(slurmFile,"w")
  header="\n".join(["#!/bin/bash",
  "#",
  "#SBATCH -J %(jobName)s%(jobID)i" % locals(),
  "#SBATCH -o %(jobName)s%(jobID)i.output" % locals(),
  "#SBATCH -e %(jobName)s%(jobID)i.output" % locals(),
  "#SBATCH -A %(jobName)s%(jobID)i" % locals(),
  "#SBATCH --mem %(memory)i" %locals(),
  "#"])
  print >>OUT, header
  print >>OUT, "cd "+samDir
  command="%(samtools)s view -bS %(samfile)s | samtools sort - %(bamStem)s" % locals()
  print >>OUT, command
  command="%(samtools)s index %(bamfile)s" % locals()
  print >>OUT, command
  OUT.close()
  jobID=jobID+1

