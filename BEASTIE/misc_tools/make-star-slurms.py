#!/usr/bin/env python
import glob
import os
import re
import sys

# Global variables
jobName="STAR"
samplesDir="/data/chilab/RNAseq_2015-07"
slurmDir="/data/chilab/bill/slurm-STAR"
starIndex="/data/chilab/bill/STAR-index"
fastqFiles="/data/chilab/bill/STAR"
outputDir="/data/chilab/bill/sam"
sjdbOverhang=125
numThreads=8
memory=40000
STAR="/data/reddylab/software/STAR_2.4.2a/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR";

# Make output directories
if(not os.path.exists(slurmDir)):
    os.makedirs(slurmDir)
if(not os.path.exists(outputDir)):
    os.makedirs(outputDir);

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
  "#SBATCH -J %(jobName)s%(jobID)i" % locals(),
  "#SBATCH -o %(jobName)s%(jobID)i.output" % locals(),
  "#SBATCH -e %(jobName)s%(jobID)i.output" % locals(),
  "#SBATCH -A %(jobName)s%(jobID)i" % locals(),
  "#SBATCH --mem %(memory)i" %locals(),
  "#SBATCH --cpus-per-task=%(numThreads)s" %locals(),
  "#"])
  print >>OUT, header
  #print >>OUT, "cd "+outputDir
  print >>OUT, "cd "+outputDir

  # Process each file
  files=glob.glob(sample+"/*.fastq.gz")
  for file in files:
    match=re.search("([^/]+)\s*$",file);
    if(match is None): sys.exit("can't parse filename")
    fileNoPath=match.group(1)
    match=re.search("(\S+_R)([12])(_\S+.fastq.gz)",fileNoPath);
    if(match is None): sys.exit("can't parse paired file indicator: "+fileNoPath)
    prefix=match.group(1)
    R=int(match.group(2))
    suffix=match.group(3)
    if(R!=1): continue
    firstFile=fastqFiles+"/"+fileNoPath
    secondFile=fastqFiles+"/"+prefix+"2"+suffix
    match=re.search("(\S+).fastq.gz",fileNoPath)
    if(match is None): sys.exit("Can't parse filename")
    filestem=match.group(1)
    command=STAR+" --genomeLoad LoadAndKeep --genomeDir %(starIndex)s --readFilesIn %(firstFile)s %(secondFile)s --readFilesCommand zcat --outFileNamePrefix %(filestem)s --outSAMstrandField intronMotif --runThreadN %(numThreads)i" % locals()
    print >>OUT, command

  OUT.close()
  jobID=jobID+1

