#!/bin/bash
#SBATCH --get-user-env
#SBATCH --job-name=myjobR
#SBATCH -A reddylab   
#SBATCH -p new,all
#SBATCH --nice=500
#SBATCH --mem=10G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL 
#SBATCH --error=MA-1-%j.out
#SBATCH --output=MA-1-%J.out

inputBam=/data/reddylab/GSD/pacbio_MA1/r54089_20201109_134937/m54089_201109_140212.subreads.bam
outDir=/data/reddylab/scarlett/GSD/Pacbio/bam2fastq
outFilename="MA-1"

/data/reddylab/scarlett/software/bam2fastq/bam2fastq -o $outDir/${outFilename}# ${inputBam}

mv ${outFilename}_M ${outFilename}.fastq
# reference
# https://github.com/PacificBiosciences/bam2fastx
