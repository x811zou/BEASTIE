#!/bin/bash
fastq_R1=$1
fastq_R2=$2
trimmed_fastq=$3
sample=$4

mkdir -p $trimmed_fastq

java -jar /data/reddylab/software/Trimmomatic-0.33/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 16 -phred33 $fastq_R1 $fastq_R2 $trimmed_fastq/${sample}_FWD_paired.fq.gz $trimmed_fastq/${sample}_FWD_unpaired.fq.gz $trimmed_fastq/${sample}_REV_paired.fq.gz $trimmed_fastq/${sample}_REV_unpaired.fq.gz ILLUMINACLIP:/data/reddylab/gjohnson/reference_data/trimmomatic_MHPS.fa:2:30:10:8:TRUE LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:36
