#!/bin/bash
# ./Process_RNAseq_pipeline_I_trim.sh <trimmomatic jar location> <illumina adapter fasta> <sample> <fastq dir> <processing output dir>

trimmomatic=${1}
illuminaAdapters=${2}
sample=${3}
fastqDir=${4}
pipelineDir=${5:-./Process_RNAseq}

processingOutputDir="${pipelineDir}/${sample}/trimmed_fastq"
mkdir -p ${processingOutputDir}

####################################
echo "step1.1: trim fastq"
java -jar ${trimmomatic} PE -threads 16 -phred33 "${fastqDir}"/*1.fastq.gz "${fastqDir}"/*2.fastq.gz \
    ${processingOutputDir}/${sample}_FWD_paired.fq.gz \
    ${processingOutputDir}/${sample}_FWD_unpaired.fq.gz \
    ${processingOutputDir}/${sample}_REV_paired.fq.gz \
    ${processingOutputDir}/${sample}_REV_unpaired.fq.gz \
    ILLUMINACLIP:${illuminaAdapters}:2:30:10:8:TRUE \
    LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:36

