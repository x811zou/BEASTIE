#!/bin/bash
#sh Process_RNAseq_pipeline_I_trim.sh <sample> <ref> <fastq_Dir>

sample=$1
fastq_Dir=$2
fastq_R1=$fastq_Dir/*1.fastq.gz
fastq_R2=$fastq_Dir/*2.fastq.gz
pipelineDir=/Process_RNAseq
mkdir -p $pipelineDir
mkdir -p $pipelineDir/$sample
mkdir -p $pipelineDir/$sample/trimmed_fastq
####################################
echo "step1.1: trim fastq"
sh step1_trim_fastq.sh $fastq_R1 $fastq_R2 $pipelineDir/$sample/trimmed_fastq $sample
