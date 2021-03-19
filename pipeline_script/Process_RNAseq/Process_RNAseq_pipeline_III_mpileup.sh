#!/bin/bash
#sh 
module load samtools/1.3.1-gcb01

sample=$1

#PipelineDir=/data/allenlab/scarlett/result/project_ASE/Process_RNAseq
OutDir=${PipelineDir}/${sample}/star2pss_EndtoEnd_wasp/mismatch10/mpileup
mkdir -p $OutDir

ref=$2
#ref=/data/reddylab/Reference_Data/Genomes/hg19/hg19.fa
hetsMeta_mpileup=$3
inputDir=./Process_VCF/${sample}/hetsMeta
Bam=$4
RNAseqDir=${PipelineDir}/${sample}/star2pss_EndtoEnd_wasp/mismatch${N}
Bam=$RNAseqDir/aligned.bam
pileup_out=$OutDir/Allchr_hets_all_transcript.pileup
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> We start with sample: ${sample}"
cd $OutDir
samtools mpileup -d 0 -B -s -f $ref -l $inputDir/hetmeta_mpileup.tsv $Bam > $pileup_out

