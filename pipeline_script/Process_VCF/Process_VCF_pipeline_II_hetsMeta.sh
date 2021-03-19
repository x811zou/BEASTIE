#!/bin/bash

module load python/3.7.4-gcb01

sample=$1
vcfDir=$2
PipelineDir=/Process_VCF
#vcfDir=/data/allenlab/scarlett/data/VCF/GSD/DNA_vcf/${sample}.vcf.recode.vcf.gz
#PipelineDir=/data/allenlab/scarlett/result/project_ASE/Process_VCF
mkdir -p $PipelineDir
outDir=$PipelineDir/$sample
mkdir -p $outDir
echo "step1: extract het sites for mpileup"${sample}
mkdir -p $outDir/hetsMeta
mkdir -p $outDir/hetsDict
mkdir -p $outDir/hetsDict/genom_all
mkdir -p $outDir/hetsDict/trans_all
python ./Process_VCF/step1_extracthets.py $sample $vcfDir $outDir  

#### if input sample is a list
#for sample in "123375" "125249" "125260" "122687" "122698" "123667"
#do
#done


