#!/bin/bash

sample="${1}"
vcfDir="${2}"
annoGtfFile="${3}"
annoGtfDir="${4}"
pipelineDir="${5:-./Process_VCF}"

vcfFile = ${vcfDir}/${sample}.vcf.recode.vcf.gz

outDir=${pipelineDir}/${sample}
mkdir -p ${outDir}

echo "step 1: split the annotation file by chromosome"
for N in {1..22}; do
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with chr${N}"
    vcftools --vcf ${annoGtfFile} --chr chr${N} --out ${annoGtfDir}/chr${N} --recode
    mv ${annoGtfDir}/chr${N}.recode.vcf ${annoGtfDir}/chr${N}.gtf
done

echo "step 2: extract het sites for mpileup${sample}"
mkdir -p ${outDir}/hetsMeta
python ./step1_extracthets.py ${sample} ${vcfFile} ${annoGtfDir} ${outDir}

echo "step 3: remove per-chromosome annotations"
for N in {1..22}; do
    rm ${annoGtfDir}/chr${N}.gtf
done
