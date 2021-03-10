#!/bin/bash
vcfDir=$1
chr_vcf=$2
sample=$3

mkdir -p $chr_vcf
cat $vcfDir | gunzip | grep "#" > $chr_vcf/${sample}_header.txt
cat $vcfDir | gunzip | grep -E "^[^#]" > $chr_vcf/${sample}_content.vcf
sed -i -e 's/^/chr/' $chr_vcf/${sample}_content.vcf
cat $chr_vcf/${sample}_header.txt $chr_vcf/${sample}_content.vcf > $chr_vcf/${sample}_chr.vcf
