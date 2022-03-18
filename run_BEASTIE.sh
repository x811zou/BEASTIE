#!/bin/bash

input_dir=/Users/scarlett/allenlab/BEASTIE_other_example/NA06984

input_vcfgz=$input_dir/VCF/NA06984.no_chr.content.SNPs.hets.vcf.gz
sample_name_in_vcf=NA06984
input_pileup=$input_dir/mpileup/NA06984.pileup
input_simulation_pileup=$input_dir/simulation/mpileup/NA06984.pileup
input_shapeit2=$input_dir/shapeit2/NA06984.shapeit.tsv
input_hetsnp=$input_dir/hetSNP/NA06984_hetSNP.tsv
ancestry=GBR
read_length=75
LD_token=08685d209c72
output_dir=/Users/scarlett/allenlab/BEASTIE_other_example/NA06984/beastie

### in cluster
#beastie \

### in local
PYTHONPATH='.' python bin/beastie \
    --prefix test \
    --vcfgz-file $input_vcfgz \
    --vcf-sample-name $sample_name_in_vcf \
    --pileup-file $input_pileup \
    --het-snp-file $input_hetsnp \
    --shapeit2-phasing-file $input_shapeit2 \
    --simulation-pileup-file $input_simulation_pileup \
    --ancestry $ancestry \
    --read-length $read_length \
    --ld-token $LD_token \
    --chr-start 21 \
    --chr-end 22 \
    --STAN /Users/scarlett/allenlab/software/cmdstan/examples/iBEASTIE2 \
    --output-dir $output_dir