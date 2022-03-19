#!/bin/bash

base_dir=/mnt
sample=NA06984
input_dir=$base_dir/$sample

input_vcfgz=$input_dir/$sample.no_chr.content.SNPs.hets.vcf.gz
sample_name_in_vcf=$sample
input_pileup=$input_dir/$sample.pileup.gz
input_simulation_pileup=$input_dir/$sample.simulation.pileup.gz
input_shapeit2=$input_dir/$sample.shapeit.tsv
input_hetsnp=$input_dir/${sample}_hetSNP.tsv
ancestry=GBR
read_length=75
LD_token=410113891a71
output_dir=$base_dir/$sample/beastie


### in cluster
#beastie \

### in local
PYTHONPATH='.' python3 bin/beastie \
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
    --chr-start $1 \
    --chr-end $2 \
    --STAN /usr/local/bin \
    --output-dir $output_dir \
    --ldlink-cache-dir $base_dir \
    --ldlink-token-db /mnt/ldlink_tokens.db