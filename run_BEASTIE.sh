#!/bin/bash

#base_dir=/mnt
base_dir=/Users/scarlett/allenlab/BEASTIE_other_example
ref_dir=/Users/scarlett/BEASTIE/BEASTIE
sample=HG00171
sample_name_in_vcf=$sample
beastie_outdir=$base_dir/$sample/beastie_noshapeit2_random
input_vcfgz=$base_dir/$sample/$sample.no_chr.content.SNPs.hets.vcf.gz
filtered_vcfgz=$base_dir/$sample/$sample.no_chr.content.SNPs.hets.filtered.vcf.gz
input_pileup=$base_dir/$sample/$sample.pileup.gz
input_simulation_pileup=$base_dir/$sample/$sample.simulation_random.pileup.gz
input_shapeit2=$base_dir/$sample/$sample.shapeit.tsv
hetSNP_file=$base_dir/$sample/${sample}_hetSNP.tsv
genotypeEr_file=$base_dir/$sample/${sample}_genotypeEr.tsv
filtered_hetSNP_file=$base_dir/$sample/${sample}_filtered_hetSNP.tsv
ancestry=FIN
read_length=75
LD_token=410113891a71
chr_start=1
chr_end=22
genotypeEr_cutoff=0.05
binomialp_cutoff=0.05
ASE_cutoff=0.5
### in cluster
#beastie \

# PYTHONPATH='.' python3 bin/beastie \
#     extractHets \
#     --output $hetSNP_file \
#     --vcfgz-file $input_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --gencode-dir $ref_dir/reference/gencode_chr 

PYTHONPATH='.' python3 bin/beastie \
    filterGenotypingError \
    --vcfgz-file $input_vcfgz \
    --vcf-sample-name $sample_name_in_vcf \
    --pileup-file $input_pileup \
    --input-het-snp-file $hetSNP_file \
    --ancestry $ancestry \
    --read-length $read_length \
    --chr-start $chr_start \
    --chr-end $chr_end \
    --beastie-outdir $beastie_outdir \
    --af-dir $ref_dir/reference/AF \
    --outfile $genotypeEr_file \
    --genotypeEr-cutoff $genotypeEr_cutoff \
    --filtered-het-snp-file $filtered_hetSNP_file

# PYTHONPATH='.' python3 bin/beastie \
#     runModel \
#     --vcfgz-file $input_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --simulation-pileup-file $input_simulation_pileup \
#     --filtered-het-snp-file $filtered_hetSNP_file \
#     --ancestry $ancestry \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --STAN /Users/scarlett/allenlab/software/cmdstan/examples/iBEASTIE2 \
#     --output-dir $beastie_outdir \
#     --ldlink-cache-dir $base_dir \
#     --save-intermediate \
#     --binomialp-cutoff $binomialp_cutoff \
#     --ase-cutoff $ASE_cutoff \
#     --ld-token $LD_token \
#     --shapeit2-phasing-file $input_shapeit2 
#     # --ldlink-token-db /Users/scarlett/allenlab/BEASTIE_other_example/ldlink_tokens.db
