#!/bin/bash

#base_dir=/mnt
base_dir=/Users/scarlett/allenlab/BEASTIE_other_example
ref_dir=/Users/scarlett/BEASTIE/BEASTIE
sample=HG00171
input_dir=$base_dir/$sample
input_vcfgz=$input_dir/$sample.no_chr.content.SNPs.hets.vcf.gz
filtered_vcfgz=$input_dir/$sample.no_chr.content.SNPs.hets.filtered.vcf.gz
sample_name_in_vcf=$sample
input_pileup=$input_dir/$sample.pileup.gz
input_simulation_pileup=$input_dir/$sample.simulation_random.pileup.gz
input_shapeit2=$input_dir/$sample.shapeit.tsv
input_hetsnp=$input_dir/${sample}_hetSNP.tsv
ancestry=FIN
read_length=75
LD_token=410113891a71
output_dir=$input_dir/beastie_noshapeit2_random
hetSNP_filename=$base_dir/$sample/${sample}_hetSNP.tsv
chr_start=1
chr_end=22
genotypeEr_cutoff=0.05
binomialp_cutoff=0.05
ASE_cutoff=0.5
genotypeEr_filename=${sample}.genotypeEr.tsv
### in cluster
#beastie \

# PYTHONPATH='.' python3 bin/beastie \
#     extractHets \
#     --output $hetSNP_filename \
#     --vcfgz-file $input_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --gencode-dir $ref_dir/reference/gencode_chr 

# PYTHONPATH='.' python3 bin/beastie \
#     filterGenotypingError \
#     --prefix $sample \
#     --vcfgz-file $input_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --pileup-file $input_pileup \
#     --het-snp-file $hetSNP_filename \
#     --ancestry $ancestry \
#     --read-length $read_length \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --input-dir $input_dir \
#     --output-dir $output_dir \
#     --af-dir $ref_dir/reference/AF \
#     --filename $genotypeEr_filename \
#     --genotypeEr-cutoff $genotypeEr_cutoff

PYTHONPATH='.' python3 bin/beastie \
    runModel \
    --prefix $sample \
    --vcfgz-file $input_vcfgz \
    --vcf-sample-name $sample_name_in_vcf \
    --simulation-pileup-file $input_simulation_pileup \
    --ancestry $ancestry \
    --chr-start $chr_start \
    --chr-end $chr_end \
    --STAN /Users/scarlett/allenlab/software/cmdstan/examples/iBEASTIE2 \
    --output-dir $output_dir \
    --ldlink-cache-dir $base_dir \
    --save-intermediate \
    --binomialp-cutoff $binomialp_cutoff \
    --ase-cutoff $ASE_cutoff \
    --ld-token $LD_token
#     # --ldlink-token-db /Users/scarlett/allenlab/BEASTIE_other_example/ldlink_tokens.db
