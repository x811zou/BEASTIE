#!/bin/bash
base_dir=/mnt
ref_dir=/mnt
LD_token=410113891a71
genotypeEr_cutoff=0.05
binomialp_cutoff=0.05
ASE_cutoff=0.5
min_single_count=1
min_total_count=1
stan=/usr/local/bin/
sigma=0.7
set -x

########################## initialization
atacseq_dir=$base_dir/ATACseq
beastie_dir=$atacseq_dir/beastie
mkdir -p $beastie_dir
atacseq_sample=test
sample=$atacseq_sample
# unphased VCF files (bgzip -c X > X.gz; tabix -p vcf X.gz)
atacseq_input_vcfgz=$atacseq_dir/${sample}.no_chr.content.SNPs.hets.vcf.gz
input_vcfgz=$atacseq_input_vcfgz
# annotation file could be generated from peaks bed file ("chr""start_pos""end_pos""regionID""geneID")
atacseq_annotationFile=$atacseq_dir/chr1_annotation.tsv.gz
# hetSNP file could be obtained from extractHets-ATACSEQ option
atacseq_hetSNP_file=$atacseq_dir/${sample}_hetSNP.tsv

########################## option1: use local python script running
#these codes are meant to use in command line
#cd /home/scarlett/github/BEASTIE
#docker run -it --entrypoint /bin/bash -v `pwd`:`pwd` -v /home/scarlett/github/BEASTIE/BEASTIE_example:/mnt xuezou/beastie
#cd /home/scarlett/github/BEASTIE/
#sh run_BEASTIE_atacseq.sh
#exit

############ (1) extractHets-ATACSEQ
# PYTHONPATH='.' python3 bin/beastie \
#     extractHetsATACseq \
#     --vcfgz-file $input_vcfgz \
#     --output $atacseq_hetSNP_file \
#     --annotationFile $atacseq_annotationFile

############ (2) filterGenotypingError & mpileup
# hetSNP_file=$atacseq_hetSNP_file
# filtered_hetSNP_file=$atacseq_dir/${sample}_hetSNP_filtered.tsv
# genotypeEr_file=$beastie_dir/${sample}_genotypeEr.tsv
# input_pileup=$atacseq_dir/${sample}.pileup.gz
# WARMUP=5000
# KEEPER=7300
# chr_start=1
# chr_end=1
# read_length=50
# filterGenotypingError_foldername="filterGenotypingError"

# PYTHONPATH='.' python3 bin/beastie \
#     filterGenotypingError \
#     --atacseq True \
#     --vcfgz-file $input_vcfgz \
#     --pileup-file $input_pileup \
#     --input-het-snp-file $hetSNP_file \
#     --filtered-het-snp-file $filtered_hetSNP_file \
#     --sample $sample \
#     --read-length $read_length \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --out-dir $beastie_dir/$filterGenotypingError_foldername \
#     --warmup $WARMUP \
#     --keeper $KEEPER 

############ (2) runModel version3: no VCF phasing, no shapeit2 --> use beastie-fix-uniform model
# filtered_hetSNP_file=$atacseq_dir/${sample}_hetSNP_filtered.tsv
# filtered_vcfgz=$atacseq_input_vcfgz
# output_dir=runModel_phased_even${read_length}_atacseq
# chr_start=1
# chr_end=1
# read_length=50

# PYTHONPATH='.' python3 bin/beastie \
#     runModel \
#     --atacseq True \
#     --vcfgz-file $filtered_vcfgz \
#     --vcf-sample-name $sample \
#     --filtered-het-snp-file $filtered_hetSNP_file \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --min-single-cov $min_single_count \
#     --min-total-cov $min_total_count \
#     --read-length $read_length \
#     --output-dir $beastie_dir/$output_dir \
#     --ldlink-cache-dir $base_dir \
#     --save-intermediate \
#     --alignBiasP-cutoff $binomialp_cutoff \
#     --ase-cutoff $ASE_cutoff \
#     --ld-token $LD_token \
#     --sigma $sigma

########################## option2: use docker running
# docker run -v `pwd`:`pwd` -v /data2:/mnt ee1449909f33ce7c74cbc0b8f25c604b27b415e87bfe6e63af4487e6aad2ce75 \
#     extractHetsATACseq \
#     --vcfgz-file $atacseq_input_vcfgz \
#     --output $atacseq_hetSNP_file \
#     --annotationFile $atacseq_annotationFile

