#!/bin/bash

# need you to specify
base_dir=/mnt
ref_dir=/mnt

# pre-specified parameters for the example data NA12878_chr21
sample=NA12878_chr21
sample_name_in_vcf=HG001
ancestry=CEU
read_length=101
LD_token=410113891a71 # This is for testing only. You could use your own LD token.
chr_start=21
chr_end=21
input_dir=/home/scarlett/github/BEASTIE/BEASTIE_example/${sample}
original_vcfgz=$input_dir/${sample}.vcf.gz
original_bamgz=$input_dir/${sample}.bam.gz
input_simulation_pileup_gz=$input_dir/${sample}/${sample}.simulation_even_101.pileup.tsv.gz
input_shapeit2=$input_dir/${sample}.shapeit.tsv

# No need to specificy (default setting)
### directory
out_dir=$base_dir/BEASTIE_example_output/$sample
beastie_outdir=$out_dir/beastie
tmp_genotypeEr_dir=$out_dir/tmp/genotyping_error
tmp_simulation_dir=$out_dir/tmp/simulation
tmp_mpileup_dir=$out_dir/tmp/mpileup
tmp_vcf_dir=$out_dir/tmp/vcf
### output
bi_vcfgz=$out_dir/${sample}.bi.vcf.gz
bihets_vcfgz=$out_dir/${sample}.bihets.vcf.gz
out_pileup_gz=$out_dir/${sample}.pileup.tsv.gz
out_hetSNP=$out_dir/${sample}.hetSNP.tsv
out_hetSNP_filtered=$out_dir/${sample}.hetSNP.filtered.tsv
out_genotypeEr=$out_dir/${sample}.genotypeEr.tsv
out_fwd_fastqreads=$out_dir/${sample}.fwd.fastq.gz
out_rev_fastqreads=$out_dir/${sample}.rev.fastq.gz

### Test all pipeline stpes at a time
############ (0) cleanVCF
docker run  -v `pwd`:`pwd` -v /data2:/mnt xuezou/beastie
    cleanVCF \
    --tmp-dir $tmp_vcf_dir \
    --vcfgz-file $original_vcfgz \
    --filtered-bi-vcfgz $bi_vcfgz \
    --filtered-bihets-vcfgz $bihets_vcfgz

# PYTHONPATH='.' python3 bin/beastie \
#     extractHets \
#     --output $hetSNP_file \
#     --vcfgz-file $input_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --gencode-dir $ref_dir/reference/gencode_chr 

# PYTHONPATH='.' python3 bin/beastie \
#     filterGenotypingError \
#     --vcfgz-file $input_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --pileup-file $input_pileup \
#     --input-het-snp-file $hetSNP_file \
#     --ancestry $ancestry \
#     --read-length $read_length \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --beastie-outdir $beastie_outdir \
#     --af-dir $ref_dir/reference/AF \
#     --outfile $genotypeEr_file \
#     --genotypeEr-cutoff $genotypeEr_cutoff \
#     --filtered-het-snp-file $filtered_hetSNP_file

# PYTHONPATH='.' python3 bin/beastie \
#     runModel \
#     --vcfgz-file $input_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --simulation-pileup-file $input_simulation_pileup \
#     --filtered-het-snp-file $filtered_hetSNP_file \
#     --ancestry $ancestry \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --STAN /usr/local/bin/iBEASTIE \
#     --output-dir $beastie_outdir \
#     --ldlink-cache-dir $base_dir \
#     --save-intermediate \
#     --binomialp-cutoff $binomialp_cutoff \
#     --ase-cutoff $ASE_cutoff \
#     --ld-token $LD_token \
#     --shapeit2-phasing-file $input_shapeit2 
