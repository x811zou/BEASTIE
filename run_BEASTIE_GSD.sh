#!/bin/bash
base_dir=/mnt
ref_dir=/mnt
# collected_alignmentBias_file=$base_dir/$sample/alignmentBias_list.tsv
LD_token=410113891a71
chr_start=1
chr_end=22
genotypeEr_cutoff=0.05
binomialp_cutoff=0.05
ASE_cutoff=0.5
WARMUP=5000
KEEPER=7300
stan=/usr/local/bin/
model=iBEASTIE4
sigma=0.7
min_single_count=1
min_total_count=1
### GSD individual
sample=123375
sample_name_in_vcf=$sample
beastie_outdir=$base_dir/$sample/beastie
input_vcfgz=$base_dir/$sample/tmp/${sample}.no_chr.content.SNPs.hets.vcf.gz
filtered_vcfgz=$base_dir/$sample/${sample}.no_chr.content.SNPs.hets.filtered.vcf.gz
input_pileup=$base_dir/$sample/${sample}.pileup.gz
input_shapeit2=$base_dir/$sample/${sample}.shapeit.tsv
hetSNP_file=$base_dir/${sample}_hetSNP.tsv
genotypeEr_file=$base_dir/${sample}_genotypeEr.tsv
input_simulation_pileup=$base_dir/$sample/${sample}.simulation_even_100.pileup.gz
filtered_hetSNP_file=$base_dir/$sample/${sample}_hetSNP_filtered.tsv
ancestry=EUR
read_length=51
filterGenotypingError_foldername="filterGenotypingError"
mkdir -p $beastie_outdir/$filterGenotypingError_foldername/${filterGenotypingError_foldername}
set -x

########################## use local python script running
# cd /home/scarlett/github/BEASTIE
# docker run -it --entrypoint /bin/bash -v `pwd`:`pwd` -v /data2/GSD:/mnt xuezou/beastie
# cd /home/scarlett/github/BEASTIE/
# sh run_BEASTIE_GSD.sh
# exit

############ (1) extractHets
#docker run -v `pwd`:`pwd` -v /data2:/mnt xuezou/beastie \
# PYTHONPATH='.' python3 bin/beastie \
#     extractHets \
#     --vcfgz-file $input_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --output $hetSNP_file \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --gencode-dir $ref_dir/reference/gencode_chr \
#--skip-require-pass # only for UDN

# ############# filterGenotypingError
# docker run -v `pwd`:`pwd` -v /home/scarlett/data:/mnt xuezou/beastie \
# PYTHONPATH='.' python3 bin/beastie \
#     filterGenotypingError \
#     --filtered-het-snp-file $filtered_hetSNP_file \
#     --genotype-error-file $genotypeEr_file \
#     --vcfgz-file $input_vcfgz \
#     --sample $sample \
#     --pileup-file $input_pileup \
#     --input-het-snp-file $hetSNP_file \
#     --ancestry $ancestry \
#     --read-length $read_length \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --af-dir $ref_dir/reference/AF \
#     --genotypeEr-cutoff $genotypeEr_cutoff \
#     --out-dir $beastie_outdir/$filterGenotypingError_foldername/${filterGenotypingError_foldername} \
#     --warmup $WARMUP \
#     --keeper $KEEPER 

############# runModel biasSNPlist
# PYTHONPATH='.' python3 bin/beastie \
#     runModel \
#     --vcfgz-file $filtered_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --filtered-het-snp-file $filtered_hetSNP_file \
#     --ancestry $ancestry \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --read-length $read_length \
#     --STAN $stan \
#     --output-dir $beastie_outdir/beastie_shapeit2_even_biasSNPlist \
#     --ldlink-cache-dir $base_dir \
#     --save-intermediate \
#     --binomialp-cutoff $binomialp_cutoff \
#     --ase-cutoff $ASE_cutoff \
#     --ld-token $LD_token \
#     --shapeit2-phasing-file $input_shapeit2 \
#     --simulation-pileup-file $input_simulation_pileup \
#     --collected-alignmentBias-file $collected_alignmentBias_file

############# runModel version1: VCF phasing
# docker run -v `pwd`:`pwd` -v /data2:/mnt b3c496fe9dfe05a4acd08f66160aa97465bac90a99afbedf04a415767a7 \
# PYTHONPATH='.' python3 bin/beastie \
#     runModel \
#     --vcfgz-file $filtered_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --filtered-het-snp-file $filtered_hetSNP_file \
#     --ancestry $ancestry \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --min-single-cov $min_single_count \
#     --min-total-cov $min_total_count \
#     --read-length $read_length \
#     --output-dir $beastie_outdir/runModel_phased_even${read_length}_iBEASTIE4_VCF \
#     --ldlink-cache-dir $base_dir \
#     --save-intermediate \
#     --alignBiasP-cutoff $binomialp_cutoff \
#     --ase-cutoff $ASE_cutoff \
#     --ld-token $LD_token \
#     --model $model \
#     --gam-model-name iBEASTIE4_s0.7_GAM/gam4_lambdamodel.pkl \
#     --sigma $sigma \
#     --simulation-pileup-file $input_simulation_pileup
   # --STAN /home/scarlett/github/BEASTIE/BEASTIE

############# runModel version2: shapeit2 phasing
# cd /home/scarlett/github/BEASTIE
# docker run -it --entrypoint /bin/bash -v `pwd`:`pwd` -v /data2/GSD:/mnt xuezou/beastie
# cd /home/scarlett/github/BEASTIE/
# sh run_BEASTIE_GSD.sh

#docker run -v `pwd`:`pwd` -v /data2:/mnt xuezou/beastie \
PYTHONPATH='.' python3 bin/beastie \
    runModel \
    --vcfgz-file $filtered_vcfgz \
    --vcf-sample-name $sample_name_in_vcf \
    --filtered-het-snp-file $filtered_hetSNP_file \
    --ancestry $ancestry \
    --chr-start $chr_start \
    --chr-end $chr_end \
    --min-single-cov $min_single_count \
    --min-total-cov $min_total_count \
    --read-length $read_length \
    --output-dir $beastie_outdir/runModel_phased_even${read_length}_iBEASTIE4 \
    --ldlink-cache-dir $base_dir \
    --save-intermediate \
    --alignBiasP-cutoff $binomialp_cutoff \
    --ase-cutoff $ASE_cutoff \
    --ld-token $LD_token \
    --model $model \
    --gam-model-name iBEASTIE4_s0.7_GAM/gam4_lambdamodel.pkl \
    --sigma $sigma \
    --shapeit2-phasing-file $input_shapeit2 \
    --simulation-pileup-file $input_simulation_pileup
    #--STAN /home/scarlett/github/BEASTIE/BEASTIE

############# runModel version3: no VCF phasing, no shapeit2 --> use beastie-fix-uniform model
# PYTHONPATH='.' python3 bin/beastie \
#     runModel \
#     --nophasing True \
#     --vcfgz-file $filtered_vcfgz \
#     --vcf-sample-name $sample_name_in_vcf \
#     --filtered-het-snp-file $filtered_hetSNP_file \
#     --ancestry $ancestry \
#     --chr-start $chr_start \
#     --chr-end $chr_end \
#     --min-single-cov $min_single_count \
#     --min-total-cov $min_total_count \
#     --read-length $read_length \
#     --output-dir $beastie_outdir/runModel_phased_even${read_length}_BEASTIEFIXUNIFORM_nophasing \
#     --ldlink-cache-dir $base_dir \
#     --save-intermediate \
#     --alignBiasP-cutoff $binomialp_cutoff \
#     --ase-cutoff $ASE_cutoff \
#     --ld-token $LD_token \
#     --model $model \
#     --gam-model-name iBEASTIE3_s0.7_GAM/gam4_lambdamodel.pkl \
#     --sigma $sigma \
#     --simulation-pileup-file $input_simulation_pileup \
#     --STAN /home/scarlett/github/BEASTIE/BEASTIE

########################## use docker directly then use this line instead of PYTHON xxx
# docker run -v `pwd`:`pwd` -v /data2:/mnt xuezou/beastie \