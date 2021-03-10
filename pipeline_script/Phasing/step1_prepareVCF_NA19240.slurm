#!/bin/bash
#SBATCH --get-user-env
#SBATCH -A reddylab
#SBATCH -p all
#SBATCH -t 4-00:00:00
#SBATCH --mem=350G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --out=0926-NA19240-hg18-prepareVCF-%j.out
#SBATCH --error=0926-NA19240-prepareVCF-%j.error

module load vcftools/0.1.15-gcb01
module load glibc/2.14-gcb01

sample="NA19240"

vcf=/data/allenlab/scarlett/data/VCF/1KGP/YRI_trio/YRI.trio.2010_03.genotypes.vcf.gz
outDir=/data/allenlab/scarlett/data/shapeit2/1KGP/$sample/hg18_pilot
vcfDir=/data/allenlab/scarlett/data/shapeit2/1KGP/$sample/VCF

### pre-pare VCF (we only need to do this once!!!!)
vcftools --gzvcf $vcf --indv $sample --min-alleles 2 --max-alleles 2 --out $vcfDir/biallelic_${sample} --remove-filtered-all \
--recode-INFO-all --recode --stdout | bgzip -c > $vcfDir/biallelic_${sample}.vcf.gz
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Finish preparing the initital VCF"

for N in {1..22}
do
    cd $vcfDir
    now=$(date +"%T")
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with extracting Individual "${sample}" chr"${N}" -"${now}
    vcftools --gzvcf $vcfDir/biallelic_${sample}.vcf.gz --indv $sample --chr ${N} --recode
############ 0. only extract genotype
    cat $vcfDir/out.recode.vcf |grep -E "^[^#]" |awk '{split($10,a,/:/);$10=a[1]}1' | bgzip  > tmp.biallelic_${sample}_chr${N}.clear.vcf.gz
    rm out.recode.vcf
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> after cleaning, original VCF has this many combinations"
    cat $vcfDir/tmp.biallelic_${sample}_chr${N}.clear.vcf.gz| gunzip | grep -E "^[^#]" |awk '{split($10,a,/:/);$10=a[1]}1' | awk '{print $10}' | sort | uniq -c
    cat $vcfDir/tmp.biallelic_${sample}_chr${N}.clear.vcf.gz| gunzip |grep "#" > header.txt
    now=$(date +"%T")
    echo ">> 1.modifying VCF -"${now}
############ 1. only keeps SNPs, filter out indels
    cat $vcfDir/tmp.biallelic_${sample}_chr${N}.clear.vcf.gz | gunzip | grep -E "^[^#]" | awk 'length($4)==1 && length($5)==1' > $vcfDir/tmp.${sample}_chr${N}.SNPs.vcf
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> only keep biallelic sites, original VCF has this many combinations"
    cat $vcfDir/tmp.${sample}_chr${N}.SNPs.vcf | awk '{print $10}' | sort | uniq -c
    cat $vcfDir/tmp.${sample}_chr${N}.SNPs.vcf | grep "|" > $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.vcf
    cat $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.vcf | awk '{print $10}' | sort | uniq -c
    cat $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.vcf | awk '{ if($10=="1|0" || $10=="0|1"){print} }' > $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.hets.vcf
    #935566 0|1
    #947443 1|0
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> only keep phased hets, filtered VCF has this many combinations"
    cat $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.hets.vcf | awk '{print $10}' | sort | uniq -c
    cat $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.hets.vcf |  awk '{ sub(/\|/, "/", $10) }1' > $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> change phased to unphased, filtered VCF has this many combinations"
    cat $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf | awk '{print $10}' | sort | uniq -c
    #68492 0/1
    #69226 1/0
    wc -l $vcfDirtmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf
    #137718
    ############ 2. only keep hets, add 'chr'
    sed -i -e 's/^/chr/' $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf
    cat $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf | awk '{print $10}' | sort | uniq -c
    #68492 0/1
    #69226 1/0 
    #echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> after adding headers"
    #cat tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf | head -n 10 
    #chr1 748254 . C G . PASS . GT:GQ:DP 1/0
    #chr1 748289 . C G . PASS . GT:GQ:DP 1/0
    cat $vcfDir/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf > $vcfDir/${sample}_chr${N}.shapeit_input.vcf
   ############ 3. only select few columns
    input_vcf=$vcfDir/${sample}_chr${N}.shapeit_input.vcf
    wc -l ${input_vcf}
    sed -i '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA19240' $input_vcf
    sed -i '1i ##FILTER=<ID=PASS,Description="Passed variant FILTERs">' $input_vcf
    sed -i '1i ##source=pseq' $input_vcf
    sed -i '1i ##fileformat=VCFv4.1' $input_vcf
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> after adding chr and headers"
    cat $input_vcf | head -n 10
    rm $vcfDir/tmp*
    rm $vcfDir/out*
    rm *log
    now=$(date +"%T")
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> finish step1: chr"${N}" -"${now}
done
