#!/bin/bash
#sh step2_phasing.sh <sample>
module load vcftools/0.1.15-gcb01
module load glibc/2.14-gcb01

sample=$1
vcf=$2
outDir=$3
ancestry=$4

if $ancestry == "EUR":
    echo "EUR" > $outDir/group_${ancestry}.list
else:
echo "AFR" > $outDir/group_${ancestry}.list

### pre-pare VCF
cd $outDir
mkdir -p ${sample}
cd $outDir/${sample}
#mkdir -p chr${N} 
#cd $outDir/${sample}/chr${N}
 #   mkdir -p vcf
    #mkdir -p phase_withoutseq
  #  now=$(date +"%T")
    #cd $outDir
#now=$(date +"%T")
#echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with extracting Individual "${sample}" chr"${N}" -"${now}
cd $outDir/${sample}
vcftools --gzvcf $vcf --indv $sample --min-alleles 2 --max-alleles 2 --out $outDir/${sample}/chr${N}/vcf/1KGP_biallelic_${sample} --remove-filtered-all --recode-INFO-all --recode --stdout | bgzip -c > $outDir/${sample}/1KGP_biallelic_${sample}.vcf.gz
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Finish preparing the initital VCF"

for N in {1..22}
do
    cd $outDir/${sample}
    mkdir -p chr${N}
    cd $outDir/${sample}/chr${N}
    mkdir -p vcf
    mkdir -p phase_withoutseq
    #now=$(date +"%T")
    #cd $outDir
    now=$(date +"%T")
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with extracting Individual "${sample}" chr"${N}" -"${now}
    cd $outDir/${sample}/chr${N}/vcf
    vcftools --gzvcf $outDir/${sample}/1KGP_biallelic_${sample}.vcf.gz --indv $sample --chr ${N} --recode
    #vcftools --gzvcf $vcf --indv $sample --min-alleles 2 --max-alleles 2 --out $outDir/${sample}/chr${N}/vcf/1KGP_biallelic_${sample}_chr${N} --remove-filtered-all --recode-INFO-all --recode --stdout | bgzip -c > $outDir/${sample}/chr${N}/vcf/1KGP_biallelic_${sample}_chr${N}.vcf.gz
    cat $outDir/${sample}/chr${N}/vcf/out.recode.vcf |grep -E "^[^#]" |awk '{split($10,a,/:/);$10=a[1]}1' | bgzip  > tmp.1KGP_biallelic_${sample}_chr${N}.clear.vcf.gz
    rm out.recode.vcf
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> after cleaning, original VCF has this many combinations"
    cat $outDir/${sample}/chr${N}/vcf/tmp.1KGP_biallelic_${sample}_chr${N}.clear.vcf.gz| gunzip | grep -E "^[^#]" |awk '{split($10,a,/:/);$10=a[1]}1' | awk '{print $10}' | sort | uniq -c
#  13595 0/1
#  90406 0|0
#  68492 0|1
#  13581 1/0
#  69226 1|0
#  83439 1|1
    cat $outDir/${sample}/chr${N}/vcf/tmp.1KGP_biallelic_${sample}_chr${N}.clear.vcf.gz| gunzip |grep "#" > header.txt
  
    #1721929 0|0
    #27250 0|1
    #23504 1|0
    #31186 1|1
    #cat $outDir/chr${N}/vcf/biallelic_${sample}_chr${N}.vcf.gz | gunzip | grep "#" > $outDir/chr${N}/vcf/biallelic_${sample}_chr${N}.header.vcf 
    now=$(date +"%T")
    echo ">> 1.modifying VCF -"${now}
    ############ 1. only keeps SNPs, filter out indels
    cat $outDir/${sample}/chr${N}/vcf/tmp.1KGP_biallelic_${sample}_chr${N}.clear.vcf.gz | gunzip | grep -E "^[^#]" | awk 'length($4)==1 && length($5)==1' > $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.vcf
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> only keep biallelic sites, original VCF has this many combinations"
    cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.vcf | awk '{print $10}' | sort | uniq -c
# 172222 0/1
#1165248 0|0
# 935566 0|1
# 172580 1/0
#      1 1/1
# 947443 1|0
#1109379 1|1
    cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.vcf | grep "|" > $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.vcf
    cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.vcf | awk '{print $10}' | sort | uniq -c
    cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.vcf | awk '{ if($10=="1|0" || $10=="0|1"){print} }' > $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.hets.vcf
 #935566 0|1
 #947443 1|0
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> only keep phased hets, filtered VCF has this many combinations"
    cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.hets.vcf | awk '{print $10}' | sort | uniq -c
    cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.hets.vcf |  awk '{ sub(/\|/, "/", $10) }1' > $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> change phased to unphased, filtered VCF has this many combinations"
    cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf | awk '{print $10}' | sort | uniq -c
# 935566 0/1
# 947443 1/0
    wc -l $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf
    #cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.switch.vcf | awk '{ sub("1/0", "0/1", $10) }1' > $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.switch.vcf
    #1739315

    ############ 2. only keep hets, add 'chr'
    #cat $outDir/${sample}/chr${N}/vcf/tmp_${sample}_chr${N}.SNPs.phased.vcf | grep -P '([0-9])\|(?!\1)[0-9]' > $outDir/${sample}/chr${N}/vcf/tmp_${sample}_chr${N}.SNPs.phased.hets.vcf
    #wc -l $outDir/${sample}/chr${N}/vcf/tmp_${sample}_chr${N}.SNPs.phased.hets.vcf
    # 44762
    sed -i -e 's/^/chr/' $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf
    #cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf | awk '{print $10}' | sort | uniq -c
    #24014 0|1
    #20748 1|0
    #echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> after adding headers"
    #cat tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf | head -n 10 
    cat $outDir/${sample}/chr${N}/vcf/tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf > $outDir/${sample}/chr${N}/vcf/${sample}_chr${N}.shapeit_input.vcf
   ############ 3. only select few columns
    #cat $outDir/${sample}/chr${N}/vcf/tmp_${sample}_chr${N}.SNPs.phased.vcf | grep -E "^[^#]" | awk '{print $1,$3,$2,$4,$5,$8}' > $outDir/${sample}/chr${N}/vcf/tmp_${sample}_chr${N}.SNPs.phased.hets.select.vcf
    #cp $outDir/${sample}/chr${N}/vcf/tmp_${sample}_chr${N}.SNPs.phased.hets.select.vcf $outDir/${sample}/chr${N}/vcf/1KGP_AF_${sample}_biSNPs_chr${N}_phased.vcf
    #wc -l $outDir/${sample}/chr${N}/vcf/1KGP_AF_${sample}_biSNPs_chr${N}_phased.vcf
    #44762
    input_vcf=$outDir/${sample}/chr${N}/vcf/${sample}_chr${N}.shapeit_input.vcf
    wc -l ${input_vcf}
    sed -i '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA19240' $input_vcf
    sed -i '1i ##FILTER=<ID=PASS,Description="Passed variant FILTERs">' $input_vcf
    sed -i '1i ##source=pseq' $input_vcf
    sed -i '1i ##fileformat=VCFv4.1' $input_vcf
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> after adding chr and headers"
    cat $input_vcf | head -n 10
    rm $outDir/${sample}/chr${N}/vcf/tmp*
    rm $outDir/${sample}/chr${N}/vcf/out*
    #rm *log
   
    now=$(date +"%T")
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> finish step1: chr"${N}" -"${now}
    now=$(date +"%T")
    echo ">> start phasing -"${now}
    refDir=/data/allenlab/scarlett/data/shapeit2/1000GP_Phase3
    input_hap=$refDir/1000GP_Phase3_chr${N}.hap.gz
    input_legend=$refDir/1000GP_Phase3_chr${N}.legend.gz
    input_sample=$refDir/1000GP_Phase3.sample
    input_genetic_map=$refDir/genetic_map_chr${N}_combined_b37.txt

    ### 2.1.1 Phase without RNA
    #echo "AFR" > group.list
    # step1: 
    /data/allenlab/scarlett/software/shapeit2/bin/shapeit \
    -check \
    -V $outDir/${sample}/chr${N}/vcf/${sample}_chr${N}.shapeit_input.vcf \
    -M $input_genetic_map \
    --input-ref $input_hap $input_legend $input_sample \
    --output-log $outDir/${sample}/chr${N}/phase_withoutseq/${sample}_chr${N}.AlignmentChecks
    # step2: 
    /data/allenlab/scarlett/software/shapeit2/bin/shapeit \
    -V $outDir/${sample}/chr${N}/vcf/${sample}_chr${N}.shapeit_input.vcf  \
    -M $input_genetic_map \
    --input-ref $input_hap $input_legend $input_sample \
    --include-grp $outDir/group_AFR.list \
    --exclude-snp $outDir/${sample}/chr${N}/phase_withoutseq/${sample}_chr${N}.AlignmentChecks.snp.strand.exclude \
    -O $outDir/${sample}/chr${N}/phase_withoutseq/${sample}_chr${N}.phased.with.ref \
    --output-log $outDir/${sample}/chr${N}/phase_withoutseq/${sample}_chr${N}.Phasing \
    --no-mcmc

    wc -l $outDir/${sample}/chr${N}/phase_withoutseq/${sample}_chr${N}.phased.with.ref.haps #


done
