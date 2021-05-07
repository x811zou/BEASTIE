# sh step1_prepareVCF.sh <sample> <./path/vcf.gz> <./outDir>

sample="${1}"
vcf="${2}"
vcfDir="${3}"

echo ">>>>>>>>>>>>>>>> Beginning Phasing, step 1: Prepare the VCFs -$(date +'%T')"

cd ${vcfDir}

echo ">>>>>>>>>>>>>>>> Preparing the initial VCF"
vcftools --gzvcf ${vcf} --indv ${sample} --min-alleles 2 --max-alleles 2 --out biallelic_${sample} --remove-filtered-all \
--recode-INFO-all --recode --stdout | bgzip -c > biallelic_${sample}.vcf.gz
echo ">>>>>>>>>>>>>>>> Finished preparing the initial VCF"

for N in {1..22}
do
    echo ">>>>>>>>>>>>>>>> Start with extracting individual ${sample} chr${N} -$(date +'%T')"
    vcftools --gzvcf biallelic_${sample}.vcf.gz --indv $sample --chr ${N} --recode

    ############ 0. Extract genotype
    echo ">>>>>>>>>>>>>>>> Extracting genotype from sample info"
    awk '/^[^#]/ {split($10,a,/:/);$10=a[1]}1' out.recode.vcf | bgzip > tmp.biallelic_${sample}_chr${N}.clear.vcf.gz
    rm out.recode.vcf
    echo ">>>>>>>>>>>>>>>> Original VCF has this many combinations"
    gunzip --stdout tmp.biallelic_${sample}_chr${N}.clear.vcf.gz | awk '/^[^#]/ {split($10,a,/:/); print a[1]}' | sort | uniq -c

    ############ 1. Only keep SNPs, filter out indels
    echo ">>>>>>>>>>>>>>>> Phasing 1.${N}.1. Modifying VCF to keep SNPs only (filtering out indels) -$(date +'%T')"

    gunzip --stdout tmp.biallelic_${sample}_chr${N}.clear.vcf.gz | awk '/^[^#]/ && length($4)==1 && length($5)==1' > tmp.${sample}_chr${N}.SNPs.vcf
    echo ">>>>>>>>>>>>>>>> Only keep biallelic sites, original VCF has this many combinations"
    awk '{print $10}' tmp.${sample}_chr${N}.SNPs.vcf | sort | uniq -c
    grep "|" tmp.${sample}_chr${N}.SNPs.vcf > tmp.${sample}_chr${N}.SNPs.phased.vcf
    awk '{print $10}' tmp.${sample}_chr${N}.SNPs.phased.vcf | sort | uniq -c
    awk '{ if($10=="1|0" || $10=="0|1"){print} }' tmp.${sample}_chr${N}.SNPs.phased.vcf > tmp.${sample}_chr${N}.SNPs.phased.hets.vcf
    echo ">>>>>>>>>>>>>>>> Only keep phased hets, filtered VCF has this many combinations"
    awk '{print $10}' tmp.${sample}_chr${N}.SNPs.phased.hets.vcf | sort | uniq -c
    awk '{ sub(/\|/, "/", $10) }1' tmp.${sample}_chr${N}.SNPs.phased.hets.vcf > tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf
    echo ">>>>>>>>>>>>>>>> Change phased to unphased, filtered VCF has this many combinations"
    awk '{print $10}' tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf | sort | uniq -c

    ############ 2. keep hets, add 'chr'
    echo ">>>>>>>>>>>>>>>> Phasing 1.${N}.2. Prefixing chromosome numbers with 'chr' -$(date +'%T')"
    sed -i -e 's/^/chr/' tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf
    cp tmp.${sample}_chr${N}.SNPs.phased.hets.switch.vcf ${sample}_chr${N}.shapeit_input.vcf

    ############ 3. Select a few columns
    echo ">>>>>>>>>>>>>>>> Phasing 1.${N}.3. Adding VCF headers -$(date +'%T')"
    input_vcf=${sample}_chr${N}.shapeit_input.vcf
    echo ">>>>>>>>>>>>>>>> Adding headers"
    sed -i "1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample}" ${input_vcf}
    sed -i '1i ##FILTER=<ID=PASS,Description="Passed variant FILTERs">' ${input_vcf}
    sed -i '1i ##source=pseq' ${input_vcf}
    sed -i '1i ##fileformat=VCFv4.1' ${input_vcf}
    rm tmp*
    rm out*
    rm *log
    echo ">>>>>>>>>>>>>>>> Finished preparing chr${N} vcf -$(date +'%T')"
done

echo ">>>>>>>>>>>>>>>> Finished phasing, step 1: Prepare the VCFs -$(date +'%T')"