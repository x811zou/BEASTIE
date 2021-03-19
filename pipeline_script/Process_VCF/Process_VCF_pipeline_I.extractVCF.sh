#!/bin/bash

module load vcftools/0.1.15-gcb01

#######################
#cat HG00096.vcf | grep -E "^[^#]" |awk '{print $1}' | uniq | wc -l
#--recode-INFO <string>
#--recode-INFO-all
#These options can be used with the above recode options to define an INFO key name to keep in the output file. This option can be used multiple times to keep more of the INFO fields. The second option is used to keep all INFO values in the original file.
#cat HG00096.with_chr.het.vcf | grep -E '0\|1|1\|0|0\|2|0\|3|0\|4|0\|5|0\|6|1\|2|1\|3|1\|4|1\|5|2\|0|2\|1|2\|3|2\|4|2\|5|3\|0|3\|1|3\|2|3\|4|3\|5|4\|0|4\|1|4\|2|4\|3|4\|5|5\|0|5\|1|5\|2|5\|3|5\|4|6\|0|6\|1|7\|0|7\|3|7\|5' > hets.txt
#([0-9])\|(?!\1)+[0-9]
#cat HG00096.info.only_hets.vcf |  grep -E "^[^#]" | awk '{print $10}' | sort |uniq -
########################

sample=$1
vcfDir=$2
pipelineDir=/Process_VCF
#vcfDir=/data/common/1000_genomes/VCF/20130502/bgzip
#pipelineDir=/data/allenlab/scarlett/result/project_ASE/Process_VCF

mkdir -p $pipelineDir/${sample}
outDir=$pipelineDir/${sample}/VCF
mkdir -p $outDir
N=1
cd $outDir
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>. start with chr"${N}
vcftools --gzvcf $vcfDir/ALL.chr${N}.*.vcf.gz --indv $sample --out $outDir/${sample}_chr${N} --recode-INFO-all --recode
echo "1. original file has lines:"
wc -l $outDir/${sample}_chr${N}.recode.vcf

cat $outDir/${sample}_chr${N}.recode.vcf | grep "#" > $outDir/header.vcf
echo "2. header lines number: "
wc -l $outDir/header.vcf

cat $outDir/${sample}_chr${N}.recode.vcf | grep -E "^[^#]" > $outDir/content.vcf
echo "3. content lines number: "
wc -l $outDir/content.vcf

echo "add 'chr' to first column"
sed -i -e 's/^/chr/' $outDir/content.vcf
cat  $outDir/header.vcf $outDir/content.vcf >  $outDir/tmp.vcf

echo "4. CHECK: concatenate files with lines number: "
wc -l $outDir/tmp.vcf

echo ">>>>>>>> Done with chr"${N}

rm $outDir/${sample}_chr${N}.recode.vcf

for N in {2..22}
do
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>. start with chr"${N}
    vcftools --gzvcf $vcfDir/ALL.chr${N}.*.vcf.gz --indv $sample --out $outDir/${sample}_chr${N} --recode-INFO-all --recode
    echo "1. original file has lines:"
    wc -l $outDir/${sample}_chr${N}.recode.vcf

    cat $outDir/${sample}_chr${N}.recode.vcf | grep -E "^[^#]" > $outDir/tmp${N}.vcf
    echo "2. content lines number: "
    wc -l $outDir/tmp${N}.vcf

    echo "3. add 'chr' to first column"
    sed -i -e 's/^/chr/' $outDir/tmp${N}.vcf
    wc -l $outDir/tmp${N}.vcf

    echo "4. CHECK: concatenate previous files, and now  with total lines number: "
    wc -l $outDir/tmp.vcf
    cat $outDir/tmp.vcf $outDir/tmp${N}.vcf > $outDir/combine.vcf
    mv $outDir/combine.vcf $outDir/tmp.vcf
    wc -l $outDir/tmp.vcf

    rm $outDir/${sample}_chr${N}.recode.vcf
    rm $outDir/tmp${N}.vcf
    echo ">>>>>>>> Done with chr"${N}
done

echo "Done"

mv $outDir/tmp.vcf $outDir/${sample}.with_chr.het.vcf

cat $outDir/${sample}.with_chr.het.vcf | grep -E "^[^#]" > $outDir/${sample}.content.vcf
cat $outDir/${sample}.with_chr.het.vcf | grep "#" > $outDir/${sample}_header.vcf
sed -e 's/chr//' $outDir/${sample}.content.vcf > ${sample}.remove_chr.content.vcf
cat $outDir/${sample}_header.vcf $outDir/${sample}.remove_chr.content.vcf > $outDir/${sample}.het.vcf
cat $outDir/${sample}.het.vcf | bgzip > $outDir/${sample}.het.vcf.gz
cat $outDir/${sample}.with_chr.het.vcf | bgzip > $outDir/${sample}.with_chr.het.vcf.gz
rm *log
rm *header.vcf
rm *content.vcf

#for sample in "HG00097" "HG00250" "HG00353" "NA20801" "HG00353" "NA19247"
#do
#done
