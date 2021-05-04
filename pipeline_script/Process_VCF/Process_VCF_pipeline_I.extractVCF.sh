#!/bin/bash

#######################
# The following options can be used with the recode command to define an INFO key name to keep in the output file.
# These options can be used multiple times to keep more of the INFO fields. The second option is used to keep all
# INFO values in the original file.
#   --recode-INFO <string>
#   --recode-INFO-all
########################

sample="${1}"
vcfDir="${2}"
pipelineDir="${3:-./Process_VCF}"

outDir="${pipelineDir}/${sample}/VCF"
mkdir -p ${outDir}

headerFile=${outDir}/header.vcf
contentFile=${outDir}/content.vcf
tempContentFile=${outDir}/tmp_content.vcf

N=1
echo ">>>>>>>>>>>>>>>> start with chr${N}"
vcftools --gzvcf ${vcfDir}/ALL.chr${N}.*.vcf.gz --indv ${sample} --out ${outDir}/${sample}_chr${N} --recode-INFO-all --recode
echo "${N}.1. Original file has lines:"
wc -l ${outDir}/${sample}_chr${N}.recode.vcf

grep "#" ${outDir}/${sample}_chr${N}.recode.vcf > ${headerFile}
echo "${N}.2. Chr${N} header line count: "
wc -l ${headerFile}

grep -E "^[^#]" ${outDir}/${sample}_chr${N}.recode.vcf > ${contentFile}
echo "${N}.3. Chr${N} content line count: "
wc -l ${contentFile}

echo "${N}.4. CHECK: Concatenated files line count (should match step ${N}.1 above): "
cat ${headerFile} ${contentFile} | wc -l

rm ${outDir}/${sample}_chr${N}.recode.vcf
echo ">>>>>>>>>>>>>>>> Done with chr${N}"

for N in {2..22}
do
    echo ">>>>>>>>>>>>>>>> start with chr${N}"
    vcftools --gzvcf ${vcfDir}/ALL.chr${N}.*.vcf.gz --indv ${sample} --out ${outDir}/${sample}_chr${N} --recode-INFO-all --recode
    echo "${N}.1. ${sample}_chr${N}.recode.vcf line count:"
    wc -l ${outDir}/${sample}_chr${N}.recode.vcf

    grep -E "^[^#]" ${outDir}/${sample}_chr${N}.recode.vcf > ${tempContentfile}
    echo "${N}.2. Chr${N} content line count:"
    wc -l ${tempContentfile}

    echo "${N}.3. Current total content:"
    wc -l ${contentFile}

    echo "${N}.4. Append chr${N} content to total content:"
    cat ${tempContentfile} >> ${contentFile}
    echo "${N}.4.a CHECK: New total content line count (should be sum of steps ${N}.2 and ${N}.3):"
    wc -l ${contentFile}

    rm ${outDir}/${sample}_chr${N}.recode.vcf
    echo ">>>>>>>>>>>>>>>> Done with chr${N}"
done

cat ${headerFile} ${contentFile} > ${outDir}/${sample}.het.vcf

echo "Create a new file with 'chr' prepended to value in the first column."
sed -e 's/^[^#]/chr&/' ${outDir}/${sample}.het.vcf > ${outDir}/${sample}.with_chr.het.vcf

bgzip ${outDir}/${sample}.het.vcf
bgzip ${outDir}/${sample}.with_chr.het.vcf

rm ${headerFile}
rm ${contentFile}
rm ${tempContentfile}

echo "Done"
