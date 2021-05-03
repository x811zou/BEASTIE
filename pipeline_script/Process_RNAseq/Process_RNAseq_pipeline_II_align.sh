#!/bin/bash

sample="${1}"
ref="${2}"
star_ind="${3}"
mismatch_N="${4}"
picard="${5}"
pipelineDir="${6:-./Process_RNAseq}"

fastqDir="${pipelineDir}/${sample}/trimmed_fastq"
VCF="${pipelineDir}/${sample}/chr_vcf/${sample}_chr.vcf"
outDir="${pipelineDir}/${sample}/star2pass_EndtoEnd_wasp/mismatch${mismatch_N}"

echo "====================================================================="

mkdir -p ${outDir}

STAR --twopassMode Basic --runThreadN 24 --genomeDir ${star_ind} \
  --readFilesIn ${fastqDir}/${sample}_FWD_paired.fq.gz ${fastqDir}/${sample}_REV_paired.fq.gz \
  --alignEndsType EndToEnd --waspOutputMode SAMtag --varVCFfile ${VCF} --outFilterMismatchNmax ${mismatch_N} \
  --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
  --outSAMattributes NH HI NM MD AS nM jM jI XS vA vG vW --readFilesCommand "gunzip -c" \
  --outFileNamePrefix ${outDir}/

echo "Finish STAR 2pass EndtoEnd WASP alignment"
java -jar ${picard} MarkDuplicates I=${outDir}/Aligned.sortedByCoord.out.bam \
  O=${outDir}/Aligned.sortedByCoord.out.picard_markdup.bam \
  M=${outDir}/picard.marked_dup_metrics.txt

samtools index ${outDir}/Aligned.sortedByCoord.out.picard_markdup.bam
cat ${outDir}/Aligned.sortedByCoord.out.picard_markdup.bam | \
  samtools view -h | grep -vE "vW:i:2|vW:i:3|vW:i:4|vW:i:5|vW:i:6|vW:i:7" | \
  samtools view -bS -o ${outDir}/Aligned.sortedByCoord.out.picard_markdup_filter.bam

echo "Finish filtering"
rm ${outDir}/Aligned.sortedByCoord.out.picard_markdup.bam*

echo "Start samtools flagstat for picard"
samtools flagstat ${outDir}/Aligned.sortedByCoord.out.picard_markdup_filter.bam > ${outDir}/flagstat_markdup_picard.txt
samtools index ${outDir}/Aligned.sortedByCoord.out.picard_markdup_filter.bam

echo "Start counting unique read pair"
samtools view -bh -q 30 -f 3 -F 2316 -c ${outDir}/Aligned.sortedByCoord.out.picard_markdup_filter.bam
samtools view -bh -q 30 -f 3 -F 3340 -c ${outDir}/Aligned.sortedByCoord.out.picard_markdup_filter.bam

# if a list of samples, then
#for sample in "NA12878"
#do
#done