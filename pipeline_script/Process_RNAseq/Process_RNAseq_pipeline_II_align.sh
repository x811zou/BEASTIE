#!/bin/bash
# sh Process_RNAseq_pipeline_II_align.sh 

####### specify the below directory
module load samtools/1.9-gcb01
module load STAR/2.7.2b-gcb01
module load bedtools2/2.25.0-fasrc01

sample=$1
ref=$2
#ref=/data/reddylab/Reference_Data/Genomes/hg19
star_ind=$3
#star_ind=/data/allenlab/scarlett/data/STARIndex
AnnoDir=$4
#AnnoDir=/data/allenlab/Reference_Data/Genomes/hg19/annotations
mismatch_N=$5
pipelineDir=/Process_RNAseq
fastqDir=$pipelineDir/$sample/trimmed_fastq
OutDir=$pipelineDir/$sample
VCF=$pipelineDir/$sample/chr_vcf/${sample}_chr.vcf

N=$mismatch_N
echo "====================================================================="${sample}
fastqDir=$pipelineDir/$sample/trimmed_fastq
OutDir=$pipelineDir/$sample
VCF=$pipelineDir/$sample/chr_vcf/${sample}_chr.vcf
cd $OutDir
folder=mismatch${N}
tool2="star2pass_EndtoEnd_wasp"
OutDir=$pipelineDir/$sample/$tool2
mkdir -p $OutDir
mkdir -p $OutDir/$folder
cd $OutDir/$folder
STAR --twopassMode Basic --runThreadN 24 --genomeDir $star_ind \
--readFilesIn $fastqDir/${sample}_FWD_paired.fq.gz $fastqDir/${sample}_REV_paired.fq.gz \
--alignEndsType EndToEnd --waspOutputMode SAMtag --varVCFfile $VCF --outFilterMismatchNmax $N \
--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
--outSAMattributes NH HI NM MD AS nM jM jI XS vA vG vW --readFilesCommand "gunzip -c" \
--outFileNamePrefix $OutDir/$folder/

echo "Finish STAR 2pass EndtoEnd WASP alignment"
java -jar /data/reddylab/software/picard/dist/picard.jar MarkDuplicates I=$OutDir/$folder/Aligned.sortedByCoord.out.bam O=$OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup.bam M=$OutDir/$folder/picard.marked_dup_metrics.txt
cd $OutDir/$folder

samtools index $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup.bam
cat $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup.bam | samtools view -h | grep -vE "vW:i:2|vW:i:3|vW:i:4|vW:i:5|vW:i:6|vW:i:7" | samtools view -bS -o $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam
echo "Finish filtering"
rm $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup.bam*
echo "Start samtools flagstat for picard"
samtools flagstat $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam > $OutDir/$folder/flagstat_markdup_picard.txt
samtools index $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam
echo "Start counting unique read pair"
samtools view -bh -q 30 -f 3 -F 2316 -c $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam
samtools view -bh -q 30 -f 3 -F 3340 -c $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam

# if a list of samples, then 
#for sample in "NA12878"
#do
#done