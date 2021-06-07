#!/bin/bash
#SBATCH --get-user-env
#SBATCH -A reddylab
#SBATCH -p new,all
#SBATCH --nice=500
#SBATCH --mem=350G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL 
#SBATCH --out=NA12878-II-alignment-%j.out 
#SBATCH --error=NA12878-II-alignment-%j.out

####### specify the below directory
module load samtools/1.9-gcb01
module load STAR/2.7.2b-gcb01
module load bedtools2/2.25.0-fasrc01

#sample="NA12878"

ref=/data/reddylab/Reference_Data/Genomes/hg19
star_ind=/data/allenlab/scarlett/data/STARIndex
AnnoDir=/data/allenlab/Reference_Data/Genomes/hg19/annotations
pipelineDir=/data/allenlab/scarlett/result/project_ASE/Process_RNAseq
#fastqDir=$pipelineDir/$sample/trimmed_fastq
#OutDir=$pipelineDir/$sample
#VCF=$pipelineDir/$sample/chr_vcf/${sample}_chr.vcf
N=10

for sample in "NA12878"
do
    echo "====================================================================="${sample}
    #for N in 10 5 0
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

done


#echo "====================================================================="${sample}
#echo "=================================================== For mismatch"${N}
#folder=mismatch${N}

########### STAR WASP
#cd $OutDir
#mkdir -p $OutDir/$folder
#cd $OutDir/$folder
#
#tool2="star_2pass_WASP"
#cd $OutDir/$folder
#STAR --twopassMode Basic --runThreadN 24 --genomeDir $star_ind \
#--readFilesIn $fastqDir/${sample}_FWD_paired.fq.gz $fastqDir/${sample}_REV_paired.fq.gz \
#--alignEndsType EndToEnd --waspOutputMode SAMtag --varVCFfile $VCF --outFilterMismatchNmax $N \
#--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx \
#--outSAMattributes NH HI NM MD AS nM jM jI XS vA vG vW --readFilesCommand "gunzip -c" \
#--outFileNamePrefix $OutDir/$folder/
#
#echo "Finish STAR 2pass EndtoEnd WASP alignment"
#java -jar /data/reddylab/software/picard/dist/picard.jar MarkDuplicates I=$OutDir/$folder/Aligned.sortedByCoord.out.bam O=$OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup.bam M=$OutDir/$folder/picard.marked_dup_metrics.txt
#cd $OutDir/$folder
## now = $(date +"%T")
## echo "1. Start converting BED: "${now}
# bedtools bamtobed -i Aligned.sortedByCoord.out.bam > Aligned.sortedByCoord.out.bed
# now = $(date +"%T")
# echo "2. Start sorting BED: "${now}
# sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 Aligned.sortedByCoord.out.bed > Aligned.sortedByCoord.out.sorted.bed
# rm Aligned.sortedByCoord.out.bed
#samtools flagstat Aligned.sortedByCoord.out.picard_markdup.bam > flagstat_markdup_picard.txt
#samtools flagstat Aligned.sortedByCoord.out.bam > flagstat_nomarkdup.txt
#samtools index $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup.bam
#cat $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup.bam | samtools view -h | grep -vE "vW:i:2|vW:i:3|vW:i:4|vW:i:5|vW:i:6|vW:i:7" | samtools view -bS -o $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam
#echo "Finish filtering"
#rm $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup.bam*
#echo "Start samtools flagstat for picard"
#samtools flagstat $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam > $OutDir/$folder/flagstat_markdup_picard.txt
#samtools index $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam
#echo "Start counting unique read pair"
#samtools view -bh -q 30 -f 3 -F 2316 -c $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam
#samtools view -bh -q 30 -f 3 -F 3340 -c $OutDir/$folder/Aligned.sortedByCoord.out.picard_markdup_filter.bam

#    done
#done
