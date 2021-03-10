#!/bin/bash
#SBATCH --get-user-env
#SBATCH -A reddylab
#SBATCH -p new,all
#SBATCH --nice=500
#SBATCH --mem=300G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL 
#SBATCH --out=1217-pbmm2-%j.out 
#SBATCH --error=1217-pbmm2-%j.error

refDir=/data/allenlab/scarlett/data/reference/hg38/hg38.fa
fastqDir=/data/reddylab/scarlett/GSD/Pacbio/bam2fastq
bamDir=/data/reddylab/scarlett/GSD/Pacbio/pbmm2

cd $bamDir
pbmm2_idx=$bamDir/ref.mmi
pbmm2 index $refDir $pbmm2_idx
echo "finish generating index"

for sample in "control" "CP-0416" "JP-06"
do
	cd $fastqDir
	mkdir -p $bamDir/$sample
	cd $bamDir/$sample
	input_fastq=$fastqDir/${sample}.fastq
	output_bam=${bamDir}/$sample/${sample}_GRCh38.bam
        output_sort_bam=${bamDir}/$sample/${sample}_GRCh38.sort.bam
	pbmm2 align $pbmm2_idx $input_fastq $output_bam
	samtools sort $output_bam > $output_sort_bam
	samtools index $output_sort_bam
done
