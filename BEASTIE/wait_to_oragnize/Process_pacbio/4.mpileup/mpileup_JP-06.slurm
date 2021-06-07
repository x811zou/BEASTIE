#!/bin/bash
#SBATCH --get-user-env
#SBATCH -A reddylab
#SBATCH -p new,all
#SBATCH --nice=500
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=1214-mpileup-JP-06.out
#SBATCH --output=1214-mpileup-JP-06.out

module load samtools/1.3.1-gcb01

#for sample in "MA-1" ""
sample="JP-06"
sample2="125249"
RNAseqDir=/data/allenlab/scarlett/result/project_ASE/Process_Pacbio/pbmm2
OutDir=/data/allenlab/scarlett/result/project_ASE/Process_Pacbio/mpileup/${sample}
mkdir -p $OutDir
ref=/data/reddylab/Reference_Data/Genomes/hg19/hg19.fa
inputDir=/data/allenlab/scarlett/result/hetsMeta_GSD/$sample2/for_mpileup
inputDir2=/data/allenlab/scarlett/result/project_ASE/Process_VCF/$sample2/hetsMeta/for_pileup

#PipelineDir=/data/allenlab/scarlett/result/project_ASE/Process_RNAseq
#OutDir=/data/reddylab/scarlett/GSD/Pacbio/mpileup/${sample}
#mkdir -p $OutDir
#ref=/data/reddylab/Reference_Data/Genomes/hg19/hg19.fa
#inputDir=/data/allenlab/scarlett/result/project_ASE/Process_VCF/125260/hetsMeta
#inputDir2=
#RNAseqDir=/data/reddylab/scarlett/GSD/Pacbio/pbmm2/${sample}

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> We start with sample: ${sample}"
echo ">>>>>>>>>> For mismatch"${N}
cd $OutDir
#cat $inputDir/Allchr_hets_all_transcript.tsv | awk -F';' '{print $1,$7}' > $inputDir/Allchr_hets_all_transcript_mpileup.tsv
# remove header
samtools mpileup -d 0 -B -s -f $ref -l $inputDir2/Allchr_hets_all_transcript_mpileup.tsv $RNAseqDir/$sample/${sample}.aligned.subreads.sorted.bam > $OutDir/for_pileup/Allchr_hets_all_transcript.pileup
echo "target gene region:"

# 125260 : chr16: PHKG2
echo ">> We start with sample: ${sample}"
samtools mpileup -d 0 -B -s -f $ref -r "chr16:30759591-30772490" $RNAseqDir/$sample/${sample}.aligned.subreads.sorted.bam > $OutDir/for_pileup/${sample}_PHKG2.pileup
samtools mpileup -d 0 -B -s -f $ref -l $OutDir/for_check/125249_PHKG2.tsv $RNAseqDir/$sample/${sample}.aligned.subreads.sorted.bam > $OutDir/for_check/${sample}_hets_PHKG2.pileup
echo "done"
