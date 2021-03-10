#!/bin/bash
#SBATCH --get-user-env
#SBATCH -A reddylab
#SBATCH -p new,all
#SBATCH --nice=500
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=1211-mpileup-MA-1.out
#SBATCH --output=1211-mpileup-MA-1.out

module load samtools/1.3.1-gcb01

#for sample in "MA-1" ""
sample="MA-1"

RNAseqDir=/data/allenlab/scarlett/result/project_ASE/Process_Pacbio/pbmm2
OutDir=/data/allenlab/scarlett/result/project_ASE/Process_Pacbio/mpileup/${sample}
mkdir -p $OutDir
ref=/data/reddylab/Reference_Data/Genomes/hg19/hg19.fa
inputDir=/data/allenlab/scarlett/result/hetsMeta_GSD/123375/for_mpileup

#RNAseqDir=/data/reddylab/scarlett/GSD/Pacbio/pbmm2/${sample}

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> We start with sample: ${sample}"
echo ">>>>>>>>>> For mismatch"${N}
cd $OutDir
samtools mpileup -d 0 -B -s -f $ref -l $inputDir/Allchr_hets_all_transcript.tsv $RNAseqDir/$sample/MA-1.aligned.ccs.sorted.bam > $OutDir/for_pileup/Allchr_hets_all_transcript.pileup
echo "target gene region:"

# 123375 : chr1 : AGL
#sample=123375
echo ">> We start with sample: ${sample}"
samtools mpileup -d 0 -B -s -f $ref -r "chr1:100315640-100389579" $RNAseqDir/$sample/MA-1.aligned.ccs.sorted.bam > $OutDir/for_pileup/${sample}_all_AGL.pileup
samtools mpileup -d 0 -B -s -f $ref -l $OutDir/for_check/123375_AGL.tsv $RNAseqDir/$sample/MA-1.aligned.ccs.sorted.bam > $OutDir/for_check/${sample}_hets_AGL.pileup
echo "done"
