#!/bin/bash
#SBATCH --get-user-env
#SBATCH -A reddylab
#SBATCH -p new,all
#SBATCH --nice=500
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=1201-mpileup-125249.stderr
#SBATCH --output=1201-mpileup-125249.out

module load samtools/1.3.1-gcb01

sample="125249"
N=10

PipelineDir=/data/allenlab/scarlett/result/project_ASE/Process_RNAseq
OutDir=${PipelineDir}/${sample}/star2pss_EndtoEnd_wasp/mismatch10/mpileup
mkdir -p $OutDir
ref=/data/reddylab/Reference_Data/Genomes/hg19/hg19.fa
inputDir=/data/allenlab/scarlett/result/project_ASE/Process_VCF/${sample}/hetsMeta
RNAseqDir=${PipelineDir}/${sample}/star2pss_EndtoEnd_wasp/mismatch${N}

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> We start with sample: ${sample}"
echo ">>>>>>>>>> For mismatch"${N}
cd $OutDir
samtools mpileup -d 0 -B -s -f $ref -l $inputDir/Allchr_hets_all_transcript_mpileup.tsv $RNAseqDir/Aligned.sortedByCoord.out.picard_markdup_filter.bam > $OutDir/Allchr_hets_all_transcript.pileup

