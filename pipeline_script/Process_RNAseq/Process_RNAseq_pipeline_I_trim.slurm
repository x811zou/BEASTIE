#!/bin/bash
#SBATCH --get-user-env
#SBATCH -A reddylab
#SBATCH -p new,all
#SBATCH --nice=500
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL 
#SBATCH --out=1KGP-I-trim-%j.out 
#SBATCH --error=1KGP-I-trim-%j.out

# Please specify below
#fastq_R1=/data/allenlab/scarlett/data/fastq/${sample}/*1.fastq.gz
#fastq_R2=/data/allenlab/scarlett/data/fastq/${sample}/*2.fastq.gz
#vcfDir=/data/allenlab/scarlett/data/VCF/GIAB/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
#missN=10

ref=/data/reddylab/Reference_Data/Genomes/hg19
star_ind=/data/allenlab/scarlett/data/STARIndex
AnnoDir=/data/allenlab/Reference_Data/Genomes/hg19/annotations
pipelineDir=/data/allenlab/scarlett/pipeline/Process_RNAseq
outDir=/data/allenlab/scarlett/result/project_ASE/Process_RNAseq

for sample in "HG00250" "HG00353" "NA19247" "NA20801" "NA12878"
do
    fastq_R1=/data/allenlab/scarlett/data/fastq/${sample}/*1.fastq.gz
    fastq_R2=/data/allenlab/scarlett/data/fastq/${sample}/*2.fastq.gz
    mkdir -p $outDir/$sample
    mkdir -p $outDir/$sample/trimmed_fastq
####################################
    echo "step1.1: trim fastq"
    sh $pipelineDir/step1_trim_fastq.sh $fastq_R1 $fastq_R2 $outDir/$sample/trimmed_fastq $sample
#python /data/allenlab/scarlett/python/ASE_pipeline_SZ.py trim_fastq

###################################
    #echo "step1.2: add chr to VCF"
    #sh step2_add_chr_vcf.sh $vcfDir $outDir/chr_vcf $sample
done
###################################
#echo "step2: alignment"
#python /data/allenlab/scarlett/python/ASE_pipeline_SZ.py align_bam

###################################
#echo "step3: mpileup"
#python /data/allenlab/scarlett/python/ASE_pipeline_SZ.py pileup_reads

###################################
#echo "step4: generate ASE data"
#python /data/allenlab/scarlett/python/ASE_pipeline_SZ.py ase_Data
