#!/bin/bash
#SBATCH --get-user-env
#SBATCH -A reddylab
#SBATCH -p new,all
#SBATCH --nice=500
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=log-GSD-hetsMeta-%j.out
#SBATCH --output=log-GSD-hetsMeta-%j.out

module load python/3.7.4-gcb01


for sample in "123375" "125249" "125260" "122687" "122698" "123667"
do
    vcfDir=/data/allenlab/scarlett/data/VCF/GSD/DNA_vcf/${sample}.vcf.recode.vcf.gz
    PipelineDir=/data/allenlab/scarlett/result/project_ASE/Process_VCF
    mkdir -p $PipelineDir
    outDir=$PipelineDir/$sample
    mkdir -p $outDir
    echo "step1: extract het sites for mpileup"${sample}
    mkdir -p $outDir/hetsMeta
    mkdir -p $outDir/hetsDict
    mkdir -p $outDir/hetsDict/genom_all
    mkdir -p $outDir/hetsDict/trans_all
    python /data/allenlab/scarlett/pipeline/Process_VCF/step1_extracthets.py $sample $vcfDir $outDir  
done

#PipelineDir=/data/allenlab/scarlett/result/project_ASE/Process_VCF
#mkdir -p $PipelineDir
#outDir=$PipelineDir/$sample
#mkdir -p $outDir
######################################
#echo "step1: extract het sites for mpileup"
#mkdir -p $outDir/hetsMeta
#mkdir -p $outDir/hetsDict
#mkdir -p $outDir/hetsDict/genom_all
#mkdir -p $outDir/hetsDict/trans_all
#python step1_extracthets.py $sample $vcfDir $outDir  
#python step1_extracthets.py NA12878 /data/allenlab/scarlett/data/VCF/GIAB/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz /data/allenlab/scarlett/result/project_ASE/Process_VCF/NA12878
echo "step1 done!"
