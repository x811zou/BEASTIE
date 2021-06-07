#!/bin/bash
#SBATCH --get-user-env
#SBATCH -A reddylab
#SBATCH -p new,all
#SBATCH --nice=500
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=xz195@duke.edu
#SBATCH --mail-type=END,FAIL 
#SBATCH --error=stringtie_MA1_hg38_minimap2.0126.out
#SBATCH --output=stringtie_MA1_hg38_minimap2.0126.out

module load python/3.7.4-gcb01
module load samtools/1.9-gcb01
module load StringTie/2.1.1-gcb01

## generate for all GSD cases
ref=/data/reddylab/Reference_Data/Genomes/hg38
#ref_ind=/data/allenlab/scarlett/data/STARIndex
AnnoDir=/data/allenlab/scarlett/data/reference/hg38/gencode.v24.annotation.gtf

############################################# short read data: start 2 pass
alignment="minimap2"
stringTieDir=/data/reddylab/scarlett/GSD/Pacbio/stringTie/${alignment}
bamDir=/data/reddylab/scarlett/GSD/Pacbio/${alignment}

mkdir -p $stringTieDir
cd $stringTieDir
rm -r *

echo ">>>>>>>>>>>>>>>>>>>>>>>>> pacbio"
for sample in "MA-1" "CP-0416" "JP-06" "control"
do
    echo "step1 ref_guided >>>>>>>>>>>>>>>>>>>>>>>>> "${sample}
    ### step1: stringTie reference-guided mode
    mkdir -p $stringTieDir/ref_guided/
    mkdir -p $stringTieDir/ref_guided/$sample
    cd $stringTieDir/ref_guided/
    stringtie -p 8 -L -t -G ${AnnoDir} -o $stringTieDir/ref_guided/$sample/transcripts.gtf -l $sample $bamDir/${sample}.minimap2.hg38.sorted.bam
    ### step2: de novo mode
    echo "step2 de novo >>>>>>>>>>>>>>>>>>>>>>>>> "${sample}
    mkdir -p $stringTieDir/de_novo/
    mkdir -p $stringTieDir/de_novo/$sample
    cd $stringTieDir/de_novo/
    stringtie -p 8 -L -t -o $stringTieDir/de_novo/$sample/transcripts.gtf -l $sample $bamDir/${sample}.minimap2.hg38.sorted.bam
    echo "done!"
done
    ### step3.1: ref_guided Merge transcripts from all samples; examine how the transcripts compare with the reference annotation
echo "step3.1 ref_guided merge >>>>>>>>>>>>>>>>>>>>>>>>> "
cd $stringTieDir/ref_guided/
ls -1 */transcripts.gtf > assembly_GTF_list.txt
stringtie --merge -p 8 -G ${AnnoDir} -o $stringTieDir/ref_guided/stringtie_merged.gtf assembly_GTF_list.txt
/data/allenlab/scarlett/software/gffcompare/gffcompare -r $AnnoDir -o $stringTieDir/ref_guided/gffcompare $stringTieDir/ref_guided/stringtie_merged.gtf

    ### step3.2: de novo Merge transcripts from all samples; examine how the transcripts compare with the reference annotation
echo "step3.2 de novo merge >>>>>>>>>>>>>>>>>>>>>>>>> "
cd $stringTieDir/de_novo/
ls -1 */transcripts.gtf > assembly_GTF_list.txt
stringtie --merge -p 8 -G ${AnnoDir} -o $stringTieDir/de_novo/stringtie_merged.gtf assembly_GTF_list.txt
/data/allenlab/scarlett/software/gffcompare/gffcompare -r $AnnoDir -o $stringTieDir/de_novo/gffcompare $stringTieDir/de_novo/stringtie_merged.gtf

    ### step4: estimate transcript abundance and create table counts for Ballgown:
cd $stringTieDir/
mkdir ref_guided_merged
cd ref_guided_merged

for sample in "CP-0416" "JP-06" "control" "MA-1"
do
	mkdir -p $stringTieDir/ref_guided_merged/$sample
    echo "step4 estiamte abudance ref_guided merged >>>>>>>>>>>>>>>>>>>>>>>>> "${sample}
    ### step1: stringTie
    stringtie -L -e -B -p 8 -G $stringTieDir/ref_guided/stringtie_merged.gtf -o $stringTieDir/ref_guided_merged/${sample}/transcripts.gtf $bamDir/${sample}.minimap2.hg38.sorted.bam
    echo "done!"
done

mkdir -p /data/allenlab/scarlett/data/ballgown
mkdir -p /data/allenlab/scarlett/data/ballgown/GSD/
mkdir -p /data/allenlab/scarlett/data/ballgown/GSD/ref_guided_merged/minimap2/
mkdir -p /data/allenlab/scarlett/data/ballgown/GSD/ref_guided_merged/pacbio/minimap2/
cd /data/allenlab/scarlett/data/ballgown/GSD/ref_guided_merged/pacbio/minimap2/
printf "\"ids\",\"type\",\"path\"\n\"CP-0416\",\"GSD_9c\",\"$stringTieDir/ref_guided_merged/CP-0416\"\n\"control\",\"control\",\"$stringTieDir/ref_guided_merged/control\"\n" > CP-0416_vs_control.csv
printf "\"ids\",\"type\",\"path\"\n\"JP-06\",\"GSD_9c\",\"$stringTieDir/ref_guided_merged/JP-06\"\n\"control\",\"control\",\"$stringTieDir/ref_guided_merged/control\"\n" > JP-06_vs_control.csv
printf "\"ids\",\"type\",\"path\"\n\"MA-1\",\"GSD_9c\",\"$stringTieDir/ref_guided_merged/MA-1\"\n\"control\",\"control\",\"$stringTieDir/ref_guided_merged/control\"\n" > MA-1_vs_control.csv

cd $stringTieDir/
mkdir de_novo_merged
cd de_novo_merged
for sample in "CP-0416" "JP-06" "control" "MA-1"
do
    echo "step4 estiamte abudance do novo merged >>>>>>>>>>>>>>>>>>>>>>>>> "${sample}
    mkdir -p $stringTieDir/de_novo_merged/$sample
    stringtie -p 8 -L -G $stringTieDir/de_novo/stringtie_merged.gtf -e -B -o $stringTieDir/de_novo_merged/${sample}/transcripts.gtf $bamDir/${sample}.minimap2.hg38.sorted.bam
    echo "done!"
done
# Run Ballgown using the reference guided, merged transcripts
mkdir -p /data/allenlab/scarlett/data/ballgown/GSD/de_novo_merged/pacbio/minimap2/
cd /data/allenlab/scarlett/data/ballgown/GSD/de_novo_merged/pacbio/minimap2/
printf "\"ids\",\"type\",\"path\"\n\"CP-0416\",\"GSD_9c\",\"$stringTieDir/de_novo_merged/CP-0416\"\n\"control\",\"control\",\"$stringTieDir/de_novo_merged/control\"\n" > CP-0416_vs_control.csv
printf "\"ids\",\"type\",\"path\"\n\"JP-06\",\"GSD_9c\",\"$stringTieDir/de_novo_merged/JP-06\"\n\"control\",\"control\",\"$stringTieDir/de_novo_merged/control\"\n" > JP-06_vs_control.csv
printf "\"ids\",\"type\",\"path\"\n\"MA-1\",\"GSD_9c\",\"$stringTieDir/de_novo_merged/MA-1\"\n\"control\",\"control\",\"$stringTieDir/de_novo_merged/control\"\n" > MA-1_vs_control.csv
