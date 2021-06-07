#https://github.com/PacificBiosciences/pbmm2
#https://anaconda.org/bioconda/pbmm2

conda create -n Python27 python=2.7
conda activate Python27
conda config --add channels bioconda
conda config --add channels bioconda/linux-64
conda install -c bioconda/linux-64 pbmm2
#conda deactivate


refDir=/data/reddylab/Reference_Data/Genomes/hg38/hg38.fa
fastqDir=/data/reddylab/scarlett/GSD/Pacbio/bam2fastq
bamDir=/data/reddylab/scarlett/GSD/Pacbio/pbmm2
pipelineDir=/data/allenlab/scarlett/pipeline/Process_pacbio/2.pbmm2

############ 1. Generate index file for reference and reuse it to align reads
#cd $pipelineDir
#pbmm2_idx=$pipelineDir/indx/ref.mmi
#pbmm2 index $refDir $pbmm2_idx
#echo "finish generating index"

mkdir -p $bamDir/${sample}_GRCh38 
cd $bamDir/$sample

sample="control"
subreads_bam=/data/reddylab/GSD/pacbio_control_CP0416_JP06/Bam/lima.bc1001_BAK8A_OA--bc1001_BAK8A_OA.bam
mkdir -p $bamDir/${sample}_GRCh38

sample="CP-0416"
subreads_bam=/data/reddylab/GSD/pacbio_control_CP0416_JP06/Bam/lima.bc1002_BAK8A_OA--bc1002_BAK8A_OA.bam	
mkdir -p $bamDir/${sample}_GRCh38

sample="JP-06"
subreads_bam=/data/reddylab/GSD/pacbio_control_CP0416_JP06/Bam/lima.bc1003_BAK8A_OA--bc1003_BAK8A_OA.bam
mkdir -p $bamDir/${sample}_GRCh38

sample="MA-1"
subreads_bam=/data/reddylab/GSD/pacbio_MA1/r54089_20201109_134937/m54089_201109_140212.subreads.bam
mkdir -p $bamDir/${sample}_GRCh38
#####################
aligned_subreads_sorted_bam=${bamDir}/${sample}_GRCh38/${sample}.aligned.subreads.sorted.bam
pbmm2 align $refDir $subreads_bam $aligned_subreads_sorted_bam --sort --log-level INFO
     
############ 2. Align CCS fastq input and sort output for MA-1
ccs_fastq=/data/reddylab/GSD/pacbio_MA1/CCS_Files/*.fastq
aligned_ccs_sorted_bam=${bamDir}/${sample}_GRCh38/${sample}.aligned.ccs.sorted.bam
pbmm2 align $refDir $ccs_fastq $aligned_ccs_sorted_bam --preset CCS --sort --rg '@RG\tID:myid\tSM:MA-1' --log-level INFO

conda deactivate
