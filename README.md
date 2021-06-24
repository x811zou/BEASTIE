# BEASTIE: A bioinformatics method for gene level ASE estimation
BEASTIE (Bayesian Estimation of Allele Specific Transcription Integrating across Exons) is a software suite for identifying allele-specific-expression (ASE) from regulatory variants from RNA-seq and WGS data.
BEASTIE uses a Bayesian hierarchical model to integrate prior information with read count data and genetic data. Using MCMC (Markov Chain Monte Carlo), BEASTIE efficiently performs posterior inference to estimate effect sizes of ASE. <br>

The BEASTIE workflow is currently set up for gene-level ASE estimation. This method has been tested in Europeran individual NA12878, and African individual NA19247 from 1000 Genome Project. 

BEASTIE has been found to be substantially more accurate than other tests based on the binomial distribution.

BEASTIE is free for academic and non-profit use.

## Installation
### Prerequisites
The following are required to install and run BEASTIE directly on your system:
* BEASTIE has been tested on **Linux**. It may or may not work on other UNIX systems.
* [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) must be installed.  This is the command-line interface to the STAN statistical programming language.
* [Python](https://www.python.org/downloads/release/python-360/) version 3.6 or higher is required.
* [htslib-1.12](http://www.htslib.org/download/)
* [bedtools2.25](https://bedtools.readthedocs.io/en/latest/content/installation.html)
* [picard](https://broadinstitute.github.io/picard/) - set location as $picard_path
* [samtools1.9](https://github.com/samtools/samtools)
* [STAR2.7](https://github.com/alexdobin/STAR)
* [Trimmomatic](https://github.com/usadellab/Trimmomatic) - set location as $trimmomatic_path
* [vcftools0.1.15](https://vcftools.github.io/)

### Installation steps

Using a Python 3.8 VirtualEnv:
```python
s = "example code"
```
Using Singularity with Docker Image:
```
```
Customized installation

Installing and Compiling BEASTIE source code
First download BEASTIE, copy its files into your working directory.
Then, install [CmdStan](https://mc-stan.org/users/interfaces/cmdstan), and set the environment variable $STAN to the directory where CmdStan has been installed. 
```bash
cd $STAN
mkdir BEASTIE #download BEASTIE.stan in this directory
make $STAN/BEASTIE/BEASTIE
```

## Workflow
### Summary of steps
Multiple steps are needed to identify gene level ASE. Broadly, these steps are:

* Pre-step: Gene-level pileup read counts generation. We recommend using STAR 2Pass EndtoEnd alignment mode with WASP filtering for RNAseq fastq data alignment, and extract heterozygous sites from VCF files for samtools mpile up. Extract allele frequency information for each heterozygous variant, and obtain linkage equlibirum information for each pair of consecutive heterozygous variants. 
* Pipeline: (i) BEASTIE model input data preparation. Parsing pileup read counts by our faster version python script originally adopting from [ASEreadCounter](https://github.com/gimelbrantlab/ASEReadCounter_star). (ii) Identification of genes with ASE. Parsing BEASTIE model output with customized significance cutoff.

![alt text](image/step.png "steps")

### Summary of workflow

Functionally, these above steps are accomplished by individual Python3 scripts, alongside the prior listed dependencies. This workflow is summarized in the below figure:

![alt text](image/workflow_V4.png "workflow")
 
This workflow is summarized step-by-step below. 
  
----------------------------------------
0. process raw data (optional pre-step with provided commands)
(a) processes trim raw RNAseq fastq reads
```
java -jar $trimmomatic_path/trimmomatic-0.33.jar PE -threads 16 -phred33 $fastq_R1 $fastq_R2 \
   $trimmed_fastq/${sample}_FWD_paired.fq.gz $trimmed_fastq/${sample}_FWD_unpaired.fq.gz \
   $trimmed_fastq/${sample}_REV_paired.fq.gz $trimmed_fastq/${sample}_REV_unpaired.fq.gz \
   ILLUMINACLIP:$trimmomatic_reference/trimmomatic_MHPS.fa:2:30:10:8:TRUE LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:36
```
The parameters are:
* $trimmed_fastq: saving trimmed fastq output path
* $sample: sample name or prefix for output


(b) align reads with STAR<br>
We have done extensive comparison on RNAseq alignes reference allele mapping bias, and found that the best one with high efficiency and minimal bias is splice-aware aligner STAR with 2pass EndtoEnd alignment mode and WASP filtering : https://github.com/alexdobin/STAR. If you prefer to use aligned BAM files, you can directly use that as input. 
```
STAR --twopassMode Basic --runThreadN 24 --genomeDir $star_ind \
    --readFilesIn $fastqDir/${sample}_FWD_paired.fq.gz $fastqDir/${sample}_REV_paired.fq.gz \
    --alignEndsType EndToEnd \
    --waspOutputMode SAMtag \
    --varVCFfile $VCF \
    --outFilterMismatchNmax 10 \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --outSAMattributes NH HI NM MD AS nM jM jI XS vA vG vW \
    --readFilesCommand "gunzip -c" \
    --outFileNamePrefix $output_prefix
    
java -jar $picardDir/picard.jar MarkDuplicates \
   I=$output_prefix/Aligned.sortedByCoord.out.bam \
   O=$output_prefix/Aligned.sortedByCoord.out.picard_markdup.bam
```
The parameters are:
* $fastqDir: input path for trimmed fastq
* $star_ind: index for STAR aligner
* $VCF: VCF file for sample. VCF file has 'chr' in the chromosome column.
* $sample: sample name or prefix for output
* $output_prefix: output path with output prefix for aligned BAM


(c) pile up reads for each variant.
```
samtools mpileup -d 0 -B -s -f $ref -l $het_sites_for_mpileup $bam > ${prefix}.pileup
```
The parameters are:
* $ref: annotation reference file
* $het_sites_for_mpileup: heterozygous sites extracted from VCF. Format as "chr | position"
```
chr1 | 11111
```
* $bam: Star_aligned_sortedByCoord_picard_markdup_filter.bam
* $prefix: prefix for output


4. self-annotation step to prepare for logistic regression (step1)
```
chr | pos | rsid | min_EUR_AF | diff_min_AF | log10_distance | r2 | d |
```
* Exon information: exon start and exon end information for each gene, which is provided in ./BEASTIE_example. Users could use it to filter exonic sites. 
* AF information: Allele frequency for sites among people with different ancestries extracted from 1000GP VCF, which is provided in ./BEASTIE_example. Users could also extract customized AF information from $vcf.
* LD information: users could use R package XX to annotate d,r2 for each pair of consecutive hets sites. 

----------------------------------------
1. input files required

* sample.remove_chr.content.SNPs.hets.header.vcf
* sample.pileup
* sample_logisticRegression_input.tsv

----------------------------------------
2. run BEASTIE pipeline

Parameters can be specificed in parameters.cfg files.
The model (BEASTIE.stan) must be run in the $STAN directory.
```
python run_config.py
```

