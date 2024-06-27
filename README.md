# BEASTIE: A bioinformatics method for gene level ASE estimation
BEASTIE (Bayesian Estimation of Allele Specific Transcription Integrating across Exons) is a software suite for identifying allele-specific-expression (ASE) from regulatory variants from RNA-seq and WGS data.
BEASTIE uses a Bayesian hierarchical model to integrate prior information with read count data and genetic data. Using MCMC (Markov Chain Monte Carlo), BEASTIE efficiently performs posterior inference to estimate effect sizes of ASE. <br>

The BEASTIE workflow is currently set up for gene-level ASE estimation. This method has been tested in Europeran individual NA12878, and African individual NA19247 from 1000 Genome Project.

BEASTIE has been found to be substantially more accurate than other tests based on the binomial distribution.

BEASTIE is free for academic and non-profit use.


## Installation options
### Using Docker Image locally:
Docker does not allow the container any access to the host system by default. To allow the BEASTIE container to access necessary
files on the host system you need use the `-v` flag to "mount" a host directory at a particular mount point inside the container.

If all your data and outputs exist under the current directory, the following template should work for running BEASTIE. Just replace "<BEASTIE_IMAGE>" with the docker image and "<config_file>" with the path to your BEASTIE configuration file.

```bash
% docker run -v "`pwd`":"`pwd`" xuezou/beastie -c <config_file>
```

### Using Singularity in cluster:
We don't build a singularity image directly, but you can build one using the docker image. 

1. pull docker image from Dockerhub, and create beastie.tar file
```bash
% docker save --output beastie.tar xuezou/beastie
```

2. copy the beastie.tar file to your cluster using scp

3. Create a singularity image, `beastie.sif` that can potentially be run like
```bash
% singularity build -s beastie.sif docker-archive://beastie.tar
```

4. Run
```bash
% singularity run --bind <directory> beastie.sif -c <config_file>
```


### Customized installation:
#### Software prerequisites
Only if you want to use our recommended pipeline to align RNAseq reads and modify the pipeline:
* [bedtools2.25](https://bedtools.readthedocs.io/en/latest/content/installation.html)
* [picard](https://broadinstitute.github.io/picard/) - set location as $picard_path
* [samtools1.9](https://github.com/samtools/samtools)
* [htslib-1.12](http://www.htslib.org/download/)
* [STAR2.7](https://github.com/alexdobin/STAR)  - set location as $STAR
* [Trimmomatic](https://github.com/usadellab/Trimmomatic) - set location as $trimmomatic_path
* [vcftools0.1.15](https://vcftools.github.io/)
* [tabix](https://github.com/samtools/htslib)

The following tools are required to install and run BEASTIE directly on your system:
* BEASTIE has been tested on **Linux**. It may or may not work on other UNIX systems.
* [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) must be installed.  This is the command-line interface to the STAN statistical programming language. Set the location as $STAN. Installing and Compiling BEASTIE source code.
* [Python 3.6](https://www.python.org/downloads/release/python-360/) version 3.6 or higher is required.
* [R 4.0](https://cran.r-project.org/bin/macosx/)

The following Python packages are required to install in your system:
* os, sys, configparser, subprocess, pandas, re, pickle, numpy

The following R packages are required to install in your system:
* [LDlinkR](https://github.com/CBIIT/LDlinkR)
* [reshape2](https://github.com/hadley/reshape)
* [data.table](https://github.com/Rdatatable/data.table)
* [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html)
* [pasilla](https://bioconductor.org/packages/release/data/experiment/html/pasilla.html)
* [readr](https://cran.r-project.org/web/packages/readr/readme/README.html)
* [glmnetUtils](https://www.rdocumentation.org/packages/glmnetUtils/versions/1.1.8)

Most of these pakages can be installed from CRAN using the `install.packages` R function. However, "pasilla" is a Bioconductor package and must be installed using the
[Bioconductor Manager](https://cran.r-project.org/web/packages/BiocManager/index.html) package. Once BiocManager is installed run `BiocManager::install("pasilla")` to install it.
Git clone our BEASTIE scripts and example data in your working directory ($workdir)
```bash
% git clone --recurse-submodules https://github.com/x811zou/BEASTIE.git
```
Installing [CmdStan](https://mc-stan.org/users/interfaces/cmdstan), and set the environment variable $STAN to the directory where CmdStan has been installed.

Register [LD token](https://ldlink.nci.nih.gov/?tab=apiaccess)

Compiling BEASTIE stan model
```bash
% cd $STAN
% mkdir iBEASTIE2
% mv $workdir/BEASTIE/BEASTIE/iBEASTIE2.stan $STAN/iBEASTIE2/.
% make /iBEASTIE2/iBEASTIE2
% cp /iBEASTIE2/iBEASTIE2 /usr/local/bin # or some other directory in your PATH
```

Download reference data and unzip it, and set the environment variable $refdir to the directory where reference folder has been unzipped.
```
https://drive.google.com/file/d/106HEFKV_0TAzoekXNEWYf-XqdqAYZoh4/view?usp=sharing
```

Then either install BEASTIE:
```bash
% cd $workdir/BEASTIE
% make install
```

or install BEASTIE into a python virtual environment:
```bash
% cd $workdir/BEASTIE
% virtualenv -p python3 venv # create a new virtual environment with at least python 3.8
% ./venv/bin/activate
% make install
```

## Workflow
### Summary of steps
Multiple steps are needed to identify gene level ASE. Broadly, these steps are:

Preparation-step:

Gene-level pileup read counts generation. We recommend using STAR 2Pass EndtoEnd alignment mode with WASP filtering for RNAseq fastq data alignment to generate BAM files. Extract allele frequency information for each heterozygous variant from 1000 Genome VCF file for corresponding ancestry (We provide AF_1_22.tsv in reference folder for all ancestry data).

Pipeline-step:

step1: Model input data preparation.
* Extract heterozygous sites from gencode reference for samtools mpileup (We provide pre-split gencode v19 for all 22 chromosome in reference folder, users are free to use their own version of gencode reference and use vcftools tools to split it).
* Parse pileup read counts by our faster version python script originally adopting from [ASEreadCounter](https://github.com/gimelbrantlab/ASEReadCounter_star).
* Thinning reads by read length. One read only count once.
* Annotate AF and LD for bi-allelic het SNPs pairs

step2: Identification of genes with ASE. Parsing BEASTIE model output with customized significance cutoff.
* Convert data in format for model input
* Predict phasing error
* Update model input with phasing error
* Run BEASTIE model
* Generate gene list with user-defined cutoff


### Summary of workflow

Functionally, these above steps are accomplished by individual Python3 scripts, alongside the prior listed dependencies. This workflow is summarized in the below figure:

![alt text](image/newBEASTIE_workflow.png "workflow")

This workflow is summarized step-by-step below.

Preparation-step: process raw data (optional with provided commands)
----------------------------------------

a. trim raw RNAseq fastq reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
```bash
# install at $trimmomatic_path
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip

% java -jar $trimmomatic_path/trimmomatic-0.33.jar PE -threads 16 -phred33 $fastq_R1 $fastq_R2 \
     $trimmed_fastq/${sample}_FWD_paired.fq.gz $trimmed_fastq/${sample}_FWD_unpaired.fq.gz \
     $trimmed_fastq/${sample}_REV_paired.fq.gz $trimmed_fastq/${sample}_REV_unpaired.fq.gz \
     ILLUMINACLIP:$trimmomatic_reference/trimmomatic_MHPS.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
The parameters are:
* $trimmed_fastq: saving trimmed fastq output path
* $sample: sample name or prefix for output


b. align reads with STAR<br>
We have done extensive comparison on RNAseq alignes reference allele mapping bias, and found that the best one with high efficiency and minimal bias is splice-aware aligner STAR with 2pass EndtoEnd alignment mode and WASP filtering : https://github.com/alexdobin/STAR. If you prefer to use aligned BAM files, you can directly use that as input.
```bash
% STAR --twopassMode Basic --runThreadN 24 --genomeDir $star_ind \
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

% java -jar $picardDir/picard.jar MarkDuplicates \
     I=$output_prefix/Aligned.sortedByCoord.out.bam \
     O=$output_prefix/Aligned.sortedByCoord.out.picard_markdup.bam
```
The parameters are:
* $fastqDir: input path for trimmed fastq
* $star_ind: index for STAR aligner
* $VCF: VCF file for sample. VCF file has 'chr' in the chromosome column.
* $sample: sample name or prefix for output
* $output_prefix: output path with output prefix for aligned BAM


c. pile up reads for each variant.
```bash
% samtools mpileup -d 0 -B -s -f $ref -l $het_sites_for_mpileup $bam > ${prefix}.pileup
```
The parameters are:
* $ref: annotation reference file
* $het_sites_for_mpileup: heterozygous sites extracted from VCF. Format as "chr | position"
```
chr1 | 11111
```
* $bam: Star_aligned_sortedByCoord_picard_markdup_filter.bam
* $prefix: prefix for output
----------------------------------------

Pipeline-step: prepare input and run model (required with example config scripts)
----------------------------------------

a. input files required
* sample.vcf
* sample.pileup

b. run BEASTIE pipeline
Parameters can be specificed in parameters.cfg files.
The model (BEASTIE.stan) must be run in the $STAN directory.
```bash
% beastie -c <config_file>
```
----------------------------------------
