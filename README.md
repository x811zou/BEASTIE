# BEASTIE: A bioinformatics method for estimation of gene level ASE
BEASTIE (Bayesian Estimation of Allele Specific Transcription Integrating across Exons) is a software suite for identifying allele-specific-expression (ASE) from regulatory variants from RNA-seq and WGS data.
BEASTIE uses a Bayesian hierarchical model to integrate prior information with read count data and genetic data. Using MCMC (Markov Chain Monte Carlo), BEASTIE efficiently performs posterior inference to estimate effect sizes of ASE. <br>
BEASTIE has been found to be substantially more accurate than other tests based on the binomial distribution.

BEASTIE is free for academic and non-profit use.

## Installation
### Prerequisites
The following are required to install and run BEASTIE directly on your system:
* BEASTIE has been tested on **Linux**. It may or may not work on other UNIX systems.
* [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) must be installed.  This is the command-line interface to the STAN statistical programming language.
* [Python](https://www.python.org/downloads/release/python-360/) version 3.6 or higher is required.
* [htslib](https://www.htslib.org/)
* [picard](https://broadinstitute.github.io/picard/) - put location of .jar file in parameters.cfg
* [samtools](https://github.com/samtools/samtools)
* [STAR](https://github.com/alexdobin/STAR)
* [Trimmomatic](https://github.com/usadellab/Trimmomatic) - put location of .jar file in parameters.cfg
* [vcftools](https://vcftools.github.io/)

### Installing and Compiling BEASTIE source code
First download BEASTIE, copy its files into your working directory.
```python
s = "example code"
```
Then, install [CmdStan](https://mc-stan.org/users/interfaces/cmdstan), and set the environment variable $STAN to the directory where CmdStan has been installed. 
```python
s = "example code"
```

## Workflow
### Summary of steps
Multiple steps are needed to identify gene level ASE. Broadly, these steps are:

* Step1: Gene-level pileup read counts generation. Using STAR 2Pass EndtoEnd alignment mode with WASP filtering for RNAseq fastq data alignment, and extract heterozygous sites from VCF files for samtools mpile up. 
* Step2: BEASTIE model input data preparation. Parsing pileup read counts by using the faster adopted python script from [ASEreadCounter](https://github.com/gimelbrantlab/ASEReadCounter_star). 
* Step3: Identification of genes with ASE by running BEASTIE.

![alt text](workflow_figure/steps.png "steps")

### Summary of workflow
<place holder for workflow figure>
 
  
Before the BEASTIE model can be run, you must create a file containing the read counts for each allele of a gene.  The format of this required file is described below.
```
gene_ID | ALT1 | REF1 | ALT2 | REF2 | pred_prob
```
  
0. input files
  
The model must be run in the $STAN directory.  The following command will run the model on a set of variants:
```
BEASTIE.py A B C D > out.txt
```
The parameters are:
* A
* B
* C
* D

1. XX.py
  
This step XX
```
mkdir 
python /BEASTIE/XX.py 
```
  
2. XX.py
  
This step XX
```
python /BEASTIE/XX.py 
```
  

n. Output
The output consists of XXX columns: ???

