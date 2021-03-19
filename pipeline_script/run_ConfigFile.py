###############################################
# usage: python run_ConfigFile.py 
###############################################
from typing import ForwardRef
import ConfigFile
import os

###############################################
# read in parameters defined in parameter.txt 
###############################################
configFile=ConfigFile("parameter.txt")
sample_name=configFile.lookup("Sample")
fastq_path=configFile.lookup("fastqPath")
vcf_file=configFile.lookup("vcfFile")
model_input_path=configFile.lookup("modelInputPath")
model_input=configFile.lookup("modelInput")
sigma=configFile.lookup("Sigma")
model_output_folder=configFile.lookup("modelOutputFolder")
###############################################
# data processing step by step
###############################################
# step0. check data requirement
### check the existence of RNAseq fastq file and VCF file path
### fastq path requirement: fastq path should have files containing "*1.fastq.gz" and "*2.fastq.gz"
R1_fastq=""
R2_fastq=""
VCF=""
if os.path.exists(fastq_path):
    for filename in os.listdir(fastq_path):
        if "R1.fastq" in filename:
            R1_fastq = fastq_path+"/"+filename
        if "R2.fastq" in filename:
            R2_fastq = fastq_path+"/"+filename        
else:
    print("Oops! fastq path not existed. Try again ...")
    break
if not (R1_fastq or R2_fastq):
    print("Oops! fastq files could not be located. Try again ...")
    break

if os.path.isfile(vcf_file):
    VCF = vcf_file
else:
    print("Oops! VCF file not existed. Try again ...")
    break
###############################################
# Process RNAseq & VCF
###############################################
#### step1.1 trimming fastq
cmd = "sh ./Process_RNAseq/Process_RNAseq_pipeline_I_trim.sh %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
os.system(cmd)
#### step1.2 RNAseq fastq file alignment
cmd = "sh ./Process_RNAseq/Process_RNAseq_pipeline_II_align.sh %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
os.system(cmd)
#### step1.3 clean VCF files
cmd = "sh ./Process_VCF/Process_VCF_pipeline_I.extractVCF.sh %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
os.system(cmd)
#### step1.4 extract het sites information and store it in a file for mpileup
cmd = "sh ./Process_VCF/Process_VCF_pipeline_II_hetsMeta.sh %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
os.system(cmd)
#### step1.5 mpileup
cmd = "sh ./Process_RNAseq/Process_RNAseq_pipeline_III_mpileup.sh %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
os.system(cmd)

###############################################
# Phasing
###############################################
#### step2.1 prepare VCF
cmd = "sh ./Phasing/step1_prepareVCF.sh %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
os.system(cmd)
#### step2.2 Phasing
cmd = "sh ./Phasing/step2_phasing.slurm %s" % (ref,star_ind,AnnoDir,pipelineDir,outDir,sample_name)
os.system(cmd)

###############################################
# BEASTIE
###############################################
### (for now) requires user to make the input file in a format that can be used for BEASTIE
#### step3.1 run BEASTIE
cmd = "python ./BEASTIE/wrapper.py %s" % (model_input,sigma,0,"BEASTIE",model_input_path,model_output_folder)
os.system(cmd)

