#!/usr/bin/env python
#=========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
#=========================================================================
import os
from sys import argv
#from sklearn import preprocessing

#============================specify parameters============================
#============================ please customize
prefix="NA12878"     #step2-3-4
#prefix="HG00096"     #step1
#prefix="NA19247"     #step1
min_total_cov = 1
min_single_cov = 0
model="BEASTIE"
STAN="/Users/scarlett/allenlab/software/cmdstan/examples/"
out="./output/"
if not os.path.exists(out):
    os.makedirs(out)
sigma=0.5
cutoff=0.5
alpha=0.05

#============================ do not change below this line please
local_dir="/Users/scarlett/Box/Allen_lab/Backup/BEASTIE3/github_example/"
model=str(STAN)+str(model)+"/"+str(model)
vcfgz_file=local_dir+"BEASTIE_example/"+str(prefix)+".remove_chr.content.SNPs.hets.header.vcf"     #NA12878/HG00096/NA19247
#CHROM   | POS  |    ID     |  REF   |  ALT  |  QUAL |  FILTER |     INFO      |           FORMAT             |                HG001
#1       837214  rs72631888      G       C       50      PASS    platforms=3;   GT:DP:ADALL:AD:GQ:IGT:IPS:PS    0|1:558:122,138:135,166:581:0/1:.:PATMAT
pileup_file=local_dir+"BEASTIE_example/"+str(prefix)+".pileup"                                     #HG00096/NA19247
#chr10	323283	A	1	g	E	~
hetSNP_file=local_dir+"BEASTIE_example/"+str(prefix)+"_hetSnps.tsv"                                #HG00096/NA19247
#chr   | chrN  |     geneID     | genomicCoord_pos | transcriptCoord |    SNP_id |  genotype
#chr1    1       ENSG00000227232.4       14930             -1         rs75454623      1|0     
meta_file=local_dir+"BEASTIE_example/"+str(prefix)+"_logisticRegression_input.tsv"            #NA12878
#"chr"   "pos"   "lag_pos"       "rsid"  "exon_start"    "exon_end"      "geneID"        "error" "distance"      "log10_distance"        "r2"    "d"     "EUR_MAF"       "lag_EUR_MAF"   "min_EUR_MAF"   "diff_EUR_MAF"
#"chr1"  "   916549"     "   914940"     "rs6660139"     "   916516"     "   916553"     "ENSG00000187642"       "0"     "  1609"        "3.2065560"     "0.503" "1.000" "0.2416"        "0.4105"        "0.2416"        "0.1689"

#========================= run scripts for BEASTIE ==========================
#cmd = "python stan_wrapper_2.py %s %s %s %s %s %s %s" % (str(vcfgz_dir),str(pileup_dir),str(hetSNP_dir),min_total_cov,min_single_cov,str(prefix),out,alpha,stanModel)
cmd = "python stan_wrapper.py %s %s %s %s %s %d %d %s %s %s %f %f %f" % (str(local_dir),str(vcfgz_file),str(pileup_file),str(hetSNP_file),str(meta_file),min_total_cov,min_single_cov,str(prefix),str(out),str(model),alpha,sigma,cutoff)
os.system(cmd)

print("BEASTIE finish! Please check your output!")

