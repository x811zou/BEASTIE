#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

import os
import subprocess
import sys
from datetime import datetime

def extract_file_name(input_path):
    # Extract the file name from the path
    file_name_with_extension = os.path.basename(input_path)
    # Remove the .vcf.gz extension
    file_name = file_name_with_extension.replace(".vcf.gz", "")
    return file_name

def check_and_cleanup(bi_vcfgz, bihet_vcfgz):
    if os.path.exists(bi_vcfgz) and os.path.exists(bihet_vcfgz):
        bi_vcfgz_lines = subprocess.run(f"zcat {bi_vcfgz} | head -n20 | wc -l", shell=True, capture_output=True, text=True)
        bihet_vcfgz_lines = subprocess.run(f"zcat {bihet_vcfgz} | head -n20 | wc -l", shell=True, capture_output=True, text=True)
        if int(bi_vcfgz_lines.stdout.strip()) > 11 and int(bihet_vcfgz_lines.stdout.strip()) > 11:
            now = datetime.now().strftime("%T")
            print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Both {bi_vcfgz} and {bihet_vcfgz} exist: {now}")
            print("both output files exist- skipping")
            sys.exit()
        else:
            subprocess.run(f"rm {bihet_vcfgz}*", shell=True)
            subprocess.run(f"rm {bi_vcfgz}*", shell=True)

def filter_vcf(input_vcfgz, tmp_dir, sample, bihet_vcfgz, threads, pass_flag):
    os.makedirs(tmp_dir, exist_ok=True)
    if pass_flag:
        print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VCF file has PASS, filter by PASS flag in 7th column")
        command = f"bgzip -@ {threads} -cd {input_vcfgz} | grep -v '^#' | awk -v OFS='\t' 'sub(/:.*/,\"\",$10) && length($4)==1 && length($5)==1 && $7==\"PASS\" && ($10==\"1/0\" || $10==\"0/1\" || $10==\"1|0\" || $10==\"0|1\")' | awk '{{gsub(/\\chr/, \"\")}}1' > {tmp_dir}/tmp.content.vcf"
        header_command = f"zcat {input_vcfgz} | grep '^#' > {tmp_dir}/tmp.header.vcf"
        combined_command = f"cat {tmp_dir}/tmp.header.vcf {tmp_dir}/tmp.content.vcf > {tmp_dir}/{sample}.no_chr.content.SNPs.hets.vcf"
    else:
        print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VCF file does not have PASS, filter by quality score in 6th column > 10")
        command = f"bgzip -@ {threads} -cd {input_vcfgz} | grep -v '^#' | awk -v OFS='\t' 'sub(/:.*/,\"\",$10) && length($4)==1 && length($5)==1 && $6>10 && ($10==\"1/0\" || $10==\"0/1\" || $10==\"1|0\" || $10==\"0|1\")' | awk '{{gsub(/\\chr/, \"\")}}1' > {tmp_dir}/tmp.content.vcf"
        header_command = f"zcat {input_vcfgz} | grep '^#' > {tmp_dir}/tmp.header.vcf"
        combined_command = f"cat {tmp_dir}/tmp.header.vcf {tmp_dir}/tmp.content.vcf > {tmp_dir}/{sample}.no_chr.content.SNPs.hets.vcf"

    subprocess.run(command, shell=True)
    subprocess.run(header_command, shell=True)
    subprocess.run(combined_command, shell=True)
    compress_and_index(f"{tmp_dir}/{sample}.no_chr.content.SNPs.hets.vcf", bihet_vcfgz)

    if pass_flag:
        all_command = f"bgzip -@ {threads} -cd {input_vcfgz} | grep -v '^#' | awk -v OFS='\t' 'sub(/:.*/,\"\",$10) && length($4)==1 && length($5)==1 && $7==\"PASS\"' > {tmp_dir}/tmp.content.all.vcf"
    else:
        all_command = f"bgzip -@ {threads} -cd {input_vcfgz} | grep -v '^#' | awk -v OFS='\t' 'sub(/:.*/,\"\",$10) && length($4)==1 && length($5)==1 && $6>10' > {tmp_dir}/tmp.content.all.vcf"

    subprocess.run(all_command, shell=True)

def compress_and_index(input_vcf, output_vcfgz):
    compress_command = f"bgzip -@2 -c -f {input_vcf} > {output_vcfgz} && tabix -fp vcf {output_vcfgz}"
    subprocess.run(compress_command, shell=True)

def compress_and_index_vcf(tmp_dir, sample, bi_vcfgz):
    header_command = f"cat {tmp_dir}/tmp.header.vcf {tmp_dir}/tmp.content.all.vcf > {tmp_dir}/{sample}.no_chr.content.SNPs.vcf"
    subprocess.run(header_command, shell=True)
    compress_and_index(f"{tmp_dir}/{sample}.no_chr.content.SNPs.vcf", bi_vcfgz)

# python cleanVCF.py NA12878_chr21 /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/HG001_GRCh37_GIAB.chr21.vcf.gz /data2/BEASTIE_example_output/NA12878_chr21 HG001_GRCh37_GIAB.chr21.bi.vcf.gz HG001_GRCh37_GIAB.chr21.bihets.vcf.gz true

def main():
    vcf_sample = sys.argv[1]
    input_vcfgz, out_dir, bi_vcfgz, bihet_vcfgz, pass_flag = sys.argv[2:7]

    now = datetime.now().strftime("%T")
    print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with {vcf_sample} : {now}")

    # Extract file name without path and .vcf.gz extension
    file_name = extract_file_name(input_vcfgz)
    os.makedirs(out_dir, exist_ok=True)
    bi_vcfgz = f"{out_dir}/{file_name}.bi.vcf.gz"
    bihet_vcfgz = f"{out_dir}/{file_name}.bihets.vcf.gz"

    # Check if the output files already exist
    check_and_cleanup(bi_vcfgz, bihet_vcfgz)

    # if output files do not exist
    tmp_dir = os.path.join(out_dir, "tmp")
    filter_vcf(input_vcfgz, tmp_dir, vcf_sample, bihet_vcfgz, threads=1, pass_flag=True)
    
    compress_and_index_vcf(tmp_dir, vcf_sample, bi_vcfgz)

    # check the existence and non-emptyness of bihet_vcfgz, bi_vcfgz file, and then remove the tmp directory
    check_and_cleanup(bi_vcfgz, bihet_vcfgz)
    subprocess.run(f"rm -r {tmp_dir}", shell=True)

    now = datetime.now().strftime("%T")
    print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> done with {vcf_sample} : {now}")


if __name__ == "__main__":
    main()