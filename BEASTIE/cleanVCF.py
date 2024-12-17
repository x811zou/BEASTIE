#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

import os
import subprocess
import sys
import logging
import shutil

def extract_file_name(input_path):
    file_name_with_extension = os.path.basename(input_path)
    file_name = file_name_with_extension.replace(".vcf.gz", "")
    return file_name

def check_and_cleanup(bi_vcfgz, bihet_vcfgz, tmp_dir):
    if os.path.exists(bi_vcfgz) and os.path.exists(bihet_vcfgz):
        bi_vcfgz_lines = subprocess.run(f"zcat {bi_vcfgz} | head -n20 | wc -l", shell=True, capture_output=True, text=True)
        bihet_vcfgz_lines = subprocess.run(f"zcat {bihet_vcfgz} | head -n20 | wc -l", shell=True, capture_output=True, text=True)
        if int(bi_vcfgz_lines.stdout.strip()) > 11 and int(bihet_vcfgz_lines.stdout.strip()) > 11:
            logging.info(f"..... Check: Both {os.path.basename(bi_vcfgz)} and {os.path.basename(bihet_vcfgz)} exist -- skipping VCF file generation")
            return True
        else:
            subprocess.run(f"rm {bihet_vcfgz}*", shell=True)
            subprocess.run(f"rm {bi_vcfgz}*", shell=True)
    else:
        logging.info(f"..... Check: {os.path.basename(bi_vcfgz)} && {os.path.basename(bihet_vcfgz)} does not exist -- generating VCF files")
    if os.path.exists(tmp_dir):
        subprocess.run(f"rm -r {tmp_dir}", shell=True)
    return False

def filter_vcf(input_vcfgz, tmp_dir, sample, vcf_sample_name, bihet_vcfgz, threads, cutoff, pass_flag):
    os.makedirs(tmp_dir, exist_ok=True)
    logging.info("..... Filtering VCF file")

    # Shared command templates
    filter_condition = '$7=="PASS"' if pass_flag else f"$6>={cutoff}"
    het_condition = '($10=="1/0" || $10=="0/1" || $10=="1|0" || $10=="0|1")'

    # Command to filter VCF based on the PASS flag or quality score
    command = f"""
    bgzip -@ {threads} -cd {input_vcfgz} | \
    bcftools view -s '{vcf_sample_name}' | \
    grep -v '^#' | \
    awk -v OFS='\\t' 'sub(/:.*/,"",$9) && sub(/:.*/,"",$10) && length($4)==1 && length($5)==1 && {filter_condition} && {het_condition}' | \
    awk '{{gsub(/\\chr/, ""); print}}' > {tmp_dir}/tmp.content.vcf
    """

    # Command to extract header
    header_command = f"zcat {input_vcfgz} | grep '^#' > {tmp_dir}/tmp.header.vcf"

    # Combine header and filtered content
    combined_command = f"cat {tmp_dir}/tmp.header.vcf {tmp_dir}/tmp.content.vcf > {tmp_dir}/{sample}.no_chr.content.SNPs.hets.vcf"

    # Execute the commands
    subprocess.run(command, shell=True, check=True)
    subprocess.run(header_command, shell=True, check=True)
    subprocess.run(combined_command, shell=True, check=True)
    logging.info("..... Finished filtering")

    # Compress and index the final VCF file
    compress_and_index(f"{tmp_dir}/{sample}.no_chr.content.SNPs.hets.vcf", bihet_vcfgz)

    # Extract all filtered SNPs without het filtering
    all_command = f"""
    bgzip -@ {threads} -cd {input_vcfgz} | \
    bcftools view -s '{vcf_sample_name}' | \
    grep -v '^#' | \
    awk -v OFS='\\t' 'sub(/:.*/,"",$9) && sub(/:.*/,"",$10) && length($4)==1 && length($5)==1 && {filter_condition}' > {tmp_dir}/tmp.content.all.vcf
    """
    subprocess.run(all_command, shell=True, check=True)

def compress_and_index(input_vcf, output_vcfgz):
    compress_command = f"bgzip -@2 -c -f {input_vcf} > {output_vcfgz} && tabix -fp vcf {output_vcfgz}"
    subprocess.run(compress_command, shell=True)

def compress_and_index_vcf(tmp_dir, sample, bi_vcfgz):
    header_command = f"cat {tmp_dir}/tmp.header.vcf {tmp_dir}/tmp.content.all.vcf > {tmp_dir}/{sample}.no_chr.content.SNPs.vcf"
    subprocess.run(header_command, shell=True)
    compress_and_index(f"{tmp_dir}/{sample}.no_chr.content.SNPs.vcf", bi_vcfgz)

# python cleanVCF.py /data2/BEASTIE_example_output/NA12878_chr21/tmp/vcf /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/HG001_GRCh37_GIAB.chr21.vcf.gz /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21.bi.vcf.gz /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21.bihets.vcf.gz 0 true
def cleaning(vcf_sample_name, tmp_dir, input_vcfgz, bi_vcfgz, bihet_vcfgz, cutoff, pass_flag):
    logging.basicConfig(
        format="%(asctime)-15s [%(levelname)s] %(message)s",
        level=logging.INFO,
    )
    # Extract the filename from the path
    filename = os.path.basename(bi_vcfgz)
    # Extract the part before the first dot
    sample = filename.split('.')[0]
    out_dir = os.path.dirname(bi_vcfgz)

    logging.info(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CleanVCF: Start sample {sample}")

    os.makedirs(out_dir, exist_ok=True)

    if check_and_cleanup(bi_vcfgz, bihet_vcfgz, tmp_dir):
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        return

    filter_vcf(input_vcfgz, tmp_dir, sample, vcf_sample_name, bihet_vcfgz, threads=1, cutoff=cutoff, pass_flag=pass_flag)
    
    compress_and_index_vcf(tmp_dir, sample, bi_vcfgz)

    if check_and_cleanup(bi_vcfgz, bihet_vcfgz, tmp_dir):
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
    
    logging.info(f"..... Done sample {sample}")

def main(): 
    if len(sys.argv) != 7:
        print("Usage: cleanVCF.py <tmp_dir> <in_vcf> <bi_vcf> <bihet_vcf> <vcf_sample_coln> <quality_score_min> <pass_flag>")
        sys.exit(1)
    
    tmp_dir = sys.argv[1]
    input_vcfgz = sys.argv[2]
    bi_vcfgz = sys.argv[3]
    bihet_vcfgz = sys.argv[4]
    vcf_sample_name = sys.argv[5]
    quality_score_min = sys.argv[6]
    skip_require_pass = False

    cleaning(vcf_sample_name, tmp_dir, input_vcfgz, bi_vcfgz, bihet_vcfgz, quality_score_min, not skip_require_pass)

if __name__ == "__main__":
    main()
