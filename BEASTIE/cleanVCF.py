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

def filter_vcf(input_vcfgz, tmp_dir, sample, bihet_vcfgz, threads, cutoff, pass_flag):
    os.makedirs(tmp_dir, exist_ok=True)
    if pass_flag:
        logging.info(f"..... VCF file has PASS, filter by PASS flag in 7th column")
        command = f"bgzip -@ {threads} -cd {input_vcfgz} | grep -v '^#' | awk -v OFS='\t' 'sub(/:.*/,\"\",$9) && sub(/:.*/,\"\",$10) && length($4)==1 && length($5)==1 && $7==\"PASS\" && ($10==\"1/0\" || $10==\"0/1\" || $10==\"1|0\" || $10==\"0|1\")' | awk '{{gsub(/\\chr/, \"\")}}1' > {tmp_dir}/tmp.content.vcf"
    else:
        logging.info(f"..... VCF file does not have PASS, filter by quality score in 6th column >= {cutoff}")
        command = f"bgzip -@ {threads} -cd {input_vcfgz} | grep -v '^#' | awk -v OFS='\t' 'sub(/:.*/,\"\",$9) && sub(/:.*/,\"\",$10) && length($4)==1 && length($5)==1 && $6>={cutoff} && ($10==\"1/0\" || $10==\"0/1\" || $10==\"1|0\" || $10==\"0|1\")' | awk '{{gsub(/\\chr/, \"\")}}1' > {tmp_dir}/tmp.content.vcf"

    header_command = f"zcat {input_vcfgz} | grep '^#' > {tmp_dir}/tmp.header.vcf"
    combined_command = f"cat {tmp_dir}/tmp.header.vcf {tmp_dir}/tmp.content.vcf > {tmp_dir}/{sample}.no_chr.content.SNPs.hets.vcf"

    subprocess.run(command, shell=True)
    subprocess.run(header_command, shell=True)
    subprocess.run(combined_command, shell=True)
    logging.info(f"..... Finished filtering")
    compress_and_index(f"{tmp_dir}/{sample}.no_chr.content.SNPs.hets.vcf", bihet_vcfgz)

    if pass_flag:
        all_command = f"bgzip -@ {threads} -cd {input_vcfgz} | grep -v '^#' | awk -v OFS='\t' 'sub(/:.*/,\"\",$10) && length($4)==1 && length($5)==1 && $7==\"PASS\"' > {tmp_dir}/tmp.content.all.vcf"
    else:
        all_command = f"bgzip -@ {threads} -cd {input_vcfgz} | grep -v '^#' | awk -v OFS='\t' 'sub(/:.*/,\"\",$10) && length($4)==1 && length($5)==1 && $6>={cutoff}' > {tmp_dir}/tmp.content.all.vcf"

    subprocess.run(all_command, shell=True)

def compress_and_index(input_vcf, output_vcfgz):
    compress_command = f"bgzip -@2 -c -f {input_vcf} > {output_vcfgz} && tabix -fp vcf {output_vcfgz}"
    subprocess.run(compress_command, shell=True)

def compress_and_index_vcf(tmp_dir, sample, bi_vcfgz):
    header_command = f"cat {tmp_dir}/tmp.header.vcf {tmp_dir}/tmp.content.all.vcf > {tmp_dir}/{sample}.no_chr.content.SNPs.vcf"
    subprocess.run(header_command, shell=True)
    compress_and_index(f"{tmp_dir}/{sample}.no_chr.content.SNPs.vcf", bi_vcfgz)

# python cleanVCF.py /data2/BEASTIE_example_output/NA12878_chr21/tmp/vcf /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/HG001_GRCh37_GIAB.chr21.vcf.gz /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21.bi.vcf.gz /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21.bihets.vcf.gz 0 true
def cleaning(tmp_dir, input_vcfgz, bi_vcfgz, bihet_vcfgz, cutoff, pass_flag):
    logging.basicConfig(
        format="%(asctime)-15s [%(levelname)s] %(message)s",
        level=logging.INFO,
    )
    # Extract the filename from the path
    filename = os.path.basename(bi_vcfgz)
    # Extract the part before the first dot
    vcf_sample = filename.split('.')[0]
    out_dir = os.path.dirname(bi_vcfgz)

    logging.info(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CleanVCF: Start sample {vcf_sample}")

    os.makedirs(out_dir, exist_ok=True)

    if check_and_cleanup(bi_vcfgz, bihet_vcfgz, tmp_dir):
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        return

    filter_vcf(input_vcfgz, tmp_dir, vcf_sample, bihet_vcfgz, threads=1, cutoff=cutoff, pass_flag=pass_flag)
    
    compress_and_index_vcf(tmp_dir, vcf_sample, bi_vcfgz)

    if check_and_cleanup(bi_vcfgz, bihet_vcfgz, tmp_dir):
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
    
    logging.info(f"..... Done sample {vcf_sample}")

def main(): 
    if len(sys.argv) != 7:
        print("Usage: cleanVCF.py <tmp_dir> <in_vcf> <bi_vcf> <bihet_vcf> <cutoff> <pass_flag>")
        sys.exit(1)
    
    tmp_dir = sys.argv[1]
    input_vcfgz = sys.argv[2]
    bi_vcfgz = sys.argv[3]
    bihet_vcfgz = sys.argv[4]
    cutoff = sys.argv[5]
    skip_require_pass = False

    cleaning(tmp_dir, input_vcfgz, bi_vcfgz, bihet_vcfgz, cutoff, not skip_require_pass)

if __name__ == "__main__":
    main()
