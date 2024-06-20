#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

import os
import subprocess
import sys
from datetime import datetime
import tempfile

def mpileup(tmp_dir, save_dir, hetSNP, bam, output_file, ref, sample, chrom_start, chrom_end):
    now = datetime.now().strftime("%T")
    print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with {sample} : {now}")
    
    # Check and create tmp/out dir
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(save_dir, exist_ok=True)

    # mpileup for each chromosome
    for chr in range(chrom_start, chrom_end + 1):
        tmp_output = f"{tmp_dir}/{sample}_chr{chr}.pileup"
        run_mpileup_chrom(hetSNP, bam, tmp_output, chr, ref)
        
    # Combine all pileup files
    combined_tmp = f"{save_dir}/combined.gz"
    chrom_range = " ".join([f"{tmp_dir}/{sample}_chr{chr}.pileup" for chr in range(chrom_start, chrom_end + 1)])
    run_command(f"cat {chrom_range} | bgzip -@ {os.cpu_count()} -l3 -c > {combined_tmp}")
    os.rename(combined_tmp, output_file)
    
    now = datetime.now().strftime("%T")
    print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> done with {sample} : {now}")

def run_mpileup_chrom(hetSNP, bam, output, chr, ref):
    # Create a temporary file to hold the processed hetSNP content
    with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_hetSNP:
        # Process hetSNP content with tail, awk, and grep
        awk_grep_command = f"cat {hetSNP} | tail -n+2 | awk '{{print $1,$3}}' | grep 'chr{chr}\\s'"
        try:
            result = subprocess.run(awk_grep_command, shell=True, check=True, stdout=temp_hetSNP, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Error running command: {awk_grep_command}\n"
                               f"Command output: {e.stdout}\n"
                               f"Command stderr: {e.stderr.decode()}") from e

    try:
        # Run samtools mpileup command
        mpileup_command = f"samtools mpileup -d 0 -s -f {ref} -r chr{chr} -l {temp_hetSNP.name} {bam} > {output}"
        
        try:
            subprocess.run(mpileup_command, shell=True, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Error running command: {mpileup_command}\n"
                               f"Command output: {e.stdout}\n"
                               f"Command stderr: {e.stderr.decode()}") from e
    finally:
        # Clean up the temporary file
        os.remove(temp_hetSNP.name)

def run_command(command):
    subprocess.run(command, shell=True, check=True)

#python pileupREADS.py /data2/BEASTIE_example_output/NA12878_chr21/tmp/mpileup /data2/BEASTIE_example_output/NA12878_chr21 /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21_hetSNP.tsv /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/NA12878_chr21.bam /data2/BEASTIE_example_output/NA12878_chr21/NA12878.pileup.gz /data2/reference/hg19.fa NA12878_chr21 21 21

def main():
    if len(sys.argv) != 10:
        print("Usage: pileupREADS.py <tmp_dir> <save_dir> <hetSNP> <bam> <output_file> <ref> <sample_name> <chrom_start> <chrom_end>")
        sys.exit(1)
    
    tmp_dir = sys.argv[1]
    save_dir = sys.argv[2]
    hetSNP = sys.argv[3]
    bam = sys.argv[4]
    output_file = sys.argv[5]
    ref = sys.argv[6]
    sample = sys.argv[7]
    chrom_start = int(sys.argv[8])
    chrom_end = int(sys.argv[9])

    mpileup(tmp_dir, save_dir, hetSNP, bam, output_file, ref, sample, chrom_start, chrom_end)

if __name__ == "__main__":
    main()

