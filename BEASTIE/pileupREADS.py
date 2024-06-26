#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import sys
import os
import subprocess
import tempfile
from datetime import datetime
import logging
import shutil

def mpileup(tmp_dir, hetSNP, bam_gz, output_file, ref, chrom_start, chrom_end):
    base_name = os.path.basename(output_file)
    file_name_without_extension = os.path.splitext(base_name)[0]
    sample = os.path.splitext(file_name_without_extension)[0]

    # Check if output files exist and are newer than the SAM file
    if ((os.path.exists(output_file) and os.path.getsize(output_file) > 0) and 
        (os.path.getmtime(output_file) > os.path.getmtime(hetSNP))):
        print(f"{output_file} exists and non-empty, and is newer than {hetSNP}")
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        return
    
    now = datetime.now().strftime("%T")
    logging.info(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start sample {sample} : {now}")

    # Check and create tmp/out dir
    os.makedirs(tmp_dir, exist_ok=True)

    # Decompress the bam.gz to a temporary BAM file
    temp_bam = tempfile.NamedTemporaryFile(delete=False, suffix=".bam")
    temp_bam.close()  # Close the file so samtools can write to it
    run_command(f"gzip -dc {bam_gz} | samtools view -b -o {temp_bam.name}")

    # Create an index for the temporary BAM file
    run_command(f"samtools index {temp_bam.name}")

    # mpileup for each chromosome
    for chr in range(chrom_start, chrom_end + 1):
        logging.info(f"... start with chr{chr}")
        tmp_output = f"{tmp_dir}/{sample}_chr{chr}.pileup"
        run_mpileup_chrom(hetSNP, temp_bam.name, tmp_output, chr, ref)
        
    # Combine all pileup files
    combined_tmp = f"{tmp_dir}/combined.gz"
    chrom_range = " ".join([f"{tmp_dir}/{sample}_chr{chr}.pileup" for chr in range(chrom_start, chrom_end + 1)])
    run_command(f"cat {chrom_range} | bgzip -@ {os.cpu_count()} -l3 -c > {combined_tmp}")
    os.rename(combined_tmp, output_file)
    
    now = datetime.now().strftime("%T")
    logging.info(f"... output {output_file}")
    logging.info(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> done sample {sample} : {now}")

    # Clean up the temporary BAM file and its index
    os.remove(temp_bam.name)
    os.remove(temp_bam.name + ".bai")

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



#python pileupREADS.py /data2/BEASTIE_example_output/NA12878_chr21/tmp/mpileup /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21_hetSNP.tsv /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/NA12878_chr21.bam.gz /data2/reference/hg19/hg19.fa 21 21 /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21.pileup.gz

def main():
    if len(sys.argv) != 8:
        print("Usage: pileupREADS.py <tmp_dir> <in_hetSNP> <in_bamgz> <in_ref> <chrom_start> <chrom_end> <out_file>")
        sys.exit(1)
    
    tmp_dir = sys.argv[1]
    hetSNP = sys.argv[2]
    bamgz = sys.argv[3]
    ref = sys.argv[4]
    chrom_start = int(sys.argv[5])
    chrom_end = int(sys.argv[6])
    output_file = sys.argv[7]

    logging.basicConfig(
        format="%(asctime)-15s [%(levelname)s] %(message)s",
        level=logging.INFO,
    )

    mpileup(tmp_dir, hetSNP, bamgz, output_file, ref, chrom_start, chrom_end)

if __name__ == "__main__":
    main()

