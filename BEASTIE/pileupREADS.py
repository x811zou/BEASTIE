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

def normalize_chromosome_name(chrom):
    """
    Normalize chromosome names by removing the 'chr' prefix (if present)
    and returning the name in uppercase for consistency.
    """
    return chrom.strip().lstrip("chr").upper()


def mpileup(tmp_dir, hetSNP, bam_gz, output_file, ref, chrom_start, chrom_end):
    # Extract the filename from the path
    filename = os.path.basename(hetSNP)
    # Extract the part before the first dot
    sample = filename.split('.')[0]

    logging.info(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> pileupReads: Start sample {sample}")

    # Check if output files exist and are newer than the SAM file
    if ((os.path.exists(output_file) and os.path.getsize(output_file) > 0) and 
        (os.path.getmtime(output_file) > os.path.getmtime(hetSNP))):
        logging.info(f"..... Check: {os.path.basename(output_file)} exists and non-empty, and is newer than {os.path.basename(hetSNP)} -- skipping pileup file generation")
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        return
    else:
        logging.info(f"..... Check: {os.path.basename(output_file)} does not exists or is older than {os.path.basename(hetSNP)} -- generating pileup file")

    # Check and create tmp/out dir
    os.makedirs(tmp_dir, exist_ok=True)

    # Decompress the bam.gz to a temporary BAM file
    temp_bam = tempfile.NamedTemporaryFile(delete=False, suffix=".bam")
    temp_bam.close()  # Close the file so samtools can write to it
    run_command(f"gzip -dc {bam_gz} | samtools view -b -o {temp_bam.name}")

    # Create an index for the temporary BAM file
    run_command(f"samtools index {temp_bam.name}")

    # # mpileup for each chromosome
    # for chr in range(chrom_start, chrom_end + 1):
    #     logging.info(f"..... Start with chr{chr}")
    #     tmp_output = f"{tmp_dir}/{sample}_chr{chr}.pileup"
    #     run_mpileup_chrom(hetSNP, temp_bam.name, tmp_output, chr, ref)
    # Define the list of chromosomes to process (special handling for 2A/2B)
    target_chromosomes = [str(i) for i in range(chrom_start, chrom_end + 1)]
    
    # Handle special case for chromosome 2A and 2B
    if "2A" in hetSNP or "2B" in hetSNP:
        target_chromosomes += ["2A", "2B"]

    # mpileup for each chromosome
    for chr in target_chromosomes:
        logging.info(f"..... Start with {chr}")
        tmp_output = f"{tmp_dir}/{sample}_chr{chr}.pileup"
        run_mpileup_chrom(hetSNP, temp_bam.name, tmp_output, chr, ref)

    # Combine only the pileup files that exist
    pileup_files = [f"{tmp_dir}/{sample}_chr{chr}.pileup" for chr in target_chromosomes if os.path.exists(f"{tmp_dir}/{sample}_chr{chr}.pileup")]

    if not pileup_files:
        raise RuntimeError("No pileup files were generated. Check input data and parameters.")

    combined_tmp = f"{tmp_dir}/combined.gz"
    run_command(f"cat {' '.join(pileup_files)} | bgzip -@ {os.cpu_count()} -l3 -c > {combined_tmp}")
    os.rename(combined_tmp, output_file)

    logging.info(f"..... Output {output_file}")
    logging.info(f"..... Done sample {sample}")

    # Clean up the temporary BAM file and its index
    os.remove(temp_bam.name)
    os.remove(temp_bam.name + ".bai")

def run_mpileup_chrom(hetSNP, bam, output, chr, ref):
    chr_variants = [f"{chr}"]

    # Read the first chromosome entry from hetSNP file to determine case
    with open(hetSNP, 'r') as f:
        # Skip the header and get the first data line
        next(f)  # Skip header
        first_line = next(f).strip().split()
        chr_in_hetSNP = first_line[0]  # Get the chromosome column (chr)

    # Normalize chr2 special cases dynamically based on hetSNP content
    if chr == "2":
        if "2A" in chr_in_hetSNP or "2B" in chr_in_hetSNP:
            chr_variants = ["2A", "2B"]
        elif "2a" in chr_in_hetSNP or "2b" in chr_in_hetSNP:
            chr_variants = ["2a", "2b"]
        else:
            logging.info("..... No chr2A/chr2B found in hetSNP, skipping special chromosomes.")
            return

    for chr_variant in chr_variants:
        logging.info(f"..... Processing chr{chr_variant}")
        tmp_output = output

        # Extract relevant lines from hetSNP for the current chromosome (exact case match)
        awk_grep_command = f"cat {hetSNP} | tail -n+2 | awk '{{print $1,$3}}' | grep 'chr{chr_variant}\\s'"
        try:
            with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_hetSNP:
                subprocess.run(awk_grep_command, shell=True, check=True, stdout=temp_hetSNP, stderr=subprocess.PIPE)

            # Define and run the samtools mpileup command using the detected chromosome case
            mpileup_command = f"samtools mpileup -d 0 -s -f {ref} -r chr{chr_variant} -l {temp_hetSNP.name} {bam} > {tmp_output}"
            subprocess.run(mpileup_command, shell=True, check=True, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Error running command: {awk_grep_command}\n"
                               f"Command output: {e.stdout}\n"
                               f"Command stderr: {e.stderr.decode()}") from e

        finally:
            # Ensure the temporary file is deleted
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

