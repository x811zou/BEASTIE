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
import gzip

def setup_logging(output_file):
    """Set up logging with a dynamic log file name inside a 'log' directory."""
    output_dir = os.path.dirname(os.path.abspath(output_file))
    log_dir = os.path.join(output_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    current_date = datetime.now().strftime("%b-%d-%Y")
    log_filename = os.path.join(log_dir, f"BEASTIE-pileupReads-{current_date}.log")
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(
        filename=log_filename,
        filemode='w',
        format="%(asctime)-15s [%(levelname)s] %(message)s",
        level=logging.INFO,
    )
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(asctime)-15s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(console_handler)
    logging.info(f"Logging initialized. Log file: {log_filename}")

def normalize_chromosome_name(chrom):
    """
    Normalize chromosome names by removing the 'chr' prefix (if present)
    and returning the name in uppercase for consistency.
    """
    return chrom.strip().lstrip("chr").upper()


def mpileup(tmp_dir, hetSNP, bam_gz, output_file, ref, chrom_start, chrom_end):
    setup_logging(output_file)  # Set up logging at the beginning
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
    # Define the list of chromosomes to process (special handling for 2A/2B)
    target_chromosomes = [str(i) for i in range(chrom_start, chrom_end + 1)]
    
    summary = []

    # mpileup for each chromosome
    for chrN in target_chromosomes:
        logging.info(f"..... Start with chr{chrN}")
        tmp_output = f"{tmp_dir}/{sample}_chr{chrN}.pileup"
        chr_variants = detect_chr_variants(hetSNP) if chrN == "2" else [chrN]
        for chr_variant in chr_variants:
            chr_str, num_genes_input, num_hets_input, num_hets_output = run_mpileup_chrom(hetSNP, temp_bam.name, tmp_output, chr_variant, ref)
            summary.append((chr_str, num_genes_input, num_hets_input, num_hets_output))

    # Combine only the pileup files that exist
    pileup_files = [f"{tmp_dir}/{sample}_chr{chrN}.pileup" for chrN in target_chromosomes if os.path.exists(f"{tmp_dir}/{sample}_chr{chrN}.pileup")]

    if not pileup_files:
        logging.error("No pileup files were generated. Check input data and parameters.")
        return

    combined_tmp = f"{tmp_dir}/combined.gz"
    run_command(f"cat {' '.join(pileup_files)} | bgzip -@ {os.cpu_count()} -l3 -c > {combined_tmp}")
    os.rename(combined_tmp, output_file)

    logging.info(f".......... Output {output_file}")
    logging.info(f"..... Done sample {sample}")

    # Generate and log the summary table
    summarize_pileup(summary)

    # Clean up the temporary BAM file and its index
    os.remove(temp_bam.name)
    os.remove(temp_bam.name + ".bai")

def detect_chr_variants(hetSNP_file):
    """
    Detect whether chr2, chr2A, or chr2B (case-sensitive) exists in the hetSNP file.
    Returns a list of chromosome variants to process.
    """
    chr_variants = set()
    with open(hetSNP_file, 'r') as infile:
        for line_num, line in enumerate(infile):
            if line_num == 0:  # Skip header
                continue
            parts = line.strip().split("\t")
            if len(parts) < 1:
                continue
            chrom = parts[1]  # Keep the chromosome name as-is
            if chrom in {"2", "2A", "2B", "2a", "2b"}:
                chr_variants.add(chrom)

    # If chr2A/chr2B are present, ignore chr2
    if "2" in chr_variants and any(chrom in chr_variants for chrom in {"2A", "2B", "2a", " 2b"}):
        chr_variants.discard("2")
    if not chr_variants:              
        raise ValueError("No valid chr2 variants found in the hetSNP file.")
    return sorted(chr_variants)  # Return sorted list for consistent processing order


def run_mpileup_chrom(hetSNP, bam, output, chr_variant, ref):
    # Read the first chromosome entry from hetSNP file to determine case

    # Detect chr2 variants dynamically if chr == "chr2"

    chr_str_variant = f"chr{chr_variant}"
    print(chr_str_variant)
    positions = [(line.split()[0], line.split()[2]) for line in open(hetSNP) if line.split()[0] == chr_str_variant]

    if not positions:
        raise ValueError(f"No positions found for {chr_str_variant} in {hetSNP}.")

    logging.info(f".......... Processing chr{chr_variant}")

    # Extract relevant lines from hetSNP for the current chromosome (exact case match)
    awk_grep_command = f"cat {hetSNP} | tail -n+2 | awk '{{print $1,$3}}' | grep 'chr{chr_variant}\\s'"
    try:
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_hetSNP:
            subprocess.run(awk_grep_command, shell=True, check=True, stdout=temp_hetSNP, stderr=subprocess.PIPE)

        tmp_output = output
        # Define and run the samtools mpileup command using the detected chromosome case
        mpileup_command = f"samtools mpileup -d 0 -s -f {ref} -r chr{chr_variant} -l {temp_hetSNP.name} {bam} > {tmp_output}"
        logging.info(f".......... Running mpileup command: {mpileup_command}")
        try:
            subprocess.run(mpileup_command, shell=True, check=True, stderr=subprocess.PIPE)
            logging.info(f".......... Completed mpileup for chr{chr_variant}.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error during mpileup for chr{chr_variant}. stderr: {e.stderr.decode()}")
            raise

        num_genes_input, num_snps_input, num_snps_output = get_genes_and_snps_from_pileup(hetSNP, output, target_chrom=chr_str_variant)

        return chr_str_variant, num_genes_input, num_snps_input, num_snps_output
    
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {awk_grep_command}")
        logging.error(f"stderr: {e.stderr.decode()}")
        raise RuntimeError(f"Error running command: {awk_grep_command}") from e


    finally:
        if os.path.exists(temp_hetSNP.name):
            os.remove(temp_hetSNP.name)


def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {command}")
        logging.error(f"Error message: {e.stderr.decode() if e.stderr else 'No stderr available'}")
        raise


def get_genes_and_snps_from_pileup(hetSNP_file, pileup_file, target_chrom):
    unique_genes = set()
    num_snps_input = 0

    # Parse hetSNP_file
    with open(hetSNP_file, 'r') as infile:
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) >= 4 and parts[0] == target_chrom:
                unique_genes.add(parts[3])
                num_snps_input += 1

    num_genes_input = len(unique_genes)

    # Parse pileup_file
    def is_gzipped(filepath):
        with open(filepath, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'

    open_func = gzip.open if is_gzipped(pileup_file) else open
    with open_func(pileup_file, 'rt') as infile:
        num_snps_pileup = sum(1 for _ in infile)

    return num_genes_input, num_snps_input, num_snps_pileup


def summarize_pileup(summary):
    """Generate a final summary table for the pileup process."""

    logging.info("\nFinal Summary of Pileup by Chromosome:")
    logging.info(f"{'Chromosome':<12}{'Genes (input)':<15}{'Het SNPs (input)':<20}{'Het SNPs (output)':<20}")
    logging.info("-" * 70)

    total_genes_input = total_snps_input = total_snps_output = 0

    for chr, genes_input, snps_input, snps_output in summary:
        logging.info(f"{chr:<12}{genes_input:<15}{snps_input:<20}{snps_output:<20}")
        total_genes_input += genes_input
        total_snps_input += snps_input
        total_snps_output += snps_output

    # Totals row
    logging.info("-" * 70)
    logging.info(f"{'Total':<12}{total_genes_input:<15}{total_snps_input:<20}{total_snps_output:<20}")

    if total_snps_input > total_snps_output:
        logging.info(f"{total_snps_input - total_snps_output} heterozygous SNPs were filtered out.")
        logging.info("Filtered out SNPs during samtoos mpileup process are likely due to: low coverage (max read depth -d 8000 = Filters out regions with excessive coverage that might indicate PCR artifacts, repeat regions, or other sequencing biases), low quality (base quality score -q 20 = 1% probability of error in base calling; mapping quality -Q score 0 = a read that could not be uniquely aligned), or other filtering criteria.")
    else:
        logging.info("All SNPs were successfully output in the pileup files.")

#python pileupREADS.py /data2/BEASTIE_example_output/NA12878_chr21/tmp/mpileup /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21_hetSNP.tsv /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/NA12878_chr21.bam.gz /data2/reference/hg19/hg19.fa 21 21 /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21.pileup.gz

def validate_files(files):
    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input file {f} not found.")


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

    try:
        validate_files([hetSNP, bamgz, ref])
        mpileup(tmp_dir, hetSNP, bamgz, output_file, ref, chrom_start, chrom_end)
    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        raise



if __name__ == "__main__":
    main()