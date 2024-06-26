#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================

import os
import subprocess
import datetime
import random
import shutil
import sys
import logging

current_script_path = os.path.abspath(__file__)
current_script_dir = os.path.dirname(current_script_path)
simulator_path = os.path.join(current_script_dir, '../spliced_simulator')
sys.path.append(simulator_path)
# print("Current script path:", current_script_path)
# print("Current script directory:", current_script_dir)
# print("Simulator path:", simulator_path)
# print("Sys.path:", sys.path)
import unbiased_spliced_rna

def run_command(command):
    subprocess.run(command, shell=True, check=True)

def convert_bamtosam(bamgz, samgz, thread):
    # Run samtools mpileup command
    #conversion_command = f"samtools view -h --threads {thread} {bam} | bgzip -l3 -@ {thread} > {sam}.tmp"
    conversion_command = f"gzip -cd {bamgz} | samtools view -h --threads {thread} - | bgzip -@ {thread} -c > {samgz}.tmp"
    subprocess.run(conversion_command, shell=True, check=True, stderr=subprocess.PIPE)
    run_command(f"mv {samgz}.tmp {samgz}")
    run_command(f"tabix -p sam {samgz}")
    logging.info("... start indexing sam.gz")


def run_simulation(chr_start, chr_end, outfastq_fwd, outfastq_rev, tmp_dir, twobit_dir, genome_2bit, gff, bamgz, vcfgz):
    # Extract the filename from the path
    filename = os.path.basename(bamgz)
    # Extract the part before the first dot
    sample = filename.split('.')[0]
    samgz = os.path.join(tmp_dir, sample + ".sam.gz")
    logging.info(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> simulateReads: start sample {sample}")
    # Check if output files exist and are newer than the SAM file
    if (os.path.exists(outfastq_fwd) and os.path.getsize(outfastq_fwd) > 0 and 
        os.path.exists(outfastq_rev) and os.path.getsize(outfastq_rev) > 0 and
        os.path.getmtime(outfastq_fwd) > os.path.getmtime(bamgz) and 
        os.path.getmtime(outfastq_rev) > os.path.getmtime(bamgz)):
        logging.info(f"..... Check: both {os.path.basename(outfastq_fwd)} & {os.path.basename(outfastq_rev)} exist and non-empty, newer than BAMGZ -- skipping simulation")
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        return
    else:
        logging.info(f"..... Check: Not both {os.path.basename(outfastq_fwd)} exist or not newer than BAMGZ -- starting simulation")

    # Create the working directory if it doesn't exist
    os.makedirs(tmp_dir, exist_ok=True)

    # convert bam to sam
    logging.info(f"..... Start converting {os.path.basename(bamgz)}")
    convert_bamtosam(bamgz, samgz, thread=1)
    logging.info(f"..... Finish converting to {os.path.basename(samgz)}")

    # Run simulation for each chromosome
    seed = random.randint(0, 2**32 - 1)
    seed_file = os.path.join(tmp_dir, 'seed.txt')
    with open(seed_file, 'w') as f:
        f.write(f"{seed}\n")

    logging.info(f"..... Start with even fastq reads simulation")

    # Loop through 22 chromosomes
    for N in range(int(chr_start), int(chr_end)+1):
        outFile1name = os.path.join(tmp_dir, f'chr{N}_1.fastq.gz')
        outFile2name = os.path.join(tmp_dir, f'chr{N}_2.fastq.gz')

        logging.info(f"..... chr{N} fastq reads simulation")
        chrN = f"chr{N}"
        unbiased_spliced_rna.simulation_pipeline(gff = gff, target_chromosome = chrN, target_gene = None, twoBitDir = twobit_dir, genome2bit = genome_2bit, outFile1 = outFile1name, outFile2 = outFile2name, SNPs = True, vcf = vcfgz, samgz = samgz, random = False, read_depth = 100, max_genes = 0)

    logging.info(f"..... Finish chr{int(chr_start)} to chr{int(chr_end)} simulation")
    logging.info(f"..... start combining {sample}")

    fwd_tmp = os.path.join(tmp_dir, 'fwd.tmp.fastq.gz')
    rev_tmp = os.path.join(tmp_dir, 'rev.tmp.fastq.gz')

    with open(fwd_tmp, 'wb') as f_out:
        for N in range(int(chr_start), int(chr_end)+1):
            with open(os.path.join(tmp_dir, f'chr{N}_1.fastq.gz'), 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

    with open(rev_tmp, 'wb') as f_out:
        for N in range(int(chr_start), int(chr_end)+1):
            with open(os.path.join(tmp_dir, f'chr{N}_2.fastq.gz'), 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

    shutil.move(fwd_tmp, outfastq_fwd)
    shutil.move(rev_tmp, outfastq_rev)
    logging.info(f"..... Output {os.path.basename(outfastq_fwd)} {os.path.basename(outfastq_rev)}")
    
#python simulateREADS.py 21 21 /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21_fwd.fastq.gz /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21_rev.fastq.gz /data2/BEASTIE_example_output/NA12878_chr21/tmp/simulation /data2/reference/two_bit_linux /data2/reference/hg19/hg19.2bit /data2/reference/gencode.v19.annotation.level12.gtf /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/NA12878_chr21.bam.gz /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/NA12878_chr21.bi.vcf.gz

if __name__ == "__main__":
    if len(sys.argv) != 11:
        print("Usage: simulateREADS.py <chr_start> <chr_end> <out_fwd> <out_rev> <tmp_dir> <util_dir> <genome> <gff> <bamgz> <vcfgz>")
        sys.exit(1)
    args = sys.argv[1:]
    logging.basicConfig(
        format="%(asctime)-15s [%(levelname)s] %(message)s",
        level=logging.INFO,
    ) 
    run_simulation(*args)

