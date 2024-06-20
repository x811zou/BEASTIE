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

def run_command(command):
    subprocess.run(command, shell=True, check=True)

def convert_bamtosam(bam, sam, thread):
    # Run samtools mpileup command
    mpileup_command = f"samtools view -h --threads {thread} {bam} | bgzip -l3 -@ {thread} > {sam}.tmp"
    print(mpileup_command)
    subprocess.run(mpileup_command, shell=True, check=True, stderr=subprocess.PIPE)
    run_command(f"mv {sam}.tmp {sam}")
    run_command(f"tabix -p sam {sam}")

def run_simulation(sample, out_fwd, out_rev, tmp_dir, util_dir, genome, gff, bam, vcfgz):

    # Extract the base name (filename with extension)
    base_name = os.path.basename(bam)
    file_name_without_extension = os.path.splitext(base_name)[0]
    sam = os.path.join(tmp_dir, file_name_without_extension + ".sam")

    # Check if output files exist and are newer than the SAM file
    if (os.path.exists(out_fwd) and os.path.getsize(out_fwd) > 0 and 
        os.path.exists(out_rev) and os.path.getsize(out_rev) > 0 and 
        os.path.getmtime(out_fwd) > os.path.getmtime(sam) and 
        os.path.getmtime(out_rev) > os.path.getmtime(sam)):
        print(f"{out_fwd} exists and non-empty, newer than {sam}.")
        return

    # convert bam to sam
    convert_bamtosam(bam, sam, thread=1)

    # Create the working directory if it doesn't exist
    os.makedirs(tmp_dir, exist_ok=True)

    now = datetime.datetime.now().strftime("%T")
    print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> writing {sample} simulation seed file {seed_file}: {now}")
    seed = random.randint(0, 2**32 - 1)
    with open(seed_file, 'w') as f:
        f.write(f"{seed}\n")

    now = datetime.datetime.now().strftime("%T")
    print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> start with {sample} simulation : {now}")

    # Loop through 22 chromosomes
    job_ids = []
    for N in range(1, 23):
        job_name = f"{os.getenv('SLURM_JOB_NAME')}-chr{N}"
        job_id = subprocess.check_output([
            'sbatch', '--parsable',
            '--job-name', job_name,
            '--mem', mem,
            os.path.join(scripts_dir, '6__run_simulation_chrom.slurm'),
            sample, simulation_working_dir, util_dir, genome, gff,
            sam, vcfgz, str(depth), str(N), str(seed),
            str(simulate_allSNPs), str(haplotype)
        ]).strip().decode('utf-8')
        job_ids.append(job_id)

    # Wait for all jobs to finish
    wait_for_jobs(job_ids)

    now = datetime.datetime.now().strftime("%T")
    print(f">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> combining {sample} : {now}")

    fwd_tmp = os.path.join(simulation_working_dir, 'fwd.tmp.fastq.gz')
    rev_tmp = os.path.join(simulation_working_dir, 'rev.tmp.fastq.gz')

    with open(fwd_tmp, 'wb') as f_out:
        for N in range(1, 23):
            with open(os.path.join(simulation_working_dir, f'chr{N}_1.fastq.gz'), 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

    with open(rev_tmp, 'wb') as f_out:
        for N in range(1, 23):
            with open(os.path.join(simulation_working_dir, f'chr{N}_2.fastq.gz'), 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)

    shutil.move(fwd_tmp, out_fwd)
    shutil.move(rev_tmp, out_rev)

def run_chromosome_simulation(sample, out_path, util_dir, genome, gff, sam, vcfgz, depth, N, seed, snp, haplotype, simulator_py):
    snp_arg = "--all_snps" if snp else ""
    random_arg = "--random" if haplotype == "random" else ""

    print(f">>> chromosome: chr{N}")

    os.makedirs(out_path, exist_ok=True)

    out1 = os.path.join(out_path, f'chr{N}_1.fastq.gz')
    out2 = os.path.join(out_path, f'chr{N}_2.fastq.gz')

    cmd = [
        'python', simulator_py,
        util_dir, genome, gff, sam, vcfgz,
        '--out1', out1,
        '--out2', out2,
        '--read_depth', str(depth),
        '--chr', f'chr{N}',
        '--seed', str(seed),
        snp_arg,
        random_arg
    ]
    subprocess.check_call(cmd)

#python simulateREADS.py NA12878_chr21 /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21_fwd.fastq.gz /data2/BEASTIE_example_output/NA12878_chr21/NA12878_chr21_rev.fastq.gz /data2/BEASTIE_example_output/NA12878_chr21/tmp/simulation /data2/reference/two_bit /data2/reference/two_bit/hg19.2bit /data2/reference/gencode.v19.annotation.level12.gtf /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/NA12878_chr21.bam /home/scarlett/github/BEASTIE/BEASTIE_example/NA12878_chr21/HG001_GRCh37_GIAB.chr21.vcf.gz

if __name__ == "__main__":
    if len(sys.argv) != 10:
        print("Usage: simulateREADS.py <sample_name> <out_fwd> <out_rev> <tmp_dir> <util_dir> <genome> <gff> <bam> <vcfgz>")
        sys.exit(1)
    args = sys.argv[1:]
    run_simulation(*args)

