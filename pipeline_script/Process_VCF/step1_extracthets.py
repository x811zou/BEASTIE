import os.path
import subprocess
import sys

# We have to adjust the import paths so python can find all the file_tools modules
sys.path.append(os.path.abspath('..'))
sys.path.append(os.path.abspath('../file_tools'))

from file_tools.GffTranscriptReader import GffTranscriptReader

def is_homozygous(genotype):
    homozygous_genotypes = {"0|0","1|1","2|2","3|3","4|4","5|5","6|6","7|7","0/0","1/1","2/2"}
    return genotype in homozygous_genotypes

def shell_out(cmd):
    output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    return output.decode(sys.stdout.encoding)

def count_all_het_sites(sample, vcf_filename, annotation_dir, output_filename):
    print(f">>>>>>>>>>>>>>>> sample {sample}")

    out_file = open(output_filename, "w")
    out_file.write("chr\tchrN\tgeneID\tgenomicCoord_pos\ttranscriptCoord\tSNP_id\tgenotype\n")

    print("- all transcripts:")
    for chromosome in range(1,23):
        reader = GffTranscriptReader()
        print(f"we are working with chr{chromosome}")
        gene_list = reader.loadGenes(os.path.join(annotation_dir, f"chr{chromosome}.gtf"))
        print(f"{len(gene_list)} Genes loaded")

        for gene in gene_list:
            transcript = gene.longestTranscript()
            geneID = transcript.getGeneId()      # column 5
            chrom = transcript.getSubstrate()      # column 1
            chromN = chrom.strip("chr")
            begin = gene.getBegin()    # column 7
            end = gene.getEnd()      # column 8

            output = shell_out(f"tabix {vcf_filename} {chromN}:{begin}-{end}")

            if(len(output) <= 9):
                continue

            lines = output.split("\n")
            for line in lines:
                fields = line.split("\t")
                if(fields[6] != "PASS"):
                    continue
                pos = fields[1]
                rs = fields[2]
                genotype = fields[9].split(':')[0]

                # Ignore homozygous genotypes
                if is_homozygous(genotype):
                    continue

                transcriptCoord = transcript.mapToTranscript(int(pos))

                out_file.write(f"{chrom}\t{chromN}\t{geneID}\t{pos}\t{transcriptCoord}\t{rs}\t{genotype}\n")

        print(f"Finished with chr{chromosome}")

    out_file.close()

if __name__ == "__main__":
    sample = sys.argv[1]
    vcf_filename = sys.argv[2]
    annotation_dir = sys.argv[3]
    out_dir = sys.argv[4]
    # There is only one VCF file for GSD case
    hets_meta_output_dir = os.path.join(out_dir, "hetsMeta/")
    output_filename = os.path.join(hets_meta_output_dir, "Allchr_hets_all_transcript.tsv")
    count_all_het_sites(sample, vcf_filename, annotation_dir, output_filename)
