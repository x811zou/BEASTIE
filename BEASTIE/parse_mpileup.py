#!/usr/bin/env python
# =========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
# =========================================================================
import gzip
import logging
import os

from cyvcf2 import VCF


def if_ref(read, ref, allele):
    if read.lower() == ref or read.upper() == ref:
        return True
    else:
        return False


def if_alt(read, ref, allele):
    if read.lower() in allele.get(ref) or read.upper() in allele.get(ref):
        return True
    else:
        return False


def isHeterozygous(genotype):
    Homo = ["0|0", "1|1", "2|2", "3|3", "4|4", "5|5", "6|6", "7|7", "0/0", "1/1", "2/2"]
    if genotype in Homo:
        return False
    else:
        return True


def count_raw_depth(reads):
    qualified_reads = []
    i = 0
    for x in reads:
        if x != ">" and x != "<":
            i = i + 1
            qualified_reads.append(x)
    return i


def ref_alt_dict(REF, ALT):
    allele = {}
    allele[REF] = allele.get(REF, set())
    for j in ALT:
        allele[REF].add(str(j))
    return allele


def GATK_ParseMpileup(
    record, hets, indels, snp, counter, bi, i, pipeup_dict, min_cov, either_cov
):
    columns = (
        str(record).strip("\n").split("\t")
    )  # parse it because cyvcf2 does not support extracting the genotype]
    ref = ""
    alt = ""
    rsid = "."  # read in the dict of pileup result because it is relatively small file
    vcf_key = str(record.CHROM) + "-" + str(record.start + 1)
    pipeup_line = pipeup_dict.get(vcf_key).strip("\n").split("\t")
    alleles = ref_alt_dict(record.REF, record.ALT)
    ref = ", ".join([key for key in alleles.keys()])
    alt = ", ".join([", ".join(value) for value in alleles.values()])
    length_key = len(alleles[str(ref)])
    rsid = columns[2]
    if ref != "" and alt != "":  # print(list(pipeup_line[4]))
        reads = list(pipeup_line[4])  # read bases
        reads = [
            x
            for x in reads
            if (x != "$" and x != "^" and x != "~" and x != '"' and x != "!")
        ]  # remove this non-sense symbols
        block = 0
        out_reads = []
        for read in reads:
            if block == -1:
                i = i + 1
                if ord(read) >= 48 and ord(read) <= 57:  # up to 9 insertion/deletion
                    block_str += read
                else:
                    if block_str != "":
                        block = int(block_str) - 1
            elif read == "+" or read == "-":
                block = -1
                block_str = ""
            elif block > 0:
                i = i + 1
                block -= 1
            elif block == 0:
                out_reads.append(read)
        reads = out_reads
        baseqs = list(pipeup_line[5])  # base quality
        mapqs = list(pipeup_line[6].strip())  # mapping quality
        low_baseq = 0
        low_mapq = 0
        ref_count = 0
        alt_count = 0
        other_count = 0
        raw_depth = count_raw_depth(
            reads
        )  # number of qualified reads    #print("read: "+str(out_reads))
        for read, baseq, mapq in zip(reads, baseqs, mapqs):
            if read != "<" and read != ">":  # D-68; E-69,F-70; H-72; J-74
                basequal = ord(baseq) - 33
                mapqual = ord(mapq) - 33
                if basequal >= 0 and mapqual >= 0:  # min_mapq=93    # count this base
                    # logging.info("..... reads: {0}".format(read))
                    if (
                        read == "." or read == "," or if_ref(read, record.REF, alleles)
                    ):  # print("ref: "+str(read))
                        ref_count += 1
                        # logging.info("..... ref - {0} - {1}".format(read, ref_count))
                    elif if_alt(read, record.REF, alleles):  # print("alt: "+str(read))
                        alt_count += 1
                        # logging.info("..... alt - {0} - {1}".format(read, alt_count))
                    else:
                        other_count += 1  # print("other: "+str(read))
                if basequal < 0:
                    low_baseq += 1
                if mapqual < 0:
                    low_mapq += 1
        totalCount = ref_count + alt_count
        if totalCount > 0:
            ratio = alt_count / totalCount
        else:
            ratio = 0
        if_biallelic = "NA"
        if_Indel = "NA"
        if_sv = "NA"
        if_snp = "NA"
        if (
            totalCount >= int(min_cov)
            and ref_count >= int(either_cov)
            and alt_count >= int(either_cov)
        ):
            hets += 1
            if length_key == 1:
                if_biallelic = "Y"
                bi += 1
            else:
                if_biallelic = "N"
            if record.is_indel:
                if_Indel = "Y"
                indels += 1
            else:
                if_Indel = "N"
            if record.is_sv:
                if_sv = "Y"
            else:
                if_sv = "N"
            if record.is_snp:
                if_snp = "Y"
                snp += 1
            else:
                if_snp = "N"
        return (
            columns[0],
            columns[1],
            rsid,
            str(ref),
            str(ref_count),
            str(alt),
            str(alt_count),
            str(totalCount),
            str(ratio),
            if_Indel,
            if_sv,
            if_snp,
            if_biallelic,
            str(low_mapq),
            str(low_baseq),
            str(raw_depth),
            str(other_count),
            hets,
            counter,
        )
    raise ValueError("Empty content")


def Parse_mpileup_allChr(
    vcfFile_sample, inputGZ, pileup_file, min_cov, either_cov, output
):
    pipeup_dict = {}  # DEBUGGING
    logging.info("..... Finish loading pileup file: {0}".format(pileup_file))
    logging.info("..... Finish loading vcf file: {0}".format(inputGZ))
    out_stream = open(output, "w")
    logging.info("..... Start parsing mpile up")
    out_stream.write(
        "contig\tposition\tvariantID\trefAllele\trefCount\taltAllele\taltCount\ttotalCount\taltRatio\tif_Indel\tif_SV\tif_SNP\tif_biallelic\tlowMAPQDepth\tlowBaseQDepth\trawDepth\totherCount\n"
    )

    with gzip.open(pileup_file, "rt") if os.path.splitext(pileup_file)[
        1
    ] == ".gz" else open(pileup_file, "r") as stream_in:
        for i, line in enumerate(
            stream_in
        ):  # start reading in my pileup results line by line
            chr_info = line.split("\t")[0].strip("chr")
            start_info = line.split("\t")[1]
            key_info = chr_info + "-" + start_info
            pipeup_dict[key_info] = line  # save key: [1] value: line
        counter = 0
        hets = 0
        snp = 0
        indels = 0
        bi = 0
        vcf = VCF(inputGZ)
        vcf.set_samples(
            [vcfFile_sample]
        )  # prefix must be consistent with the sample name in VCF
        for record in vcf:
            counter += 1
            if counter % 1000000 == 0:
                logging.info("..... processed {0} lines".format(counter))
            columns = (
                str(record).strip("\n").split("\t")
            )  # parse it because cyvcf2 does not support extracting the genotype]
            vcf_key = str(record.CHROM) + "-" + str(record.start + 1)
            pipeup_line = pipeup_dict.get(vcf_key)
            if not pipeup_line:
                continue
            if not isHeterozygous(columns[-1]):
                continue
            (
                chrNum,
                start,
                rsid,
                ref,
                ref_count,
                alt,
                alt_count,
                total_count,
                ratio,
                if_Indel,
                if_SV,
                if_SNP,
                if_biallelic,
                low_mapg,
                low_baseq,
                raw_depth,
                other_count,
                hets,
                counter,
            ) = GATK_ParseMpileup(
                record,
                hets,
                indels,
                snp,
                counter,
                bi,
                i,
                pipeup_dict,
                min_cov,
                either_cov,
            )
            if chrNum is None:
                continue
            out_stream.write(
                "\t".join(
                    [
                        chrNum,
                        start,
                        rsid,
                        ref,
                        ref_count,
                        alt,
                        alt_count,
                        total_count,
                        ratio,
                        if_Indel,
                        if_SV,
                        if_SNP,
                        if_biallelic,
                        low_mapg,
                        low_baseq,
                        raw_depth,
                        other_count,
                        "\n",
                    ]
                )
            )
    out_stream.close()
    stream_in.close()
    logging.info("Pileup dimension: {0}".format(len(pipeup_dict)))
    logging.info("VCF dimension: {0}".format(counter))
    logging.info(
        "{0} het sites found with total counts above {1} and either allele above {2}".format(
            hets, min_cov, either_cov
        )
    )
    logging.info("Parsed pileup reads file saved to {0}".format(output))
