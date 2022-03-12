#!/usr/bin/env python
# =========================================================================
# 2021 Xue Zou (xue.zou@duke.edu)
# =========================================================================
from collections import namedtuple
import csv
import gzip
import logging
import multiprocessing
import os
import sqlite3
import pandas as pd
from pkg_resources import resource_filename
import math
import pandas
import datetime
import requests
from io import StringIO

from .helpers import flatten, runhelper

ANNOTATION_ALGORITH = "parallel_join"  #'constant_memory'
MAX_WORKERS = 4


def annotateAF(ancestry, hetSNP, out_AF):
    if ANNOTATION_ALGORITH == "parallel_join":
        annotateAFConstantMemoryParallel(ancestry, hetSNP, out_AF)
    elif ANNOTATION_ALGORITH == "constant_memory":
        AF_file = resource_filename("BEASTIE", "reference/AF/AF_1_22_trimmed2.csv.gz")
        annotateAFConstantMemory(ancestry, hetSNP, out_AF, AF_file)
    else:
        AF_file = resource_filename("BEASTIE", "reference/AF/AF_1_22_trimmed2.csv")
        annotateAFPandas(ancestry, hetSNP, out_AF, AF_file)

    logging.info("..... finish annotating AF for SNPs, file save at {0}".format(out_AF))


def annotateAFPandas(ancestry, hetSNP, out_AF, AF_file):
    logging.info("..... start reading 1000 Genome AF annotation file")
    AF = pd.read_csv(AF_file, header=0, sep=",", engine="c", na_filter=False)
    logging.info("..... finish reading 1000 Genome AF annotation file")
    data = pd.read_csv(hetSNP, sep="\t", header=0, index_col=False)
    if ancestry == "EUR":
        AF = AF[["chr", "pos", "rsid", "EUR_AF"]]
        AF = AF.rename(columns={"EUR_AF": "AF"})
    elif ancestry == "AFR":
        AF = AF[["chr", "pos", "rsid", "AFR_AF"]]
        AF = AF.rename(columns={"AFR_AF": "AF"})
    elif ancestry == "EAS":
        AF = AF[["chr", "pos", "rsid", "EAS_AF"]]
        AF = AF.rename(columns={"EAS_AF": "AF"})
    elif ancestry == "AMR":
        AF = AF[["chr", "pos", "rsid", "AMR_AF"]]
        AF = AF.rename(columns={"AMR_AF": "AF"})
    elif ancestry == "SAS":
        AF = AF[["chr", "pos", "rsid", "SAS_AF"]]
        AF = AF.rename(columns={"SAS_AF": "AF"})
    data_AF = pd.merge(data, AF, on=["chr", "pos"], how="left")
    data_AF = data_AF.drop_duplicates()
    data_AF.to_csv(out_AF, sep="\t")


def annotateCHRLines(hetSNP_lines, chr, ancestry, hetSNP_chr_index, hetSNP_pos_index):
    AF_filename = resource_filename("BEASTIE", f"reference/AF/AF_{chr}.csv.gz")
    with gzip.open(AF_filename, mode="rt", newline="") as AF_file:
        af_reader = csv.reader(AF_file, delimiter=",", dialect="unix")
        af_header = next(af_reader)

        af_chr_index = af_header.index("chr")
        af_pos_index = af_header.index("pos")
        af_af_col_index = af_header.index(f"{ancestry}_AF")
        af_rsid_col_index = af_header.index("rsid")

        af_rows_left = True

        af_row = next(af_reader)
        af_chr, af_pos = af_row[af_chr_index], int(af_row[af_pos_index])
        af_chr_n = int(af_chr.strip("chr"))

        output = []
        for row in hetSNP_lines:
            chr, pos = row[hetSNP_chr_index], int(row[hetSNP_pos_index])
            chr_n = int(chr.strip("chr"))

            while af_rows_left and (
                (chr_n > af_chr_n) or (chr_n == af_chr_n and pos > af_pos)
            ):
                try:
                    af_row = next(af_reader)
                    af_chr, af_pos = af_row[af_chr_index], int(af_row[af_pos_index])
                    af_chr_n = int(af_chr.strip("chr"))
                except StopIteration:
                    af_row = None
                    af_rows_left = False

            if af_row is not None and chr_n == af_chr_n and pos == af_pos:
                output.append(
                    row + [af_row[af_rsid_col_index], af_row[af_af_col_index]]
                )
            else:
                output.append(row + ["", ""])
    return output


def annotateAFConstantMemoryParallel(ancestry, hetSNP, out_AF):
    logging.info("..... start annotating AF with parallel constant memory join")

    with multiprocessing.Pool(MAX_WORKERS) as pool, open(
        hetSNP, newline=""
    ) as hetSNPfile, open(out_AF, "w", newline="") as outFile:
        hetSNP_reader = csv.reader(hetSNPfile, delimiter="\t", dialect="unix")
        hetSNP_header = next(hetSNP_reader)

        hetSNP_chr_index = hetSNP_header.index("chr")
        hetSNP_pos_index = hetSNP_header.index("pos")

        output_handles = []

        current_chr = None
        current_rows = None

        for row in hetSNP_reader:
            chr = row[hetSNP_chr_index]
            if chr != current_chr:
                if current_chr is not None:
                    output_handles.append(
                        pool.apply_async(
                            annotateCHRLines,
                            (
                                current_rows,
                                current_chr,
                                ancestry,
                                hetSNP_chr_index,
                                hetSNP_pos_index,
                            ),
                        )
                    )
                current_chr = chr
                current_rows = []
            current_rows.append(row)
        if current_rows is not None and len(current_rows) > 0:
            output_handles.append(
                pool.apply_async(
                    annotateCHRLines,
                    (
                        current_rows,
                        current_chr,
                        ancestry,
                        hetSNP_chr_index,
                        hetSNP_pos_index,
                    ),
                )
            )

        pool.close()

        out_writer = csv.writer(
            outFile, delimiter="\t", dialect="unix", quoting=csv.QUOTE_MINIMAL
        )
        out_writer.writerow(hetSNP_header + ["rsid", "AF"])

        for handle in output_handles:
            handle.wait()
            output = handle.get()
            out_writer.writerows(output)


def annotateAFConstantMemory(ancestry, hetSNP, out_AF, AF_file):
    logging.info("..... start annotating AF with constant memory join")

    with gzip.open(AF_file, mode="rt", newline="") as affile, open(
        hetSNP, newline=""
    ) as hetSNPfile, open(out_AF, "w", newline="") as outFile:
        af_reader = csv.reader(affile, delimiter=",", dialect="unix")
        hetSNP_reader = csv.reader(hetSNPfile, delimiter="\t", dialect="unix")

        af_header = next(af_reader)
        hetSNP_header = next(hetSNP_reader)

        hetSNP_chr_index = hetSNP_header.index("chr")
        hetSNP_pos_index = hetSNP_header.index("pos")

        af_chr_index = af_header.index("chr")
        af_pos_index = af_header.index("pos")
        af_af_col_index = af_header.index(f"{ancestry}_AF")
        af_rsid_col_index = af_header.index("rsid")

        af_rows_left = True

        af_row = next(af_reader)
        af_chr, af_pos = af_row[af_chr_index], int(af_row[af_pos_index])
        af_chr_n = int(af_chr.strip("chr"))

        rows_written = 0
        rows_read = 0

        out_writer = csv.writer(
            outFile, delimiter="\t", dialect="unix", quoting=csv.QUOTE_MINIMAL
        )
        out_writer.writerow(hetSNP_header + ["rsid", "AF"])

        for row in hetSNP_reader:
            chr, pos = row[hetSNP_chr_index], int(row[hetSNP_pos_index])
            chr_n = int(chr.strip("chr"))

            while af_rows_left and (
                (chr_n > af_chr_n) or (chr_n == af_chr_n and pos > af_pos)
            ):
                try:
                    af_row = next(af_reader)
                    af_chr, af_pos = af_row[af_chr_index], int(af_row[af_pos_index])
                    af_chr_n = int(af_chr.strip("chr"))
                    rows_read += 1
                except StopIteration:
                    af_row = None
                    af_rows_left = False

            if af_row is not None and chr_n == af_chr_n and pos == af_pos:
                out_writer.writerow(
                    row + [af_row[af_rsid_col_index], af_row[af_af_col_index]]
                )
            else:
                out_writer.writerow(row + ["", ""])
            rows_written += 1
    logging.info("..... end annotating AF")


def annotateLD(
    prefix,
    ancestry,
    file_for_LDannotation,
    out_path,
    LD_token,
    chr_start,
    chr_end,
    meta,
):
    annotateLD_cache(file_for_LDannotation, meta, ancestry, LD_token)

    annotate_ld_new = resource_filename("BEASTIE", "annotate_LD_new.R")
    beastie_wd = resource_filename("BEASTIE", ".")

    if not os.path.isfile(meta):
        cmd = f"Rscript --vanilla {annotate_ld_new} {prefix} {ancestry} {file_for_LDannotation} {out_path} {LD_token} {chr_start} {chr_end} {meta} {beastie_wd}"
        runhelper(cmd)
        logging.info(f"..... finish annotating LD for SNP pairs, file save at {meta}")
    else:
        logging.info(
            f"..... skip annotating LD for SNP pairs, file already saved at {meta}"
        )


def annotateLD_cache(input_path, out_path, pop, ldlink_token):
    logging.info(f"Using annotateLD with python API and sqlite cache")

    df = pandas.read_csv(input_path, sep="\t", header=0)
    df = df[["chr", "pos", "geneID", "rsid", "AF"]]
    df.set_index(df["chr"] + ":" + df["pos"].map(str), inplace=True)

    pairs = []
    chrpos_to_rsid = {}
    rsid_to_chrpos = {}
    # TODO make this part faster?
    prev = df.iloc[0]
    for i in range(1, len(df)):
        cur = df.iloc[i]
        if prev.geneID == cur.geneID:
            pairs.append([cur.name, prev.name])
            chrpos_to_rsid[cur.name] = cur.rsid
            chrpos_to_rsid[prev.name] = prev.rsid
            rsid_to_chrpos[cur.rsid] = cur.name
            rsid_to_chrpos[prev.rsid] = prev.name
        prev = cur

    ldlink_infos = fetch_ldpairs(pairs, pop, ldlink_token, chrpos_to_rsid)

    df[["pair_pos", "r2", "d"]] = "NA"

    for info in ldlink_infos:
        df.loc[info.pair[0], "pair_pos"] = info.pair[1]
        df.loc[info.pair[0], "r2"] = info.r2
        df.loc[info.pair[0], "d"] = info.d

    df.reset_index()
    df.to_csv(out_path, index=False, sep="\t")


LDPairInfo = namedtuple("LDPairInfo", ["pair", "r2", "d"])


def fetch_ldpairs(pairs, pop, ldlink_token, chrpos_to_rsid):
    db = get_cache_con("test_cache.db")
    ldlink_infos = []
    pairs_to_fetch = []
    cur = db.cursor()
    for pair in pairs:
        cur.execute(
            "SELECT r2, d FROM ldpairs WHERE chrpos1 = ? AND chrpos2 = ?",
            (pair[0], pair[1]),
        )
        row = cur.fetchone()
        if not row:
            pairs_to_fetch.append(pair)
        else:
            ldlink_infos.append(LDPairInfo(pair, row[0], row[1]))
    cur.close()

    logging.debug(
        f"Got {len(ldlink_infos)} hits and {len(pairs_to_fetch)} misses from {len(pairs)} pairs"
    )

    BATCH_SIZE = 400
    batches = get_batches(pairs_to_fetch, BATCH_SIZE)
    logging.debug(f"{len(batches)} batches to fetch run for {len(pairs)} pairs")

    fetched_ldpairs = []
    for i, batch in enumerate(batches):
        logging.debug(f"Fetching batch {i+1} / {len(batches)} - {len(batch)} pairs")
        batch_pairs = fetch_ldpairs_from_api(batch, pop, ldlink_token, chrpos_to_rsid)
        fetched_ldpairs.extend(batch_pairs)
        for pair in batch_pairs:
            ldlink_infos.append(pair)

    if len(fetched_ldpairs):
        logging.debug(f"inserting {len(fetched_ldpairs)} values into cache")
        cur = db.cursor()
        cur.executemany(
            "INSERT INTO ldpairs VALUES (?, ?, ?, ?)",
            [(info.pair[0], info.pair[1], info.r2, info.d) for info in fetched_ldpairs],
        )
        cur.close()

    db.close()
    return ldlink_infos


def get_cache_con(db_path):
    con = sqlite3.connect(db_path)

    cur = con.cursor()
    cur.execute(
        "CREATE TABLE IF NOT EXISTS ldpairs (chrpos1 TEXT, chrpos2 TEXT, r2 REAL, d REAL)"
    )
    con.commit()
    cur.close()

    return con


def get_batches(pairs, pairs_per_batch):
    batches = []

    cur_batch = []
    cur_chr = None

    def finish_batch():
        nonlocal batches, cur_batch, cur_chr
        if len(cur_batch) > 0:
            batches.append(cur_batch)
        cur_batch = []
        cur_chr = None

    for pair in pairs:
        chr = pair[0].split(":")[0]
        if len(cur_batch) == pairs_per_batch or chr != cur_chr:
            finish_batch()
        cur_chr = chr
        cur_batch.append(pair)

    finish_batch()

    return batches


def unique_snps_from_pairs(pairs):
    s = set()
    ret = []
    for pair in pairs:
        for snp in pair:
            if snp not in s:
                s.add(snp)
                ret.append(snp)
    return ret


def fetch_ldpairs_from_api(pairs, pop, ldlink_token, chrpos_to_rsid):
    snps = unique_snps_from_pairs(pairs)

    url = f"https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?token={ldlink_token}"

    # TODO error handling.  This includes verifying the number of returned snps is as expected
    r2_req = requests.post(
        url, json={"snps": "\n".join(snps), "pop": pop, "r2_d": "r2"}
    )
    r2_matrix = pandas.read_csv(StringIO(r2_req.text), sep="\t", header=0, index_col=0)

    d_req = requests.post(url, json={"snps": "\n".join(snps), "pop": pop, "r2_d": "d"})
    d_matrix = pandas.read_csv(StringIO(d_req.text), sep="\t", header=0, index_col=0)

    ret = []
    for pair in pairs:
        rsid1 = chrpos_to_rsid[pair[0]]
        rsid2 = chrpos_to_rsid[pair[1]]
        if (
            rsid1 in r2_matrix
            and rsid2 in r2_matrix
            and rsid1 in d_matrix
            and rsid2 in d_matrix
        ):
            ret.append(
                LDPairInfo(pair, r2_matrix[rsid1][rsid2], d_matrix[rsid1][rsid2])
            )
        else:
            print(f"WARN: information {pair} was not returned from ldmatrix request")
            ret.append(LDPairInfo(pair, "NA", "NA"))

    return ret
