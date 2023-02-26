#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
from collections import namedtuple
from io import StringIO
import logging
import os
import sqlite3
from pathlib import Path
import sys
import time
from contextlib import contextmanager

import pandas
import requests

from BEASTIE.ldlink_token_db import acquire_ldlink_token


def annotateLD(
    ancestry, input_path, ldlink_token, out_path, ldlink_cache_dir, ldlink_token_db
):
    logging.info(f"Using annotateLD with python API and sqlite cache")

    df = pandas.read_csv(input_path, sep="\t", header=0)
    df = df[["chr", "pos", "geneID", "rsid", "AF"]]
    df.set_index(df["chr"] + ":" + df["pos"].map(str), inplace=True)

    pairs = []
    chrpos_to_rsid = {}
    # TODO make this part faster?
    prev = df.iloc[0]
    for i in range(1, len(df)):
        cur = df.iloc[i]
        if prev.geneID == cur.geneID:
            pairs.append([cur.name, prev.name])
            chrpos_to_rsid[cur.name] = cur.rsid
            chrpos_to_rsid[prev.name] = prev.rsid
        prev = cur

    with acquire_ldlink_token(ldlink_token, ldlink_token_db) as acquired_token:
        try:
            ldlink_infos = fetch_ldpairs(
                pairs,
                ancestry,
                chrpos_to_rsid,
                acquired_token,
                ldlink_cache_dir,
            )
        except Exception as e:
            logging.error(f"Error fetching LD pairs: {e}")
            raise e

    for info in ldlink_infos:
        df.loc[info.pair[0], "pair_pos"] = info.pair[1].split(":")[1]
        df.loc[info.pair[0], "r2"] = info.r2
        df.loc[info.pair[0], "d"] = info.d

    df.to_csv(out_path, index=False, sep="\t")


LDPairInfo = namedtuple("LDPairInfo", ["pair", "r2", "d"])


def pair_key(pair, pop):
    return f"{pop}_{pair[0]}_{pair[1]}"


def fetch_ldpairs(pairs, pop, chrpos_to_rsid, ldlink_token, cache_dir):
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    db = get_cache_con(os.path.join(cache_dir, "ldlink_cache.db"))
    ldlink_infos = []
    pairs_to_fetch = []
    cur = db.cursor()
    in_clause = ",".join([f"'{pair_key(pair, pop)}'" for pair in pairs])
    rows = cur.execute(
        f"SELECT key, r2, d FROM ldpairs2 WHERE key IN ({in_clause})",
    ).fetchall()
    data_by_key = {row[0]: row[1:] for row in rows}
    for pair in pairs:
        row = data_by_key.get(pair_key(pair, pop))
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
        with db:
            cur = db.cursor()
            cur.executemany(
                "INSERT OR IGNORE INTO ldpairs2 VALUES (?, ?, ?)",
                [
                    (pair_key(info.pair, pop), info.r2, info.d)
                    for info in fetched_ldpairs
                ],
            )
            logging.debug(f"Inserted {cur.rowcount} rows into cache")

    db.close()
    return ldlink_infos


def get_cache_con(db_path):
    logging.debug(f"Loading LDLink cache from {db_path}")
    con = sqlite3.connect(db_path, timeout=120)

    with con:
        cur = con.cursor()

        cur.execute(
            """
            SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='ldpairs2';
        """
        )
        r = cur.fetchone()
        if r[0] == 0:
            logging.debug(f"Creating LDLink cache table ldpairs2")
            cur.execute(
                """
        CREATE TABLE IF NOT EXISTS ldpairs2 (
            key TEXT NOT NULL PRIMARY KEY,
            r2 REAL,
            d REAL
        );"""
            )
        cur.close()

    logging.debug(f"Loaded LDLink cache from {db_path}")

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
