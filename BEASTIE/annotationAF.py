#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
from collections import namedtuple
import csv
import gzip
import logging
import multiprocessing
import os
import pandas as pd
from pkg_resources import resource_filename


def annotateAF(af_path, ancestry, hetSNP, out_AF):
    ancestry_map = {
        "EUR": "EUR",
        "CEU": "EUR",
        "TSI": "EUR",
        "FIN": "EUR",
        "GBR": "EUR",
        "IBS": "EUR",
        "AFR": "AFR",
        "YRI": "AFR",
        "LWK": "AFR",
        "GWD": "AFR",
        "MSL": "AFR",
        "ESN": "AFR",
        "ASW": "AFR",
        "ACB": "AFR",
        "SAS": "SAS",
        "GIH": "SAS",
        "PJL": "SAS",
        "BEB": "SAS",
        "STU": "SAS",
        "ITU": "SAS",
        "EAS": "EAS",
        "CHB": "EAS",
        "JPT": "EAS",
        "CHS": "EAS",
        "CDX": "EAS",
        "KHV": "EAS",
        "AMR": "AMR",
        "MXL": "AMR",
        "PUR": "AMR",
        "CLM": "AMR",
        "PEL": "AMR",
    }

    ancestry_group = ancestry_map[ancestry]

    annotateAFConstantMemoryParallel(af_path, ancestry_group, hetSNP, out_AF)

    logging.info("..... finish annotating AF for SNPs, file save at {0}".format(out_AF))


def annotateCHRLines(
    hetSNP_lines, AF_filename, ancestry, hetSNP_chr_index, hetSNP_pos_index
):
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


def annotateAFConstantMemoryParallel(af_path, ancestry, hetSNP, out_AF):
    logging.info("..... start annotating AF with parallel constant memory join")

    with multiprocessing.Pool() as pool, open(hetSNP, newline="") as hetSNPfile, open(
        out_AF, "w", newline=""
    ) as outFile:
        hetSNP_reader = csv.reader(hetSNPfile, delimiter="\t", dialect="unix")
        hetSNP_header = next(hetSNP_reader)

        hetSNP_chr_index = hetSNP_header.index("chr")
        hetSNP_pos_index = hetSNP_header.index("pos")

        output_handles = []

        current_chr = None
        current_rows = None
        current_af_path = None

        for row in hetSNP_reader:
            chr = row[hetSNP_chr_index]
            if chr != current_chr:
                if current_chr is not None:
                    output_handles.append(
                        pool.apply_async(
                            annotateCHRLines,
                            (
                                current_rows,
                                current_af_path,
                                ancestry,
                                hetSNP_chr_index,
                                hetSNP_pos_index,
                            ),
                        )
                    )
                current_chr = chr
                current_rows = []
                current_af_path = os.path.join(af_path, f"AF_{chr}.csv.gz")
            current_rows.append(row)
        if current_rows is not None and len(current_rows) > 0:
            output_handles.append(
                pool.apply_async(
                    annotateCHRLines,
                    (
                        current_rows,
                        current_af_path,
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
