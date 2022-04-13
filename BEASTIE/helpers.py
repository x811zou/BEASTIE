#!/usr/bin/env python
# =========================================================================
# Copyright (C) Xue Zou (xue.zou@duke.edu)
# =========================================================================
import logging
import subprocess
import sys
import tempfile

from BEASTIE.misc_tools.Pipe import Pipe


def runhelper(cmd):
    subprocess.run(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)


def flatten(lst):
    return [x for sublst in lst for x in sublst]


# Returns iterator for chromosomes from start to end inclusive
# and optionally the x chromosomee
def chrRange(start, end, include_x_chromosome):
    for i in range(start, end + 1):
        yield str(i)
    if include_x_chromosome:
        yield "X"


def tabix_regions(regions, line_processor, target_file_path, comment_char="#"):
    region_to_results = {}

    logging.info(
        f"Start tabix extraction of {len(regions)} regions from file {target_file_path}"
    )

    if len(regions) > 1000:
        with tempfile.NamedTemporaryFile(mode="w") as file:
            for region in regions:
                chr, rest = region.split(":")
                start, end = rest.split("-")
                file.write(f"{chr}\t{start}\t{end}\n")
            file.flush()
            command = (
                f"tabix --separate-regions {target_file_path} --regions {file.name}"
            )

            output = Pipe.run(command)
    else:
        region_batch = " ".join(regions)
        command = f"tabix --separate-regions {target_file_path} {region_batch}"
        output = Pipe.run(command)

    if len(output) == 0:
        return region_to_results

    lines = output.split("\n")
    records = []
    for line in lines:
        if len(line) == 0:
            continue
        # start accumulating new region
        if line.startswith(comment_char):
            region_str = line[1:]
            records = []
            region_to_results[region_str] = records
            continue

        result = line_processor(line)
        if result is not None:
            records.append(result)

    logging.info(f"Got {len(region_to_results)} / {len(regions)} regions with data")

    return region_to_results
