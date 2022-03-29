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


class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
        self.encoding = self.stdout.encoding

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()
        self.stdout.flush()

    def fileno(self):
        return self.stdout.fileno()

    def __enter__(self):
        pass

    def __exit__(self, _type, _value, _traceback):
        sys.stdout = self.stdout
        self.file.close()
