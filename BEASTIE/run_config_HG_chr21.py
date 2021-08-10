#!/usr/bin/env python

"""
usage: python run_config_HG_chr21.py
"""
import BEASTIE
if __name__ == "__main__":
    config = BEASTIE.load_configuration("parameters_HG_chr21.cfg")
    BEASTIE.run(config)
