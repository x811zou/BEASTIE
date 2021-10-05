#!/usr/bin/env python

"""
usage: python run_config_HG00096.py
"""
import BEASTIE

if __name__ == "__main__":
    config = BEASTIE.load_configuration("parameters_HG00096.cfg")
    BEASTIE.run(config)
