#!/usr/bin/env python

"""
usage: python run_config_HG00250.py
"""
import BEASTIE
if __name__ == "__main__":
    config = BEASTIE.load_configuration("parameters_HG00250.cfg")
    BEASTIE.run(config)
