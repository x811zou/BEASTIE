#!/usr/bin/env python

"""
usage: python parameters_HG00097.py
"""
import BEASTIE

if __name__ == "__main__":
    config = BEASTIE.load_configuration("parameters_HG00097.cfg")
    BEASTIE.run(config)
