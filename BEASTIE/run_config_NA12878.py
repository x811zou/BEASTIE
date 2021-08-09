#!/usr/bin/env python

###############################################
# usage: python run_config_NA12878.py
###############################################
import BEASTIE
if __name__ == "__main__":
    config = BEASTIE.load_configuration("parameters_NA12878.cfg")
    BEASTIE.run(config)
