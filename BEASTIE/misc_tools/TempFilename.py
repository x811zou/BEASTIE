#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import os
import tempfile


def generate(suffix=""):
  [fh,filename]=tempfile.mkstemp(suffix)
  os.close(fh)
  return filename


