#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import os
import sys

from . import ProgramName

# Process command line
name=ProgramName.get();
if(len(sys.argv)!=2):
    sys.exit(name+" <classname>")
className=sys.argv[1]

# Write file
filename=className+".py"
if(os.path.exists(filename)):
    sys.exit(filename+" exists");
fh=open(filename,"w")
header="\n".join(["#=========================================================================",
    "# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public",
    "# License (GPL) version 3, as described at www.opensource.org.",
    "# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).",
    "#=========================================================================",
    "from __future__ import (absolute_import, division, print_function,",
    "   unicode_literals, generators, nested_scopes, with_statement)",
    "from builtins import (bytes, dict, int, list, object, range, str, ascii,",
    "   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)",
    "\n"
    "#=========================================================================",
    "# Attributes:",
    "#   ",
    "# Instance Methods:",
    "#   "+className+"()",
    "# Class Methods:",
    "#   ",
    "#=========================================================================\n",
    ])
fh.write(header)
fh.write("class "+className+":\n")
fh.write("    \"\"\""+className+"\"\"\"\n");
fh.write("    def __init__(self):\n")
fh.write("        pass\n")
fh.write("\n\n\n")
fh.close()


