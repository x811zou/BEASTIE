#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .Pipe import Pipe

#import subprocess
#import sys

#def run_command(command):
#    command=command.split()
#    p = subprocess.Popen(command,
#                         stdout=subprocess.PIPE,
#                         stderr=subprocess.STDOUT)
#    return iter(p.stdout.readline, b'')

#for line in run_command("ls -la"):
#    line=line.decode(sys.stdout.encoding).rstrip()
#    #line=line.decode("utf-8").rstrip()
#    print("line=",line)

pipe=Pipe("ls -la | sort -g | sort -r | uniq -c")
while(True):
    line=pipe.readline()
    if(not line): break
    print(line)

