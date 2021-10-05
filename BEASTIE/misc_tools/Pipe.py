#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import subprocess
import sys


#=========================================================================
# Attributes:
#   iter
# Instance Methods:
#   pipe=Pipe(command)
#   line=pipe.readline()
# Class Methods:
#   output=Pipe.run(command)
#=========================================================================
class Pipe:
    """Pipe reads output from a shell command"""
    def __init__(self,cmd):
        p=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True,
                           stderr=subprocess.STDOUT)
        self.iter=iter(p.stdout.readline, b'')

    def readline(self):
        line=next(self.iter,None)
        if(line is None): return None
        line=line.decode(sys.stdout.encoding).rstrip()
        return line

    @classmethod
    def run(cls,cmd):
        line=subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)
        line=line.decode(sys.stdout.encoding).rstrip()
        return line

