
import subprocess
import sys

def runhelper(cmd):
  subprocess.run(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)

def flatten(lst):
  return [x for sublst in lst for x in sublst]
