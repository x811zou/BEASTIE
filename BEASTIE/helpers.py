
import subprocess
import sys

def runhelper(cmd):
  subprocess.run(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)