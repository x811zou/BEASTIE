#====================================================================
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURxCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
#====================================================================
import re

######################################################################
# Attributes:
#   hash
# Methods:
#   configFile=ConfigFile(filename)
#   value=configFile.lookup(key)
#   value=configFile.lookupOrDie(key)
# Private methods:
#   self.load(filename)
######################################################################

class ConfigFile:
    """ConfigFile stores variable assignments as key-value pairs
    in a hash
    """
    def __init__(self,filename):
        self.hash={}
        self.load(filename)

    def lookup(self,key):
        return self.hash[key]

    def lookupOrDie(self,key):
        if(self.hash[key] is None):
            raise Exception("$key not defined in config file\n")
        return self.hash[key]

    def load(self,filename):
        hash=self.hash
        with open(filename,"r") as fh:
            while(True):
                line=fh.readline()
                if(not line): break
                match=re.search("^(.*)#",line);
                if(match): line=match.group(1)
                match=re.search("^\s*(\S+)\s*=\s*(\S+)",line);
                if(match):
                    key=match.group(1)
                    value=match.group(2)
                    hash[key]=value

