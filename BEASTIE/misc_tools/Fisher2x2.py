#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from .Pipe import Pipe


#=========================================================================
# Attributes:
#   Matrix of the form:
#       x00  x01
#       x10  x11
# Instance Methods:
#   fisher=Fisher2x2(x00,x01,x10,x11)
#   P=fisher.getPvalue()
#   (exp00,exp01,exp10,exp11)=fisher.getExpectedCounts()
# Class Methods:
#
#=========================================================================
class Fisher2x2:
    """Fisher2x2 performs Fisher's exact test for 2x2 contingency tables"""
    def __init__(self,x00,x01,x10,x11):
        self.x00=x00
        self.x01=x01
        self.x10=x10
        self.x11=x11

    def getPvalue(self):
        executable=Pipe.run("which fisher-exact-test.R")
        cmd=executable+" "+str(self.x00)+" "+str(self.x01)+" "+\
            str(self.x10)+" "+str(self.x11)
        P=float(Pipe.run(cmd))
        return P

    def getExpectedCounts(self):
        x00=float(self.x00); x01=float(self.x01)
        x10=float(self.x10); x11=float(self.x11)
        N=x00+x01+x10+x11
        pTop=(x00+x01)/N
        pBottom=1.0-pTop
        leftSum=x00+x10
        rightSum=x01+x11
        exp00=int(round(pTop*leftSum,0))
        exp01=int(round(pTop*rightSum,0))
        exp10=int(round(pBottom*leftSum,0))
        exp11=int(round(pBottom*rightSum,0))
        return (exp00,exp01,exp10,exp11)



