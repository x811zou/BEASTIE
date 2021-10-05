#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import os

from . import TempFilename
from .CigarString import CigarString
from .FastaWriter import FastaWriter
from .Pipe import Pipe
from .Rex import Rex

rex=Rex()
#=========================================================================
# Attributes:
#
# Instance Methods:
#   aligner=SmithWaterman(alignerDir,matrixFile,openPenalty,extrendPenalty)
#   cigarString=aligner.align(seq1,seq2)
#=========================================================================

class SmithWaterman:
    def __init__(self,alignerDir,matrixFile,gapOpenPenalty,gapExtendPenalty):
        self.alignerDir=alignerDir
        self.matrixFile=matrixFile
        self.gapOpen=gapOpenPenalty
        self.gapExtend=gapExtendPenalty
        self.fastaWriter=FastaWriter()

    def writeFile(self,defline,seq):
        filename=TempFilename.generate("fasta")
        self.fastaWriter.writeFasta(defline,seq,filename)
        return filename

    def swapInsDel(self,cigar):
        # This is done because my aligner defines insertions and deletions
        # opposite to how they're defined in the SAM specification
        newCigar=""
        for x in cigar:
            if(x=="I"): x="D"
            elif(x=="D"): x="I"
            newCigar+=x
        return newCigar

    def align(self,seq1,seq2):
        file1=self.writeFile("query",seq1)
        file2=self.writeFile("reference",seq2)
        cmd=self.alignerDir+"/smith-waterman -q "+self.matrixFile+" "+\
            str(self.gapOpen)+" "+str(self.gapExtend)+" "+file1+" "+file2+" DNA"
        output=Pipe.run(cmd)
        os.remove(file1)
        os.remove(file2)
        if(not rex.find("CIGAR=(\S+)",output)):
            raise Exception("Can't parse aligner output: "+output)
        cigar=rex[1]
        cigar=self.swapInsDel(cigar) # because I define cigars differently
        return CigarString(cigar)



