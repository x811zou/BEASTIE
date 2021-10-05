#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# 2018 William H. Majoros (bmajoros@allumni.duke.edu)
#=========================================================================
from .Rex import Rex
from .SamMDtagParser import SamMDtagParser

rex=Rex()

#=========================================================================
# Attributes:
#   ID = read identifier
#   refName = name of reference sequence the read aligns to
#   refPos = position in reference where alignment begins
#   CIGAR = CigarString
#   seq = read sequence
#   flags = bitfield
#   tags = array of tags at end of record (MD:Z:122G25, NM:i:1, etc.)
# Instance Methods:
#   rec=SamRecord(ID,refName,refPos,cigar,seq,flags,tags)
#   ID=rec.getID()
#   cigar=rec.getCigar() # returns CigarString object
#   seq=rec.getSequence()
#   L=rec.seqLength()
#   refName=rec.getRefName()
#   refPos=rec.getRefPos()
#   tags=rec.getTags()
#   fields=rec.parseMDtag()
#   N=rec.countMismatches() # uses MD tag
#   tag=getTag("MD") # returns the third field, e.g. "122G25" in MD:Z:122G25
#   bool=rec.flag_hasMultipleSegments()
#   bool=rec.flag_properlyAligned()
#   bool=rec.flag_unmapped()
#   bool=rec.flag_nextSegmentUnmapped()
#   bool=rec.flag_revComp()
#   bool=rec.flag_nextSegmentRevComp()
#   bool=rec.flag_firstOfPair()
#   bool=rec.flag_secondOfPair()
#   bool=rec.flag_secondaryAlignment()
#   bool=rec.flag_failedFilters()
#   bool=rec.flag_PCRduplicate()
#   bool=rec.flag_supplAlignment()
# Class Methods:
#=========================================================================
class SamRecord:
    """SamRecord"""
    def __init__(self,ID,refName,refPos,CIGAR,seq,flags,tags):
        self.ID=ID
        self.refName=refName
        self.refPos=refPos
        self.CIGAR=CIGAR
        self.seq=seq
        self.flags=flags
        self.tags=tags

    def seqLength(self):
        return len(self.seq)

    def getTags(self):
        return self.tags

    def getTag(self,label):
        for tag in self.tags:
            if(not rex.find("^([^:]+):[^:]+:(\S+)",tag)):
                raise Exception("Can't parse SAM tag: "+tag)
            if(rex[1]==label): return rex[2]
        return None

    def countMismatches(self):
        N=SamMDtagParser.countMismatches(self.parseMDtag())
        return N

    def parseMDtag(self):
        md=self.getTag("MD")
        fields=[]
        if(md is None): return None
        while(len(md)>0):
            if(rex.find("^(\d+)(.*)",md)):
                fields.append(rex[1])
                md=rex[2]
            elif(rex.find("^([ACGTN])(.*)",md)):
                fields.append(rex[1])
                md=rex[2]
            elif(rex.find("^(\^[ACGTN]+)(.*)",md)):
                fields.append(rex[1])
                md=rex[2]
            else:
                raise Exception("Can't parse MD tag: "+md)
        return fields

    def getRefName(self):
        return self.refName

    def getRefPos(self):
        return self.refPos

    def getCigar(self):
        return self.CIGAR

    def getID(self):
        return self.ID

    def getSequence(self):
        return self.seq

    def flag_hasMultipleSegments(self):
        return bool(self.flags & 0x1)

    def flag_properlyAligned(self):
        return bool(self.flags & 0x2)

    def flag_unmapped(self):
        return bool(self.flags & 0x4)

    def flag_nextSegmentUnmapped(self):
        return bool(self.flags & 0x8)

    def flag_revComp(self):
        return bool(self.flags & 0x10)

    def flag_nextSegmentRevComp(self):
        return bool(self.flags & 0x20)

    def flag_firstOfPair(self):
        return bool(self.flags & 0x40)

    def flag_secondOfPair(self):
        return bool(self.flags & 0x80)

    def flag_secondaryAlignment(self):
        return bool(self.flags & 0x100)

    def flag_failedFilters(self):
        return bool(self.flags & 0x200)

    def flag_PCRduplicate(self):
        return bool(self.flags & 0x400)

    def flag_supplAlignment(self):
        return bool(self.flags & 0x800)

# FLAGS:
#   0x1 template having multiple segments in sequencing
#   0x2 each segment properly aligned according to the aligner
# > 0x4 segment unmapped
# > 0x8 next segment in the template unmapped
# > 0x10 SEQ being reverse complemented
#   0x20 SEQ of the next segment in the template being reverse complemented
# > 0x40 the first segment in the template
# > 0x80 the last segment in the template
#   0x100 secondary alignment
#   0x200 not passing filters, such as platform/vendor quality controls
# > 0x400 PCR or optical duplicate
#   0x800 supplementary alignment

# M03884:303:000000000-C4RM6:1:1101:1776:15706    99      chrX:31786371-31797409  6687    44      150M    =       6813    271     ATACTATTGCTGCGGTAATAACTGTAACTGCAGTTACTATTTAGTGATTTGTATGTAGATGTAGATGTAGTCTATGTCAGACACTATGCTGAGCATTTTATGGTTGCTATGTACTGATACATACAGAAACAAGAGGTACGTTCTTTTACA  BBBBFFFFFFFGGGGGEFGGFGHFHFFFHHHFFHHHFHFHHHGFHEDGGHFHBGFHGBDHFHFFFHHHHFHHHHHGHGFFBGGGHFHFFHHFFFFHHHHGHGFHHGFHGHHHGFHFFHHFHHFFGFFFFGGEHFFEHHFGHHHGHHHHFB  AS:i:300        XN:i:0

