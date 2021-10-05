#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# 2018 William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
ADVANCE_QUERY=set(["M","I","S","H","=","X"])
ADVANCE_REF=set(["M","D","N","=","X"])

#=========================================================================
# Attributes:
#   length : integer
#   interval1 : Interval (in sequence 1 = query)
#   interval2 : Interval (in sequence 2 = reference)
#   op : M(or =/X)/I/D/S:
#                                                              consumes
#                                                              query ref
#     M 0 alignment match (can be a sequence match or mismatch) yes yes
#     I 1 insertion to the reference                            yes  no
#     D 2 deletion from the reference                            no yes
#     N 3 skipped region from the reference                      no yes
#     S 4 soft clipping (clipped sequences present in SEQ)      yes  no
#     H 5 hard clipping (clipped sequences NOT present in SEQ)   no  no
#     P 6 padding (silent deletion from padded reference)        no  no
#     = 7 sequence match                                        yes yes
#     X 8 sequence mismatch                                     yes yes
# Instance Methods:
#   op=CigarOp("M",135)
#   bool=op.advanceInQuery()   # matches, insertions, etc.
#   bool=op.advanceInRef() # matches, deletions, etc.
#   op=op.getOp()
#   L=op.getLength()
#   interval=op.getQueryInterval() # sequence 1
#   interval=op.getRefInterval() # sequence 2
#=========================================================================
class CigarOp:
    def __init__(self,op,L):
        self.op=op
        self.length=L
        self.interval1=None
        self.interval2=None

    def getQueryInterval(self):
        return self.interval1

    def getRefInterval(self):
        return self.interval2

    def getOp(self):
        return self.op

    def getLength(self):
        return self.length

    def advanceInQuery(self):
        return self.op in ADVANCE_QUERY

    def advanceInRef(self):
        return self.op in ADVANCE_REF

