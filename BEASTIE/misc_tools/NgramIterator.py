#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
#=========================================================================

######################################################################
#
# NgramIterator.py bmajoros
#
# Attributes:
#   array ngram : array of integer indices into alphabet string
#   string alphabet
# Methods:
#   ngramIterator=NgramIterator("ATCG",N)
#   string=ngramIterator.nextString() # returns None if no more
#   ngramIterator.reset()
# Private methods:
#   self.ngramToString()
######################################################################

class NgramIterator:
#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
# ngramIterator=NgramIterator("ATCG",N)
    def __init__(self,alphabet,N):
        ngram=[]
        alphaSize=len(alphabet)
        for i in range(N-1): ngram.append(0)
        if(N>0): ngram.append(-1)
        self.alphabet=alphabet
        self.ngram=ngram
#---------------------------------------------------------------------
# ngramIterator.reset()
    def reset(self):
        ngram=self.ngram
        L=len(ngram)
        for i in range(L-1): ngram[i]=0
        ngram[L-1]=-1
#---------------------------------------------------------------------
# string=ngramIterator.nextString() # returns undef if no more
    def nextString(self):
        alphabet=self.alphabet
        ngram=self.ngram
        if(ngram is None): return None
        L=len(ngram)
        if(L==0):
            self.ngram=None
            return ""
        alphaSize=len(alphabet)
        for i in range(L-1,-1,-1): # for($i=$len-1 ; $i>=0 ; --$i)
            ngram[i]+=1
            index=ngram[i]
            if(index<alphaSize): return self.ngramToString()
            ngram[i]=0
        return None
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
# self.ngramToString()
    def ngramToString(self):
        ngram=self.ngram
        alphabet=self.alphabet
        length=len(ngram)
        string=""
        for i in range(length):
            index=ngram[i]
            string+=alphabet[index:index+1]
        return string


