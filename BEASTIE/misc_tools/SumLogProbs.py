#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import math

NEGATIVE_INFINITY=float("-inf")


def sumLogProbs_ordered(largerValue,smallerValue):
    if(smallerValue==NEGATIVE_INFINITY): return largerValue
    return largerValue+math.log(1+math.exp(smallerValue-largerValue))



def sumLogProbs2(logP,logQ):
    if(logP>logQ): return sumLogProbs_ordered(logP,logQ)
    return sumLogProbs_ordered(logQ,logP)



def sumLogProbs(x):
    n=len(x)
    if(n==1): return x[0]
    if(n==0): return NEGATIVE_INFINITY

    # Pull out the largest value
    largestValue=x[0]
    for v in x:
        if(v>largestValue): largestValue=v

    # Handle the case of all zeros separately
    if(largestValue==NEGATIVE_INFINITY): return NEGATIVE_INFINITY

    # Apply the Kingsbury-Raynor formula
    sum=0.0;
    for v in x:
        if(v==NEGATIVE_INFINITY): continue
        sum+=math.exp(v-largestValue)
    return largestValue+math.log(sum)


#v=[math.log(0.0001),math.log(0.0002),math.log(0.0003)]
#print(math.exp(sumLogProbs2(math.log(0.1),math.log(0.2))))
#print(math.exp(sumLogProbs(v)))


