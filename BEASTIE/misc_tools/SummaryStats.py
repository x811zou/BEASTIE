#====================================================================
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURxCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
#====================================================================
import math

######################################################################
# Attributes:
#
# Methods:
#    [mean,SD,min,max]=SummaryStats.summaryStats(array)
#    [mean,SD,min,max]=SummaryStats.roundedSummaryStats(array)
#    sum=SummaryStats.sum(array)
#    r=SummaryStats.correlation(array1,array2)
#    m=SummaryStats.median(array)
#    array=SummaryStats.getQuantiles(values,numQuantiles)
#    (mean,SD,Min,Max)=SummaryStats.trimmedStats(array,percent)
######################################################################

class SummaryStats:
    """SummaryStats computes simple mean, variance, and correlation
       statistics
    """

    @classmethod
    def median(self,array):
        a=[]
        for x in array: a.append(x)
        a.sort()
        n=len(a)
        if(n<1): raise Exception("median is undefined for 0 elements")
        halfN=int(n/2)
        if(n%2==1): return a[halfN]
        return (a[halfN-1]+a[halfN])/2

    @classmethod
    def getQuantiles(self,values,numQuantiles):
        a=[]
        for x in values: a.append(x)
        a.sort()
        n=len(a)
        q=[0]
        index=0
        for i in range(1,numQuantiles):
            index=int(float(i)/float(numQuantiles)*float(n))
            q.append(a[index])
        q.append(a[n-1])
        return q

    @classmethod
    def trimmedStats(self,array,percent):
        sorted=[x for x in array]
        sorted.sort()
        n=len(array)
        keep=percent*n
        omit=n-keep
        first=int(omit/2)
        keep=int(keep)
        sorted=sorted[first:(first+keep)]
        return SummaryStats.summaryStats(sorted)

    @classmethod
    def summaryStats(self,array):
        n=len(array)
        minX=None
        maxX=None
        sumX=0
        sumXX=0
        for i in range(0,n):
            x=array[i]
            sumX+=x
            sumXX+=x*x
            if(i==0): minX=maxX=x
            if(x<minX): minX=x
            if(x>maxX): maxX=x
        meanX=sumX/n
        varX=None if n<2 else (sumXX-sumX*sumX/n)/(n-1)
        if(varX is not None and varX<0): varX=0
        stddevX=math.sqrt(varX) if varX is not None else None
        return [meanX,stddevX,minX,maxX]

    @classmethod
    def roundedSummaryStats(self,array):
        [mean,stddev,min,max]=SummaryStats.summaryStats(array)
        mean=int(100.0*mean+5.0/9.0)/100.0
        stddev=int(100.0*stddev+5.0/9.0)/100.0
        min=int(100.0*min+5.0/9.0)/100.0
        max=int(100.0*max+5.0/9.0)/100.0
        return [mean,stddev,min,max]

    @classmethod
    def sum(self,array):
        s=0.0
        n=len(array)
        for i in range(0,n):
            s+=array[i]
        return s;

    @classmethod
    def correlation(self,Xs,Ys):
        sumX=0.0
        sumY=0.0
        sumXY=0.0
        sumXX=0.0
        sumYY=0.0
        n=len(Xs)
        for i in range(0,n):
            x=Xs[i]
            y=Ys[i]
            sumX+=x
            sumY+=y
            sumXY+=x*y
            sumXX+=x*x
            sumYY+=y*y
        r=(sumXY-sumX*sumY/n)/math.sqrt((sumXX-sumX*sumX/n)*(sumYY-sumY*sumY/n))
        return r

