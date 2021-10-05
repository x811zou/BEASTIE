#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================

class Integer(object):
    """Integer is a first-class mutable object, so it can be passed
    by reference, allowing functions to modify its value.
    """
    def __init__(self,x=0):
        self.value=int(x)
    def __int__(self):
        return self.value
    def __add__(self,x):
        return Integer(self.value+int(x))
    def __sub__(self,x):
        return Integer(self.value-int(x))
    def __mul__(self,x):
        return Integer(self.value*int(x))
    def __floordiv__(self,x):
        return Integer(int(self.value/int(x)))
    def __div__(self,x):
        newVal=float(self.value)/float(x)
        if(newVal==float(int(newVal))): return Integer(newVal)
        return newVal
    def __truediv__(self,x):
        newVal=float(self.value)/float(x)
        if(newVal==float(int(newVal))): return Integer(newVal)
        return newVal
    def __mod__(self,x):
        return Integer(self.value%int(x))
    def __pow__(self,x): # **
        return Integer(self.value**int(x))
    def __lshift__(self,x): # <<
        return Integer(self.value<<int(x))
    def __rshift__(self,x): # >>
        return Integer(self.value>>int(x))
    def __and__(self,x): # &
        return Integer(self.value & int(x))
    def __xor__(self,x): # ^
        return Integer(self.value ^ int(x))
    def __or__(self,x): # |
        return Integer(self.value | int(x))
    def __iadd__(self,x):
        self.value+=int(x)
        return self
    def __isub__(self,x): # -=
        self.value-=int(x)
        return self
    def __imul__(self,x): # *=
        self.value*=int(x)
        return self
    def __idiv__(self,x): # /=
        self.value=int(self.value/int(x))
        return self
    def __ifloordiv__(self,x): # //=
        self.value=int(self.value//int(x))
        return self
    def __imod__(self,x): # %=
        self.value%=int(x)
        return self
    def __ipow__(self,x): # **=
        self.value**=int(x)
        return self
    def __ilshift__(self,x): # <<=
        self.value<<=int(x)
        return self
    def __irshift__(self,x): # >>=
        self.value>>=int(x)
        return self
    def __iand__(self,x): # &=
        self.value &= int(x)
        return self
    def __ixor__(self,x): # ^=
        self.value ^= int(x)
        return self
    def __ior__(self,x): # |=
        self.value |= int(x)
        return self
    def __neg__(self): # -
        return Integer(-self.value)
    def __pos__(self): # +
        return Integer(+self.value)
    def __abs__(self): # abs()
        return Integer(abs(self.value))
    def __invert__(self): # ~
        return Integer(~self.value)
    def __complex__(self): # complex()
        return complex(self.value)
    def __long__(self): # long()
        return long(self.value)
    def __float__(self): # float()
        return float(self.value)
    def __oct__(self): # oct()
        return oct(self.value)
    def __hex__(self): # hex()
        return hex(self.value)
    def __lt__(self,x): # <
        if(type(x)==float): return self.value<x
        return self.value<int(x)
    def __le__(self,x): # <=
        if(type(x)==float): return self.value<=x
        return self.value<=x
    def __eq__(self,x): # ==
        if(type(x)==float): return self.value==x
        return self.value==int(x)
    def __ne__(self,x): # !=
        if(type(x)==float): return self.value!=x
        return self.value!=int(x)
    def __ge__(self,x): # >=
        if(type(x)==float): return self.value>=x
        return self.value>=float(x)
    def __gt__(self,x): # >
        if(type(x)==float): return self.value>x
        return self.value>int(x)
    def __str__(self):
        return str(self.value)
    def __repr__(self):
        return "Integer("+str(self.value)+")"

#=========================== TESTING CODE =========================
def testInteger_func(x):
    x+=2

def testInteger():
    x=Integer(1)
    y=Integer(3)
    x+=y+5
    testInteger_func(x)
    print(x)
    print(x>y,x<y,x==y,x>=y,x<=y,x!=y)
    y=11
    print(x>y,x<y,x==y,x>=y,x<=y,x!=y)
    y=Integer(11)
    print(x>y,x<y,x==y,x>=y,x<=y,x!=y)
    print(float(x))
    y=Integer(3)
    print(x,"%",y,"=",x%y)
    y=3
    print(x,"%",y,"=",x%y)
    print((x+4-y*2)/x)

#testInteger()

