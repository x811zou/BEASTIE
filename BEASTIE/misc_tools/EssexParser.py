#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import os

from .EssexNode import EssexNode
from .EssexScanner import EssexScanner

######################################################################
#
# A parser for the Essex language, which is a simple alternative to XML
# and looks like LISP programs ("S-Expressions" -- hence the name "Essex").
# It generates tree data structures built out of EssexNode objects.
#
# Attributes:
#   file : file handle
#   isOpen : boolean
#   scanner : EssexScanner
# Methods:
#   parser=EssexParser(filename) # filename is optional
#   parser.open(filename) # unnecessary if you gave filename to ctor
#   parser.close()
#   tree=parser.nextElem()   # returns root of the tree
#   forest=parser.parseAll() # returns an array of trees
# Class methods:
#   forest=EssexParser.loadFile(filename) # returns array of trees
######################################################################

class EssexParser:
    def __init__(self,filename=None):
        self.isOpen=False
        if(filename): self.open(filename)

    def open(self,filename):
        if(self.isOpen): self.close()
        if(not os.path.exists(filename)):
            raise Exception(filename+" does not exist")
        file=self.file=open(filename,"r")
        self.isOpen=True
        self.scanner=EssexScanner(file)

    def close(self):
        if(self.isOpen):
            self.file.close()
            self.isOpen=False
            self.scanner=None

    def parseAll(self):
        forest=[]
        while(True):
            tree=self.nextElem()
            if(not tree): break
            forest.append(tree)
        return forest

    @classmethod
    def loadFile(cls,filename):
        parser=EssexParser(filename)
        return parser.parseAll()

    def nextElem(self):
        if(not self.isOpen): raise Exception("file is not open")
        scanner=self.scanner
        token=scanner.nextToken()
        if(not token): return None
        if(token.isOpenParen()): return self.parseTuple()
        elif(token.isLiteral()): return token.getLexeme()
        else:
            lexeme=token.getLexeme()
            exit("Syntax error near \n"+token+"\n")

    def parseTuple(self):
        """PRECONDITION: a "(" has already been matched"""
        elements=[]
        scanner=self.scanner
        while(True):
            token=scanner.nextToken()
            if(not token): break
            if(token.isOpenParen()):
                elements.append(self.parseTuple())
            elif(token.isCloseParen()): break
            else: elements.append(token.getLexeme())
        return EssexNode(elements)

