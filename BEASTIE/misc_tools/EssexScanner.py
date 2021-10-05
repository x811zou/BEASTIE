#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import re

from .EssexToken import EssexToken

######################################################################
#
# A token scanner for the EssexParser.
#
# Attributes:
#   file : file handle
#   ungot : char
#   nextTok : EssexToken
# Public Methods:
#   scanner=EssexScanner(filehandle)
#   scanner.close()
#   token=scanner.nextToken()
#   token=scanner.match(tokenType)
#   token=scanner.peek()
# Private Methods:
#   getChar
#   unGetChar
#   skipWhitespace
######################################################################

class EssexScanner:
    def __init__(self,file):
        self.file=file
        self.nextTok=None
        self.ungot=None

    def close(self):
        self.file.close()

    def nextToken(self):
        token=self.peek()
        self.nextTok=None
        return token

    def match(self,tokenType):
        token=self.nextToken()
        if(token.getType()!=tokenType):
            lexeme=token.getLexeme()
            raise Exception("Syntax error near \""+lexeme+"\"")

    def peek(self):
        if(not self.nextTok):
            file=self.file
            self.skipWhitespace()
            c=self.getChar()
            if(c is None): return None
            tokenType=None
            lexeme=""
            if(c=="(" or c==")"): tokenType=c
            else:
                tokenType="L"
                lexeme=c
                while(True):
                    c=self.getChar()
                    if(not c): break
                    if(re.search("[\s\(\)]",c)): break
                    lexeme+=c
                self.unGetChar(c)
            lexeme=lexeme.replace("&lparen;","(")
            lexeme=lexeme.replace("&rparen;",")")
            lexeme=lexeme.replace("&tab;","\t")
            lexeme=lexeme.replace("&space;"," ")
            self.nextTok=EssexToken(tokenType,lexeme)
        return self.nextTok

    def getChar(self):
        c=None
        if(self.ungot):
            c=self.ungot
            self.ungot=None
        else:
            c=self.file.read(1)
            if(len(c)==0): c=None
        return c

    def unGetChar(self,c):
        self.ungot=c

    def skipWhitespace(self):
        while(True):
            c=self.getChar()
            if(c is None): break
            if(c=="#"):
                while(True):
                    c=self.getChar()
                    if(c is None or c=="\n"): break
                continue
            if(not re.search("\s",c)):
                self.unGetChar(c)
                return

