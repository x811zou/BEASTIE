#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
import copy
import re

######################################################################
#
# A node in a parse tree for an Essex file (based on S-expressions).
# These parse trees are produced by the EssexParser.  Literals (such as
# numbers and strings---i.e., the actual data) are not allocated an
# EssexNode; an EssexNode represents a parenthesized expression only.
# You can test whether something is an EssexNode via isaNode(), below.
#
# Attributes:
#   tag : string
#   elements : array
# Methods:
#   node=EssexNode([tag,elem1,elem2,...])
#   node->addElem(elem)
#   tag=node.getTag()
#   node.changeTag(newTag)
#   n=node.numElements()
#   elem=node.getIthElem(i) # 0-based, doesn't include the tag
#   elem=node[i] # 0-based, doesn't include the tag
#   node.setIthElem(i,dataOrNode)
#   elem=node.findChild(tag)
#   array=node.findChildren(tag)
#   node.dropChild(i)
#   array=node.findDescendents(tag) # always returns an array
#   elem=node.findDescendent(tag) # returns node or undef
#   bool=node.hasDescendentOrDatum(tagOrDatum)
#   count=node.countDescendentOrDatum(tagOrDatum)
#   n=node.countDescendents(tag)
#   bool=node.hasDescendent(tag)
#   string=node.getAttribute(attributeTag)
#   node.setAttribute(tag,value)
#   array=node.getElements()
#   bool=EssexNode.isaNode(datum)
#   bool=node.hasCompositeChildren()
#   node.print(filehandle)
#   node.printXML(filehandle)
#   node.recurse(visitor) # must have methods enter(node) and leave(node)
#   oneElem=node.pathQuery("report/reference-transcript/type")
######################################################################

class EssexNode:
    EXONS={"initial-exon":1, "internal-exon":1, "final-exon":1,
           "single-exon":1, "exon":1, "five-prime-UTR":1,
           "three-prime-UTR":1, "five_prime_UTR":1, "three_prime_UTR":1,
           "five-prime-UTR":1, "three-prime-UTR":1}
    SINGLETONS={"bad-annotation":1, "bad-start":1, "no-stop-codon":1,
                "NMD":1, "mapped":1, "noncoding":1, "protein-differs":1,
                "too-many-vcf-errors":1, "hypothetical-NMD":1}
    VARIANTS={"CDS-variants":1, "frameshift-variants":1, "UTR-variants":1,
              "near-splice-variants":1, "splice-site-variants":1}

    def __init__(self,parms):
        if(len(parms)>0):
            self.tag=parms.pop(0)
            if(len(parms)>0):
                self.elements=copy.deepcopy(parms)
            else: self.elements=[]
        else:
            self.tag=""
            self.elements=[]

    def dropChild(self,i):
        del self.elements[i]

    def addElem(self,elem):
        self.elements.append(elem)

    def setAttribute(self,tag,value):
        elements=self.elements
        for elem in elements:
            if(EssexNode.isaNode(elem) and elem.getTag()==tag):
                elem.setIthElem(0,value)
                return
        elements.append(EssexNode([tag,value]))

    def getTag(self):
        return self.tag

    def changeTag(self,tag):
        self.tag=tag

    def numElements(self):
        elements=self.elements
        return len(elements) if elements else 0

    def __getitem__(self,i):
        return self.elements[i]

    def getIthElem(self,i):
        return self.elements[i]

    def setIthElem(self,i,child):
        self.elements[i]=child

    def getElements(self):
        return self.elements

    @classmethod
    def isaNode(cls,x):
        return isinstance(x,EssexNode)

    def findChild(self,tag):
        elements=self.elements
        for elem in elements:
            if(EssexNode.isaNode(elem) and elem.getTag()==tag):
                return elem
        return None

    def findChildren(self,tag):
        results=[]
        elements=self.elements
        for elem in elements:
            if(EssexNode.isaNode(elem) and elem.getTag()==tag):
                results.append(elem)
        return results

    def getAttribute(self,tag):
        elements=self.elements
        for elem in elements:
            if(EssexNode.isaNode(elem) and elem.getTag()==tag):
                return elem.getIthElem(0)
        return None

    def print(self,file):
        self.printRecursive(0,file)

    def printXML(self,file):
        self.printRecursiveXML(0,file)

    def hasDescendent(self,tag):
        children=self.elements
        if(children):
            for child in children:
                if(EssexNode.isaNode(child)):
                    if(child.tag==tag or child.hasDescendent(tag)): return True
                elif(child==tag): return True
        return False

    def countDescendents(self,tag):
        array=self.findDescendents(tag)
        return len(array)

    def findDescendents(self,tag):
        array=[]
        self.findDesc(tag,array)
        return array

    def pathQuery(self,query):
        fields=query.split("/")
        if(len(fields)<2): exit("Essex query must contain at least one '/'\n")
        rootTag=fields.pop(0)
        depth=len(fields)
        candidates=[]
        if(self.getTag()==rootTag): candidates=[self]
        else: candidates=self.findDescendents(rootTag)
        n=len(candidates)
        for i in range(n):
            candidate=candidates[i]
            attr=candidate
            for j in range(depth):
                attr=attr.findChild(fields[j])
                if(attr is None): break
            if(attr is not None): return attr
        return None

    def recurse(self,visitor):
        visitor.enter(self)
        elements=self.elements
        for elem in elements:
            if(EssexNode.isaNode(elem)): elem.recurse(visitor)
            else:
                visitor.enter(elem)
                visitor.leave(elem)
        visitor.leave(self)

    def findDesc(self,tag,array):
        children=self.elements
        for child in children:
            if(EssexNode.isaNode(child)):
                if(child.tag==tag): array.append(child)
                else: child.findDesc(tag,array)

    def findDescendent(self,tag):
        children=self.elements
        for child in children:
            if(EssexNode.isaNode(child) and child.tag==tag): return child
        for child in children:
            if(EssexNode.isaNode(child)):
                desc=child.findDescendent(tag)
                if(desc): return desc
        return None

    def countDescendentOrDatum(self,tag):
        count=0
        if(self.tag==tag): count+=1
        children=self.elements
        for child in children:
            if(EssexNode.isaNode(child)):
                count+=child.countDescendentOrDatum(tag)
            else:  # not a node
                if(child==tag): count+=1
        return count

    def hasDescendentOrDatum(self,tag):
        if(self.tag==tag): return True
        children=self.elements
        for child in children:
            if(EssexNode.isaNode(child)):
                if(child.hasDescendentOrDatum(tag)): return True
            else: # not a node
                if(child==tag): return True
        return False

    def hasCompositeChildren(self):
        children=self.elements
        numChildren=len(children) if children else 0
        for i in range(numChildren):
            if(EssexNode.isaNode(children[i])): return True
        return False

    def printRecursive(self,depth,file):
        tab='   '*depth
        tag=self.tag
        file.write(tab+"("+tag);
        elements=self.elements
        n=len(elements) if elements else 0
        if(self.hasCompositeChildren()):
            for i in range(n):
                child=elements[i]
                if(EssexNode.isaNode(child)):
                    file.write("\n")
                    child.printRecursive(depth+1,file)
                else:
                    tab='   '*(depth+1)
                    file.write("\n"+tab+child)
        else:
             for i in range(n):
                 elem=elements[i]
                 file.write(" "+str(elem))
        file.write(")")

    def printExonXML(self,tag,tab,depth,file):
        begin=self.getIthElem(0)
        end=self.getIthElem(1)
        score=self.getIthElem(2)
        strand=self.getIthElem(3)
        frame=self.getIthElem(4)
        file.write(tab+"<"+tag+" begin=\""+begin+"\" end=\""+end+"\" score=\""
                   +score+"\" strand=\""+strand+"\" frame=\""+frame+"\"/>")

    def printRecursiveXML(self,depth,file):
        tab='   '*depth
        tag=self.tag
        tag=tag.replace("_","-")
        elements=self.elements
        n=len(elements) if elements else 0
        if(n==0):
            file.write(tab+"<"+tag+"/>")
            return
        if(EssexNode.EXONS.get(tag,None)):
            self.printExonXML(tag,tab,depth,file)
            return
        if(tag=="bad-donor" or tag=="bad-acceptor"):
            nuc=self.getIthElem(0)
            pos=self.getIthElem(1)
            file.write(tab+"<"+tag+" consensus=\""+nuc+"\" pos=\""+pos+"\"/>")
            return
        if(tag=="new-start" or tag=="old-start"):
            seq=self.getIthElem(0)
            score=self.getIthElem(1)
            threshold=self.getIthElem(3)
            file.write(tab+"<"+tag+" seq=\""+seq+"\" score=\""+score+
                       "\" threshold=\""+threshold+"\"/>")
            return
        if(tag=="percent-match"):
            percent=self.getIthElem(0)
            ratio=self.getIthElem(1)
            refLen=self.getIthElem(2)
            altLen=self.getIthElem(3)
            match=re.search("ref-length=(\d+)",refLen)
            if(not match): exit("expected: ref-length=\d+")
            refLen=match.group(1)
            match=re.search("alt-length=(\d+)",altLen)
            if(not match): exit("expected: alt-length=\d+")
            altLen=match.group(1)
            file.write(tab+"<"+tag+" percent=\""+percent+"\" ratio=\""+ratio+
                       "\" ref-length=\""+refLen+"\" alt-length=\""+altLen+
                       "\"/>")
            return
        if(tag=="ORF-length"):
            oldLen=self.getIthElem(0)
            newLen=self.getIthElem(2)
            file.write(tab+"<"+tag+" old=\""+oldLen+"\" new=\""+newLen+"\"/>")
            return
        if(tag=="broken-acceptor" or tag=="broken-donor"):
            pos=self.getIthElem(0)
            oldSeq=self.getIthElem(1)
            oldScore=self.getIthElem(2)
            newSeq=self.getIthElem(3)
            newScore=self.getIthElem(4)
            threshold=self.getIthElem(6)
            file.write(tab+"<"+tag+" pos=\""+pos+"\" ref-seq=\""+oldSeq+
                       "\" ref-score=\""+oldScore+"\" alt-seq=\""+newSeq+
                       "\" alt-score=\""+newScore+"\" threshold=\""+
                       threshold+"\"/>")
            return
        if(tag=="cryptic-site"):
            type=self.getIthElem(0)
            pos=self.getIthElem(1)
            seq=self.getIthElem(2)
            score=self.getIthElem(3)
            threshold=self.getIthElem(5)
            file.write(tab+"<"+tag+" type=\""+type+"\" pos=\""+pos+"\" seq=\""
                       +seq+"\" score=\""+score+"\" threshold=\""+threshold+
                       "\"/>")
            return
        if(EssexNode.VARIANTS.get(tag,None)):
            file.write(tab+"<"+tag+">\n")
            for i in range(n):
                variant=self.getIthElem(i)
                fields=variant.split(":")
                (id,chr,refPos,altPos,refSeq,altSeq)=fields
                file.write(tab+"   <variant id=\""+id+"\" chr=\""+chr+
                           "\" ref-pos=\""+refPos+"\" alt-pos=\""+altPos+
                           "\" ref-seq=\""+refSeq+"\" alt-seq=\""+altSeq+
                           "\"/>\n")
            file.write(tab+"</"+tag+">")
            return
        if(tag=="global-coords"):
            fields=self.getIthElem(0).split(":")
            (chr,coords,strand)=fields
            match=re.search("(\d+)-(\d+)",coords)
            if(not match): raise Exception(coords)
            (begin,end)=(match.group(1),match.group(2))
            file.write(tab+"<"+tag+" chr=\""+chr+"\" begin=\""+begin+
                       "\" end=\""+end+"\" strand=\""+strand+"\"/>")
            return
        file.write(tab+"<"+tag+">")
        if(self.hasCompositeChildren()):
            for i in range(n):
                child=elements[i]
                if(EssexNode.isaNode(child)):
                    file.write("\n")
                    child.printRecursiveXML(depth+1,file)
                else:
                    tab='   '*(depth+1)
                    if(EssexNode.SINGLETONS.get(child,None)):
                        file.write("\n"+tab+"<"+child+"/>")
                    else: file.write("\n"+tab+child)
            file.write("\n"+tab+"</"+tag+">")
        else:
            for i in range(n):
                elem=elements[i]
                if(EssexNode.SINGLETONS.get(elem,None)):
                    file.write("<"+elem+"/>")
                else:
                    if(i>0): file.write(" ")
                    file.write(elem)
            file.write("</"+tag+">")
