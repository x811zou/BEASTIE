#====================================================================
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURxCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
#====================================================================
from .NgramIterator import NgramIterator

######################################################################
# Attributes:
#    codon (class attribute) hash: codon -> amino acid
#    complementMap (class attribute): hash: nuc -> revcomp(nuc)
# Class Methods:
#    aaSeq=Translation.translate(nucSeq)
#    revSeq=Translation.reverseComplement(nucSeq)
#    hash=Translation.getFourfoldDegenerateCodons()
# Private methods:
#    initCodonMap()
#    initComplementMap()
#    nucleotide=complement($nucleotide)
######################################################################

class Translation:
    """Translation provides routines for translating or reverse-
    complementing nucleotide sequences
    """

    codon={}
    complementMap={}

    @classmethod
    def translate(cls,transcript):
        translation=""
        L=len(transcript)
        for i in range(0,L-3,3):
            three=transcript[i:i+3]
            residue=cls.codon[three]
            if(not residue): residue="X"
            translation+=residue
        return translation

    @classmethod
    def reverseComplement(cls,seq):
        L=len(seq)
        buffer=""
        for i in range(0,L):
            buffer+=cls.complement(seq[L-i-1:L-i])
        return buffer

    @classmethod
    def complement(cls,c):
        if(cls.complementMap.get(c,None) is None): cls.complementMap[c]="N";
        return cls.complementMap[c]

    @classmethod
    def initComplementMap(cls):
        cls.complementMap['A']='T'
        cls.complementMap['T']='A'
        cls.complementMap['G']='C'
        cls.complementMap['C']='G'
        cls.complementMap['R']='Y'
        cls.complementMap['Y']='R'
        cls.complementMap['-']='-' # for gaps in alignments
        cls.complementMap['N']='N' ### <------DEBUGGING!

    @classmethod
    def getFourfoldDegenerateCodons(cls):
        codonMap=cls.codon
        alphabet="ACGT"
        degenerate=set()
        ngramIterator=NgramIterator(alphabet,2)
        while(True):
            pair=ngramIterator.nextString()
            if(pair is None): break
            acids=set(); codons=set()
            for third in alphabet:
                codon=pair+third
                acid=codonMap[codon]
                acids.add(acid)
                codons.add(codon)
            if(len(acids)>1): continue
            for codon in codons: degenerate.add(codon)
        return degenerate


    @classmethod
    def initCodonMap(cls):
        cls.codon["GTT"]='V'
        cls.codon["GTC"]='V'
        cls.codon["GTA"]='V'
        cls.codon["GTG"]='V'
        cls.codon["GTR"]='V'
        cls.codon["GTY"]='V'
        cls.codon["GTN"]='V'
        cls.codon["GCT"]='A'
        cls.codon["GCC"]='A'
        cls.codon["GCA"]='A'
        cls.codon["GCG"]='A'
        cls.codon["GCR"]='A'
        cls.codon["GCY"]='A'
        cls.codon["GCN"]='A'
        cls.codon["GAT"]='D'
        cls.codon["GAC"]='D'
        cls.codon["GAY"]='D'
        cls.codon["GAA"]='E'
        cls.codon["GAG"]='E'
        cls.codon["GAR"]='E'
        cls.codon["GGT"]='G'
        cls.codon["GGC"]='G'
        cls.codon["GGA"]='G'
        cls.codon["GGG"]='G'
        cls.codon["GGR"]='G'
        cls.codon["GGY"]='G'
        cls.codon["GGN"]='G'
        cls.codon["TTT"]='F'
        cls.codon["TTC"]='F'
        cls.codon["TTY"]='F'
        cls.codon["TTA"]='L'
        cls.codon["TTG"]='L'
        cls.codon["TTR"]='L'
        cls.codon["CTT"]='L'
        cls.codon["CTC"]='L'
        cls.codon["CTA"]='L'
        cls.codon["CTG"]='L'
        cls.codon["CTN"]='L'
        cls.codon["YTR"]='L'
        cls.codon["TCT"]='S'
        cls.codon["TCC"]='S'
        cls.codon["TCA"]='S'
        cls.codon["TCG"]='S'
        cls.codon["TCN"]='S'
        cls.codon["AGT"]='S'
        cls.codon["AGC"]='S'
        cls.codon["AGY"]='S'
        cls.codon["TAT"]='Y'
        cls.codon["TAC"]='Y'
        cls.codon["TAY"]='Y'
        cls.codon["TAA"]='*'
        cls.codon["TAG"]='*'
        cls.codon["TAR"]='*'
        cls.codon["TAG"]='*'
        cls.codon["TGT"]='C'
        cls.codon["TGC"]='C'
        cls.codon["TGY"]='C'
        cls.codon["TGA"]='*'
        cls.codon["TGG"]='W'
        cls.codon["CCT"]='P'
        cls.codon["CCC"]='P'
        cls.codon["CCA"]='P'
        cls.codon["CCG"]='P'
        cls.codon["CCR"]='P'
        cls.codon["CCY"]='P'
        cls.codon["CCN"]='P'
        cls.codon["CAT"]='H'
        cls.codon["CAC"]='H'
        cls.codon["CAY"]='H'
        cls.codon["CAA"]='Q'
        cls.codon["CAG"]='Q'
        cls.codon["CAR"]='Q'
        cls.codon["CGT"]='R'
        cls.codon["CGC"]='R'
        cls.codon["CGA"]='R'
        cls.codon["CGG"]='R'
        cls.codon["CGN"]='R'
        cls.codon["ATT"]='I'
        cls.codon["ATC"]='I'
        cls.codon["ATA"]='I'
        cls.codon["ATH"]='I'
        cls.codon["ATG"]='M'
        cls.codon["ACT"]='T'
        cls.codon["ACC"]='T'
        cls.codon["ACA"]='T'
        cls.codon["ACG"]='T'
        cls.codon["ACN"]='T'
        cls.codon["AAT"]='N'
        cls.codon["AAC"]='N'
        cls.codon["AAY"]='N'
        cls.codon["AAA"]='K'
        cls.codon["AAG"]='K'
        cls.codon["AAR"]='K'
        cls.codon["AGG"]='R'
        cls.codon["AGA"]='R'
        cls.codon["AGY"]='R'

Translation.initCodonMap()
Translation.initComplementMap();


