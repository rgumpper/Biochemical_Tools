#Analysis of N sequence


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
#open the file of the sequence
textseq=open('N_seq.txt')
#dictionary for codons
codons={"UUU" : 0, "CUU" : 0, "AUU" : 0, "GUU" : 0, "UUC" : 0, "CUC" : 0, "AUC" : 0, "GUC" : 0, 
        "UUA" : 0, "CUA" : 0, "AUA" : 0, "GUA" : 0, "UUG" : 0, "CUG" : 0, "AUG" : 0, "GUG" : 0, 
        "UCU" : 0, "CCU" : 0, "ACU" : 0, "GCU" : 0, "UCC" : 0, "CCC" : 0, "ACC" : 0, "GCC" : 0,
        "UCA" : 0, "CCA" : 0, "ACA" : 0, "GCA" : 0, "UCG" : 0, "CCG" : 0, "ACG" : 0, "GCG" : 0, 
        "UAU" : 0, "CAU" : 0, "AAU" : 0, "GAU" : 0, "UAC" : 0, "CAC" : 0, "AAC" : 0, "GAC" : 0, 
        "UAA" : 0, "CAA" : 0, "AAA" : 0, "GAA" : 0, "UAG" : 0, "CAG" : 0, "AAG" : 0, "GAG" : 0, 
        "UGU" : 0, "CGU" : 0, "AGU" : 0, "GGU" : 0, "UGC" : 0, "CGC" : 0, "AGC" : 0, "GGC" : 0, 
        "UGA" : 0, "CGA" : 0, "AGA" : 0, "GGA" : 0, "UGG" : 0, "CGG" : 0, "AGG" : 0, "GGG" : 0}
#reading the file and changing out the T's to U's
seq=textseq.read()
seq=seq.upper()
seqlist=[]
for letter in seq:
    if letter=='T':
        seqlist.append('U')
    else:
        seqlist.append(letter)
#final sequence to work on
seqtransformed=''.join(seqlist)
codonlist=[]
for triplet in range(0, len(seqtransformed)+1, 3):
    a=seqtransformed[triplet:triplet+3]
    codonlist.append(a)

d=[]
for triplet in codonlist:
    if triplet in codons:
        d.append([triplet, codons[triplet]+=1])


df=pd.DataFrame.from_dict(data=codons, orient='index', dtype=None)

df2=pd.DataFrame()

x=codonlist
y=["UUU", "CUU", "AUU", "GUU", "UUC", "CUC", "AUC", "GUC", 
   "UUA", "CUA", "AUA", "GUA", "UUG", "CUG", "AUG", "GUG",
   "UCU", "CCU", "ACU", "GCU", "UCC", "CCC", "ACC", "GCC", 
   "UCA", "CCA", "ACA", "GCA", "UCG", "CCG", "ACG", "GCG", 
   "UAU", "CAU", "AAU", "GAU", "UAC", "CAC", "AAC", "GAC", 
   "UAA", "CAA", "AAA", "GAA", "UAG", "CAG", "AAG", "GAG", 
   "UGU", "CGU", "AGU", "GGU", "UGC", "CGC", "AGC", "GGC", 
   "UGA", "CGA", "AGA", "GGA", "UGG", "CGG", "AGG", "GGG"]
z=[]
for triplet in x:
    for trial in y:
        count=0
        if triplet==trial:
            count=count+1
            z.append(count)
    
    
    