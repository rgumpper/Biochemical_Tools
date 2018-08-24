# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 09:21:44 2018

@author: Ryan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
#get the sequence for NNS virus
seq=input('Please input your sequence for analysis:')

#number of Nucleotides in each nucleocapsid important for analysis
num=eval(input('How many nucleotides does each nucleocapsid contain?'))
seq=seq.replace(' ', '').replace('\n', '').replace('\r','')
#calculate # of A,G,C,T in sequence
A=0
G=0
C=0
T=0
for let in seq:
    if let=='a':
        A=A+1
    if let=='g':
        G=G+1
    if let=='c':
        C=C+1
    if let=='t':
        T=T+1
#calculate overall percentage and normalize to between 0 and 1
total=A+G+C+T
perA=A/total
perG=G/total
perC=C/total
perT=T/total

baselst=[perA, perG, perC, perT]

normA=(perA-min(baselst))/(max(baselst)-min(baselst))
normG=(perG-min(baselst))/(max(baselst)-min(baselst))
normC=(perC-min(baselst))/(max(baselst)-min(baselst))
normT=(perT-min(baselst))/(max(baselst)-min(baselst))

#sliding window to represent nucleotides encapsidated by each N subunit
seqlst=[]
for numb in range(0, len(seq), num):
    window=seq[numb:numb+9]
    seqlst.append(window)
    
#loop over and calculate weighted average for each subunit
weightlst=[]
for chunk in seqlst:
    score=0
    for letter in chunk:
        if letter=='a':
            score=score+normA
        if letter=='g':
            score=score+normG
        if letter=='c':
            score=score+normC
        if letter=='t':
            score=score+normT
    score=score/len(chunk)
    weightlst.append(score)

df=pd.DataFrame(data=weightlst)

plt.figure()
plt.xlabel('N-subunit Number')
plt.ylabel('Weighted Average Score')
plt.title('Subunit Stability as Nucleotide Abundance')
plt.plot(df.index, df[0], c='grey', linewidth=0.25)
plt.scatter(x=df.index, y=df[0], c=df[0], cmap='cool')

#for VSV to separate into codons data was taken from NCBI
#for N protein
N=[]
for i in range(63, 1332, 3):
    codon=seq[i:i+3]
    N.append(codon)
P=[]
for j in range(1395, 2193, 3):
    codon=seq[j:j+3]
    P.append(codon)
M=[]
for k in range(2249, 2939, 3):
    codon=seq[k:k+3]
    M.append(codon)
G=[]
for l in range(3077, 4613, 3):
    codon=seq[l:l+3]
    G.append(codon)
L=[]
for m in range(4732, 11061, 3):
    codon=seq[m:m+3]
    L.append(codon)
#codon dictionary for counting the codons in the different proteins  
codondict={("ttt", 'F'):0, ("ttc", 'F'):0, ("tta", 'L'):0, ("ttg", 'L'):0,
    ("tct", 'S'):0, ("tcc", 'S'):0, ("tca", 'S'):0, ("tcg", 'S'):0,
    ("tat", 'Y'):0, ("tac", 'Y'):0, ("taa", 'STOP'):0, ("tag", 'STOP'):0,
    ("tgt", 'C'):0, ("tgc", 'C'):0, ("tga", 'STOP'):0, ("tgg", 'W'):0,
    ("ctt", 'L'):0, ("ctc", 'L'):0, ("cta", 'L'):0, ("ctg", 'L'):0,
    ("cct", 'P'):0, ("ccc", 'P'):0, ("cca", 'P'):0, ("ccg", 'P'):0,
    ("cat", 'H'):0, ("cac", 'H'):0, ("caa", 'Q'):0, ("cag", 'Q'):0,
    ("cgt", 'R'):0, ("cgc", 'R'):0, ("cga", 'R'):0, ("cgg", 'R'):0,
    ("att", 'I'):0, ("atc", 'I'):0, ("ata", 'I'):0, ("atg", 'M'):0,
    ("act", 'T'):0, ("acc", 'T'):0, ("aca", 'T'):0, ("acg", 'T'):0,
    ("aat", 'N'):0, ("aac", 'N'):0, ("aaa", 'K'):0, ("aag", 'K'):0,
    ("agt", 'S'):0, ("agc", 'S'):0, ("aga", 'R'):0, ("agg", 'R'):0,
    ("gtt", 'V'):0, ("gtc", 'V'):0, ("gta", 'V'):0, ("gtg", 'V'):0,
    ("gct", 'A'):0, ("gcc", 'A'):0, ("gca", 'A'):0, ("gcg", 'A'):0,
    ("gat", 'D'):0, ("gac", 'D'):0, ("gaa", 'E'):0, ("gag", 'E'):0,
    ("ggt", 'G'):0, ("ggc", 'G'):0, ("gga", 'G'):0, ("ggg", 'G'):0}

#counting the actual codons and placing values in the dictionary
proteome=(N, P, M, G, L)

for prot in proteome:
    for cod in prot:
        for i in codondict:
            if cod in i:
                codondict[i]+=1
dictionaryvalues=list(codondict.values())
codonstdev=np.std(dictionaryvalues)
codonaverage=np.average(dictionaryvalues)

stdvalue=[]

for j in dictionaryvalues:
    x=(j-codonaverage)/codonstdev
    stdvalue.append(x)

keys=list(codondict.keys())
#creation of dataframes from dictionary for easier data manipulation and stadardization
dictdf=pd.DataFrame.from_dict(data=codondict, orient='index')
standardizeddf=(dictdf-dictdf.mean())/dictdf.std()
usagedf=pd.DataFrame(data=totallst, index=keys, dtype=float)
xrange=[]
standardizeddf.columns=['A']
for i in range(1, 65):
    xrange.append(i)

plt.figure()
plt.xlabel('Codon')
plt.ylabel('Standardized Codon Usage')
plt.title('Standardized Codon Usage Across VSV')
plt.scatter(x=xrange, y=standardizeddf['A'], c=standardizeddf['A'], cmap='cool')
plt.xticks(xrange, standardizeddf.index, rotation=90)    

newindex=[]
for index,row in standardizeddf.itertuples():
    newindex.append(index[0])
finalstanddf=standardizeddf
finalstanddf.index=newindex
normalizeddf=((dictdf-dictdf.min())/(dictdf.max()-dictdf.min()))
normalizeddf.index=newindex
normalizeddf.columns=['A']

count=0
numcodons=[]
score=[]
while count <= 51:
    count=count+1
    numcodons.append(count)
    finaladdition=[]
    for g in range(0, 51, count):
        runscore=[]
        b=N[g:g+count]
        for h in b:
            m=normalizeddf.loc[h, 'A']
            runscore.append(m)
        appendscore=np.average(runscore)
        finaladdition.append(appendscore)
    newfinalscore=finaladdition
    finalscore=np.average(finaladdition)
    score.append(finalscore)
                
plt.figure
plt.xlabel('Number of Codons')
plt.ylabel('Average Standardized Codon Score')
plt.title('Average Codon Score as a Funciton of Number of Codons')
plt.scatter(x=numcodons, y=score, c=score, cmap='cool')       
        


  
changeCodondict={"uuu":"F", "uuc":"F", "uua":"L", "uug":"L",
    "ucu":"S", "ucc":"S", "uca":"S", "ucg":"S",
    "uau":"Y", "uac":"Y", "uaa":"STOP", "uag":"STOP",
    "ugu":"C", "ugc":"C", "uga":"STOP", "ugg":"W",
    "cuu":"L", "cuc":"L", "cua":"L", "cug":"L",
    "ccu":"P", "ccc":"P", "cca":"P", "ccg":"P",
    "cau":"H", "cac":"H", "caa":"Q", "cag":"Q",
    "cgu":"R", "cgc":"R", "cga":"R", "cgg":"R",
    "auu":"I", "auc":"I", "aua":"I", "aug":"M",
    "acu":"T", "acc":"T", "aca":"T", "acg":"T",
    "aau":"N", "aac":"N", "aaa":"K", "aag":"K",
    "agu":"S", "agc":"S", "aga":"R", "agg":"R",
    "guu":"V", "guc":"V", "gua":"V", "gug":"V",
    "gcu":"A", "gcc":"A", "gca":"A", "gcg":"A",
    "gau":"D", "gac":"D", "gaa":"E", "gag":"E",
    "ggu":"G", "ggc":"G", "gga":"G", "ggg":"G"}
standardcodondict={("ttt", 'F'):0, ("ttc", 'F'):0, ("tta", 'L'):0, ("ttg", 'L'):0,
    ("tct", 'S'):0, ("tcc", 'S'):0, ("tca", 'S'):0, ("tcg", 'S'):0,
    ("tat", 'Y'):0, ("tac", 'Y'):0, ("taa", 'STOP'):0, ("tag", 'STOP'):0,
    ("tgt", 'C'):0, ("tgc", 'C'):0, ("tga", 'STOP'):0, ("tgg", 'W'):0,
    ("ctt", 'L'):0, ("ctc", 'L'):0, ("cta", 'L'):0, ("ctg", 'L'):0,
    ("cct", 'P'):0, ("ccc", 'P'):0, ("cca", 'P'):0, ("ccg", 'P'):0,
    ("cat", 'H'):0, ("cac", 'H'):0, ("caa", 'Q'):0, ("cag", 'Q'):0,
    ("cgt", 'R'):0, ("cgc", 'R'):0, ("cga", 'R'):0, ("cgg", 'R'):0,
    ("att", 'I'):0, ("atc", 'I'):0, ("ata", 'I'):0, ("atg", 'M'):0,
    ("act", 'T'):0, ("acc", 'T'):0, ("aca", 'T'):0, ("acg", 'T'):0,
    ("aat", 'N'):0, ("aac", 'N'):0, ("aaa", 'K'):0, ("aag", 'K'):0,
    ("agt", 'S'):0, ("agc", 'S'):0, ("aga", 'R'):0, ("agg", 'R'):0,
    ("gtt", 'V'):0, ("gtc", 'V'):0, ("gta", 'V'):0, ("gtg", 'V'):0,
    ("gct", 'A'):0, ("gcc", 'A'):0, ("gca", 'A'):0, ("gcg", 'A'):0,
    ("gat", 'D'):0, ("gac", 'D'):0, ("gaa", 'E'):0, ("gag", 'E'):0,
    ("ggt", 'G'):0, ("ggc", 'G'):0, ("gga", 'G'):0, ("ggg", 'G'):0}

count1=1
while count1<=64:
    easylist=[]
    








    
    
    




