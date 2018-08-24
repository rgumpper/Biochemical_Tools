# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 12:00:16 2018

@author: Ryan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches



def get_b_factors_RNA(pdb_path, y = 1):
    pdb = open(pdb_path)
    array1 = []
    for i in pdb.readlines():
        if i[0:4] == 'ATOM':
            line = i.split()
            array1.append(line)
    df = pd.DataFrame(data = array1)
    if len(df[11].isnull()) > 0:
        df[11].fillna(value = df[10], inplace = True)
        df[10].replace(['N', 'C', 'O', 'P'], value = np.nan, inplace = True)
        df[9], df[12] = df[9].str.split('.00', 1).str
        df[10].fillna(value = df[12], inplace = True)
    df=df.iloc[:,[2,3,4,5,10]]
    df=df.rename(columns={2:'Atom', 3:'Residue', 4:'Chain', 5:'Noresidue', 
                        10:'Bfactor'})
    df['Noresidue'] = pd.to_numeric(df['Noresidue'], errors = 'coerce', downcast = 'float')
    df['Bfactor'] = pd.to_numeric(df['Bfactor'], errors = 'coerce', downcast = 'float') 
    df2 = pd.pivot_table(df, index = 'Noresidue', columns = 'Chain', 
                        values = 'Bfactor', aggfunc = np.mean)
    average = df['Bfactor'].mean()
    df['Bfactor'] = (df['Bfactor']-df['Bfactor'].mean())/df['Bfactor'].std()
    if y == 1:
        df = df[(df['Residue']=='U') | (df['Residue']=='G') | (df['Residue']=='A') | (df['Residue']=='C')]
    df = pd.pivot_table(df, index = 'Noresidue', columns = 'Chain', 
                        values = 'Bfactor', aggfunc = np.mean)
    pdb.close()
    return df, average, df2

random, r, rawR = get_b_factors_RNA('2gic_original.pdb', y = 0)
A, a, rawA = get_b_factors_RNA('3ptx.pdb', y = 0)
G, g, rawG = get_b_factors_RNA('3pu1.pdb',  y = 0)
C, c, rawC = get_b_factors_RNA('3pu0.pdb', y = 0)
U, u, rawU = get_b_factors_RNA('3pu4.pdb', y = 0)


for chain in ['A', 'B', 'C', 'D', 'E', 'R']:
    plt.figure()
    plt.scatter(x = A.index, y = A[chain], color = 'red', alpha = 0.5)
    plt.scatter(x = G.index, y = G[chain], color = 'blue', alpha = 0.5)
    plt.scatter(x = C.index, y = C[chain], color = 'green', alpha = 0.5)
    plt.scatter(x = U.index, y = U[chain], color = 'purple', alpha = 0.5)
    plt.scatter(x = random.index, y = random[chain], color = 'gray', alpha = 0.5)
    red_patch = mpatches.Patch(color = 'red', label ='A', alpha = 0.5)
    blue_patch = mpatches.Patch(color = 'blue', label = 'G', alpha = 0.5)
    green_patch = mpatches.Patch(color = 'green', label = 'C', alpha = 0.5)
    purple_patch = mpatches.Patch(color = 'purple', label = 'U', alpha = 0.5)
    random_patch = mpatches.Patch(color = 'gray', label = 'Random', alpha = 0.5)
    plt.legend(handles = [red_patch, blue_patch, green_patch, 
                          purple_patch, random_patch], loc = 'upper right', ncol = 2)
    plt.plot(A.index,  A[chain], linewidth = 0.5, color = 'red', alpha = 0.5)
    plt.plot(G.index,  G[chain], linewidth = 0.5, color = 'blue', alpha = 0.5)
    plt.plot(C.index,  C[chain], linewidth = 0.5, color = 'green', alpha = 0.5)
    plt.plot(U.index,  U[chain], linewidth = 0.5, color = 'purple', alpha = 0.5)
    plt.plot(random.index,  random[chain], linewidth = 0.5, color = 'gray', alpha = 0.5)
    plt.xlabel('Residue')
    plt.ylabel('Standardized Average B-factor')
    plt.ylim((-3, 8))
    plt.title('Chain '+ chain)
    plt.show()
    
    
for chain in ['A', 'B', 'C', 'D', 'E', 'R']:
    plt.figure()
    plt.scatter(x = rawA.index, y = rawA[chain], color = 'red', alpha = 0.5)
    plt.scatter(x = rawG.index, y = rawG[chain], color = 'blue', alpha = 0.5)
    plt.scatter(x = rawC.index, y = rawC[chain], color = 'green', alpha = 0.5)
    plt.scatter(x = rawU.index, y = rawU[chain], color = 'purple', alpha = 0.5)
    plt.scatter(x = rawR.index, y = rawR[chain], color = 'gray', alpha = 0.5)
    red_patch = mpatches.Patch(color = 'red', label ='A', alpha = 0.5)
    blue_patch = mpatches.Patch(color = 'blue', label = 'G', alpha = 0.5)
    green_patch = mpatches.Patch(color = 'green', label = 'C', alpha = 0.5)
    purple_patch = mpatches.Patch(color = 'purple', label = 'U', alpha = 0.5)
    random_patch = mpatches.Patch(color = 'gray', label = 'Random', alpha = 0.5)
    plt.legend(handles = [red_patch, blue_patch, green_patch, 
                          purple_patch, random_patch], loc = 'upper right', ncol = 2)
    plt.plot(rawA.index,  rawA[chain], linewidth = 0.5, color = 'red', alpha = 0.5)
    plt.plot(rawG.index,  rawG[chain], linewidth = 0.5, color = 'blue', alpha = 0.5)
    plt.plot(rawC.index,  rawC[chain], linewidth = 0.5, color = 'green', alpha = 0.5)
    plt.plot(rawU.index,  rawU[chain], linewidth = 0.5, color = 'purple', alpha = 0.5)
    plt.plot(rawR.index,  rawR[chain], linewidth = 0.5, color = 'gray', alpha = 0.5)
    plt.xlabel('Residue')
    plt.ylabel('Raw B-factor')
    plt.title('Chain '+ chain)
    plt.show()


