# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 11:18:58 2018

@author: Ryan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches

#Program to plot (A+G)% and Kappa Index of Coincidence


def purine_percent(sequence, window):
    seq = sequence
    seq = seq.upper()
    seq = seq.replace(' ', '')
    seq = seq.replace('T', 'U')
    A = seq.count('A')
    G = seq.count('G')
    C = seq.count('C')
    U = seq.count('U')
    AGtot = (100/(A+G+U+C))*(A+G)
    
    AGsw = []
    
    for i in range(0, len(seq)-window):
        win_seq = seq[i:i+window]
        Asw = win_seq.count('A')
        Gsw = win_seq.count('G')
        Csw = win_seq.count('C')
        Usw = win_seq.count('U')
        
        score = (AGtot/(Asw + Gsw + Csw + Usw))*(Asw+Gsw)
        AGsw.append(score)
    return AGsw

def kappa_index(sequence, window):
    seq = sequence
    seq = seq.upper()
    seq = seq.replace(' ', '')
    seq = seq.replace('T', 'U')
    KIC = []
    for i in range(0, len(seq)-window):
        t=0
        win_seq = seq[i:i+window]
        for j in range(0, len(win_seq)):
            C = 0
            if j+1 < len(win_seq):
                b = win_seq[j+1:]
            elif j+1 == len(win_seq):
                b = win_seq[j]
            for u in b:
                if win_seq[j] == u:
                    C = C+1
            t = t+(C/len(b)*100)
        IC = t/(window-1)
        KIC.append(IC)
    return KIC
   
def linguistic_complexity(sequence, window, alphabet_length=4):
    seq = sequence 
    seq = seq.upper()
    seq = seq.replace(' ', '')
    seq = seq.replace('T', 'U')
    CLtot = 0
    for i in range(1, window+1):
        if alphabet_length**i <= window:
            CLtot = CLtot + alphabet_length**i
        else:
            CLtot = CLtot + (window-(i-1))
    CL = []
    for i in range(0, len(seq)-window):
        CLwin = 0
        win_seq = seq[i:i+window]
        CLwin_seq = []
        for i in range(1, window+1):
            for j in range(0, window+1, i):
                if win_seq[j:j+i] not in CLwin_seq:
                    CLwin_seq.append(win_seq[j:j+i])
                    CLwin = CLwin + 1
                else:
                    pass
        score = CLwin/CLtot
        CL.append(score)
    return CL

normal_pur = purine_percent('atgtctgttacagtcaagagaatcattgacaacacagtcatagttccaaaacttcctgcaaatgaggatccagtggaatacccggcagattacttcagaaaatcaaaggagattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaacagctacttgtatggagcattaaaggacatccggggtaagttggataaagattggtcaagtttcggaataaacatcgggaaagcaggggatacaatcggaatatttgaccttgtatccttgaaagccctggacggcgtacttccagatggagtatcggatgcttccagaaccagcgcagatgacaaatggttgcctttgtatctacttggcttatacagagtgggcagaacacaaatgcctgaatacagaaaaaagctcatggatgggctgacaaatcaatgcaaaatgatcaatgaacagtttgaacctcttgtgccagaaggtcgtgacatttttgatgtgtggggaaatgacagtaattacacaaaaattgtcgctgcagtggacatgttcttccacatgttcaaaaaacatgaatgtgcctcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgcaaaataaccggaatgtctacagaagatgtaacgacctggatcttgaaccgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgacaaggccgattcatacatgccttatttgatcgactttggattgtcttctaagtctccatattcttccgtcaaaaaccctgccttccacttctgggggcaattgacagctcttctgctcagatccaccagagcaaggaatgcccgacagcctgatgacattgagtatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgccgacttggcacaacagttttgtgttggagataacaaatacactccagatgatagtaccggaggattgacgactaatgcaccgccacaaggcagagatgtggtcgaatggctcggatggtttgaagatcaaaacagaaaaccgactcctgatatgatgcagtatgcgaaaagagcagtcatgtcactgcaaggcctaagagagaagacaattggcaagtatgctaagtcagaatttgacaaatga', window = 30)
normal_kap = kappa_index('atgtctgttacagtcaagagaatcattgacaacacagtcatagttccaaaacttcctgcaaatgaggatccagtggaatacccggcagattacttcagaaaatcaaaggagattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaacagctacttgtatggagcattaaaggacatccggggtaagttggataaagattggtcaagtttcggaataaacatcgggaaagcaggggatacaatcggaatatttgaccttgtatccttgaaagccctggacggcgtacttccagatggagtatcggatgcttccagaaccagcgcagatgacaaatggttgcctttgtatctacttggcttatacagagtgggcagaacacaaatgcctgaatacagaaaaaagctcatggatgggctgacaaatcaatgcaaaatgatcaatgaacagtttgaacctcttgtgccagaaggtcgtgacatttttgatgtgtggggaaatgacagtaattacacaaaaattgtcgctgcagtggacatgttcttccacatgttcaaaaaacatgaatgtgcctcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgcaaaataaccggaatgtctacagaagatgtaacgacctggatcttgaaccgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgacaaggccgattcatacatgccttatttgatcgactttggattgtcttctaagtctccatattcttccgtcaaaaaccctgccttccacttctgggggcaattgacagctcttctgctcagatccaccagagcaaggaatgcccgacagcctgatgacattgagtatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgccgacttggcacaacagttttgtgttggagataacaaatacactccagatgatagtaccggaggattgacgactaatgcaccgccacaaggcagagatgtggtcgaatggctcggatggtttgaagatcaaaacagaaaaccgactcctgatatgatgcagtatgcgaaaagagcagtcatgtcactgcaaggcctaagagagaagacaattggcaagtatgctaagtcagaatttgacaaatga', window = 30)
normal_ling = linguistic_complexity('atgtctgttacagtcaagagaatcattgacaacacagtcatagttccaaaacttcctgcaaatgaggatccagtggaatacccggcagattacttcagaaaatcaaaggagattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaacagctacttgtatggagcattaaaggacatccggggtaagttggataaagattggtcaagtttcggaataaacatcgggaaagcaggggatacaatcggaatatttgaccttgtatccttgaaagccctggacggcgtacttccagatggagtatcggatgcttccagaaccagcgcagatgacaaatggttgcctttgtatctacttggcttatacagagtgggcagaacacaaatgcctgaatacagaaaaaagctcatggatgggctgacaaatcaatgcaaaatgatcaatgaacagtttgaacctcttgtgccagaaggtcgtgacatttttgatgtgtggggaaatgacagtaattacacaaaaattgtcgctgcagtggacatgttcttccacatgttcaaaaaacatgaatgtgcctcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgcaaaataaccggaatgtctacagaagatgtaacgacctggatcttgaaccgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgacaaggccgattcatacatgccttatttgatcgactttggattgtcttctaagtctccatattcttccgtcaaaaaccctgccttccacttctgggggcaattgacagctcttctgctcagatccaccagagcaaggaatgcccgacagcctgatgacattgagtatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgccgacttggcacaacagttttgtgttggagataacaaatacactccagatgatagtaccggaggattgacgactaatgcaccgccacaaggcagagatgtggtcgaatggctcggatggtttgaagatcaaaacagaaaaccgactcctgatatgatgcagtatgcgaaaagagcagtcatgtcactgcaaggcctaagagagaagacaattggcaagtatgctaagtcagaatttgacaaatga', window = 30)

high_pur = purine_percent('atgtgatgatgagtcaaaagaatcattgataatacagtcatagttccaaaacttcctgcaaatgaagatccagtcgaatacccggcagattacttcagaaaatcaaaagaaattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaattcttacttgtatggagcattaaaagatatcagaggtaaattggataaagattggtcaagtttcggaataaatatcgggaaagcaggggatacaatcggaatatttgatcttgtatccttgaaagcattggatggcgtacttccagatggagtatcggatgcttccagaacatctgcagatgataaatggttgcctttgtatctacttggcttatacagagtcggcagaacacaaatgcctgaatacagaaaaaaactcatggatgggttgacaaatcaatgtaaaatgatcaatgaacaatttgaacctcttgtcccagaaggtcgtgatatttttgatgtctggggaaatgatagtaattacacaaaaattgtcgctgcagtcgatatgttcttccacatgttcaaaaaacatgaatgtgcatcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgtaaaataacaggaatgtctacagaagatgtaacgacatggatcttgaatcgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgataaagcagattcttacatgccttatttgatcgattttggattgtcttctaaatctccatattcttccgtcaaaaatcctgcattccacttctggggacaattgacagctcttttgctcagatccacaagagcaaggaatgcacgacaacctgatgatattgaatatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgcagatttggcacaacaattttgtgttggagataataaatacactccagatgatagtacaggaggattgacgactaatgcaccgccacaaggcagagatgtcgttgaatggctcggatggttcgaagatcaaaatagaaaaccgactcctgatatgatgcaatatgcgaaaagagcagtcatgtcattgcaaggcctaagagaaaaaacaattggcaaatatgctaaatcagaatttgataaatga', window = 30)
high_kap = kappa_index('atgtgatgatgagtcaaaagaatcattgataatacagtcatagttccaaaacttcctgcaaatgaagatccagtcgaatacccggcagattacttcagaaaatcaaaagaaattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaattcttacttgtatggagcattaaaagatatcagaggtaaattggataaagattggtcaagtttcggaataaatatcgggaaagcaggggatacaatcggaatatttgatcttgtatccttgaaagcattggatggcgtacttccagatggagtatcggatgcttccagaacatctgcagatgataaatggttgcctttgtatctacttggcttatacagagtcggcagaacacaaatgcctgaatacagaaaaaaactcatggatgggttgacaaatcaatgtaaaatgatcaatgaacaatttgaacctcttgtcccagaaggtcgtgatatttttgatgtctggggaaatgatagtaattacacaaaaattgtcgctgcagtcgatatgttcttccacatgttcaaaaaacatgaatgtgcatcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgtaaaataacaggaatgtctacagaagatgtaacgacatggatcttgaatcgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgataaagcagattcttacatgccttatttgatcgattttggattgtcttctaaatctccatattcttccgtcaaaaatcctgcattccacttctggggacaattgacagctcttttgctcagatccacaagagcaaggaatgcacgacaacctgatgatattgaatatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgcagatttggcacaacaattttgtgttggagataataaatacactccagatgatagtacaggaggattgacgactaatgcaccgccacaaggcagagatgtcgttgaatggctcggatggttcgaagatcaaaatagaaaaccgactcctgatatgatgcaatatgcgaaaagagcagtcatgtcattgcaaggcctaagagaaaaaacaattggcaaatatgctaaatcagaatttgataaatga', window = 30)
high_ling = linguistic_complexity('atgtgatgatgagtcaaaagaatcattgataatacagtcatagttccaaaacttcctgcaaatgaagatccagtcgaatacccggcagattacttcagaaaatcaaaagaaattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaattcttacttgtatggagcattaaaagatatcagaggtaaattggataaagattggtcaagtttcggaataaatatcgggaaagcaggggatacaatcggaatatttgatcttgtatccttgaaagcattggatggcgtacttccagatggagtatcggatgcttccagaacatctgcagatgataaatggttgcctttgtatctacttggcttatacagagtcggcagaacacaaatgcctgaatacagaaaaaaactcatggatgggttgacaaatcaatgtaaaatgatcaatgaacaatttgaacctcttgtcccagaaggtcgtgatatttttgatgtctggggaaatgatagtaattacacaaaaattgtcgctgcagtcgatatgttcttccacatgttcaaaaaacatgaatgtgcatcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgtaaaataacaggaatgtctacagaagatgtaacgacatggatcttgaatcgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgataaagcagattcttacatgccttatttgatcgattttggattgtcttctaaatctccatattcttccgtcaaaaatcctgcattccacttctggggacaattgacagctcttttgctcagatccacaagagcaaggaatgcacgacaacctgatgatattgaatatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgcagatttggcacaacaattttgtgttggagataataaatacactccagatgatagtacaggaggattgacgactaatgcaccgccacaaggcagagatgtcgttgaatggctcggatggttcgaagatcaaaatagaaaaccgactcctgatatgatgcaatatgcgaaaagagcagtcatgtcattgcaaggcctaagagaaaaaacaattggcaaatatgctaaatcagaatttgataaatga', window = 30)

low_pur = purine_percent('atgtgatgatgagtcaagcgcatcattgacaacaccgtcatagttccaaagcttcccgccaacgaggacccagtggaatacccggccgactacttccgcaagtcaaaggagattcccctttacatcaacactaccaagagtctgtcagacctacgcgggtatgtctaccagggcctcaagtccgggaacgtatcaatcatacatgtcaacagctacctgtatggggccttaaaggacatccggggtaagctggacaaggactggtcaagtttcgggataaacatcgggaaggccggggacaccatcgggatatttgaccttgtatccctgaaggccctggacggcgtacttccagacggggtatcggacgcttcccgcaccagcgccgacgacaagtggttgcccctgtatctacttggcttataccgcgtgggccgcacccagatgcccgagtaccgcaagaagctcatggacgggctgaccaaccagtgcaagatgatcaacgagcagtttgagccccttgtgccagagggtcgtgacatttttgacgtgtgggggaacgacagtaactacaccaagattgtcgctgccgtggacatgttcttccacatgttcaagaagcatgagtgcgcctcgttccgctacgggactattgtttcccgcttcaaggactgcgctgccctggccacctttgggcacctctgcaagataaccgggatgtctaccgaggacgtaacgacctggatcctgaaccgagaggttgccgacgagatggtccagatgatgcttccaggccaggagattgacaaggccgactcatacatgccctatctgatcgactttgggctgtcttctaagtctccctattcttccgtcaagaaccccgccttccacttctgggggcagctgaccgctcttctgctccgctccacccgcgccaggaacgcccgacagcccgacgacattgagtatacctctcttactaccgccggtctgctgtacgcttatgccgtagggtcctctgccgacctggcccagcagttttgcgttggggacaacaagtacactcccgacgacagtaccggggggctgacgactaacgccccgccacagggccgcgacgtggttgaatggctcggatggttcgaggaccagaaccgcaagccgactcccgacatgatgcagtatgcgaagcgcgccgtcatgtcactgcagggcctacgcgagaagaccattggcaagtatgctaagtcagagtttgacaagtga', window = 30)
low_kap = kappa_index('atgtgatgatgagtcaagcgcatcattgacaacaccgtcatagttccaaagcttcccgccaacgaggacccagtggaatacccggccgactacttccgcaagtcaaaggagattcccctttacatcaacactaccaagagtctgtcagacctacgcgggtatgtctaccagggcctcaagtccgggaacgtatcaatcatacatgtcaacagctacctgtatggggccttaaaggacatccggggtaagctggacaaggactggtcaagtttcgggataaacatcgggaaggccggggacaccatcgggatatttgaccttgtatccctgaaggccctggacggcgtacttccagacggggtatcggacgcttcccgcaccagcgccgacgacaagtggttgcccctgtatctacttggcttataccgcgtgggccgcacccagatgcccgagtaccgcaagaagctcatggacgggctgaccaaccagtgcaagatgatcaacgagcagtttgagccccttgtgccagagggtcgtgacatttttgacgtgtgggggaacgacagtaactacaccaagattgtcgctgccgtggacatgttcttccacatgttcaagaagcatgagtgcgcctcgttccgctacgggactattgtttcccgcttcaaggactgcgctgccctggccacctttgggcacctctgcaagataaccgggatgtctaccgaggacgtaacgacctggatcctgaaccgagaggttgccgacgagatggtccagatgatgcttccaggccaggagattgacaaggccgactcatacatgccctatctgatcgactttgggctgtcttctaagtctccctattcttccgtcaagaaccccgccttccacttctgggggcagctgaccgctcttctgctccgctccacccgcgccaggaacgcccgacagcccgacgacattgagtatacctctcttactaccgccggtctgctgtacgcttatgccgtagggtcctctgccgacctggcccagcagttttgcgttggggacaacaagtacactcccgacgacagtaccggggggctgacgactaacgccccgccacagggccgcgacgtggttgaatggctcggatggttcgaggaccagaaccgcaagccgactcccgacatgatgcagtatgcgaagcgcgccgtcatgtcactgcagggcctacgcgagaagaccattggcaagtatgctaagtcagagtttgacaagtga', window = 30)
low_ling = linguistic_complexity('atgtgatgatgagtcaagcgcatcattgacaacaccgtcatagttccaaagcttcccgccaacgaggacccagtggaatacccggccgactacttccgcaagtcaaaggagattcccctttacatcaacactaccaagagtctgtcagacctacgcgggtatgtctaccagggcctcaagtccgggaacgtatcaatcatacatgtcaacagctacctgtatggggccttaaaggacatccggggtaagctggacaaggactggtcaagtttcgggataaacatcgggaaggccggggacaccatcgggatatttgaccttgtatccctgaaggccctggacggcgtacttccagacggggtatcggacgcttcccgcaccagcgccgacgacaagtggttgcccctgtatctacttggcttataccgcgtgggccgcacccagatgcccgagtaccgcaagaagctcatggacgggctgaccaaccagtgcaagatgatcaacgagcagtttgagccccttgtgccagagggtcgtgacatttttgacgtgtgggggaacgacagtaactacaccaagattgtcgctgccgtggacatgttcttccacatgttcaagaagcatgagtgcgcctcgttccgctacgggactattgtttcccgcttcaaggactgcgctgccctggccacctttgggcacctctgcaagataaccgggatgtctaccgaggacgtaacgacctggatcctgaaccgagaggttgccgacgagatggtccagatgatgcttccaggccaggagattgacaaggccgactcatacatgccctatctgatcgactttgggctgtcttctaagtctccctattcttccgtcaagaaccccgccttccacttctgggggcagctgaccgctcttctgctccgctccacccgcgccaggaacgcccgacagcccgacgacattgagtatacctctcttactaccgccggtctgctgtacgcttatgccgtagggtcctctgccgacctggcccagcagttttgcgttggggacaacaagtacactcccgacgacagtaccggggggctgacgactaacgccccgccacagggccgcgacgtggttgaatggctcggatggttcgaggaccagaaccgcaagccgactcccgacatgatgcagtatgcgaagcgcgccgtcatgtcactgcagggcctacgcgagaagaccattggcaagtatgctaagtcagagtttgacaagtga', window = 30)
  

normal_pur_mean = np.mean(normal_pur)
normal_kap_mean = np.mean(normal_kap)
normal_ling_mean = np.mean(normal_ling)

high_pur_mean = np.mean(high_pur)
high_kap_mean = np.mean(high_kap)
high_ling_mean = np.mean(high_ling)

low_pur_mean = np.mean(low_pur)
low_kap_mean = np.mean(low_kap)
low_ling_mean = np.mean(low_ling)


f, ax = plt.subplots(1, 3, sharex = True, sharey = True)
 
ax[0].scatter(x = normal_pur, y = normal_kap, marker = 'o', c='gray', alpha = 0.5)
ax[0].scatter(x =normal_pur_mean, y = normal_kap_mean, marker = 'o', c = 'black', s = 50, alpha = 0.5)
ax[0].errorbar(x = normal_pur_mean, y = normal_kap_mean, xerr = np.std(normal_pur), yerr = np.std(normal_kap), fmt = 'black')
ax[0].scatter(x =high_pur_mean, y = high_kap_mean, marker = 'o', c = 'red', s = 50, alpha = 0.5)
ax[0].scatter(x =low_pur_mean, y = low_kap_mean, marker = 'o', c = 'blue', s = 50, alpha = 0.5)
ax[0].set_ylabel('Kappa Index of Coincidence')
ax[0].set_title('Normal Codon')
ax[1].scatter(x = high_pur, y = high_kap, marker = 'o', c = 'red', alpha = 0.5)
ax[1].scatter(x =high_pur_mean, y = high_kap_mean, marker = 'o', c = 'black', s = 50, alpha = 0.5)
ax[1].set_xlabel('(A+G)%')
ax[1].set_title('High Codon')
ax[2].scatter(x = low_pur, y = low_kap, marker = 'o', c = 'blue', alpha = 0.5)
ax[2].scatter(x =low_pur_mean, y = low_kap_mean, marker = 'o', c = 'black', s = 50, alpha = 0.5)
ax[2].set_title('Low Codon')
plt.show()
            
            
        
            
    
    
    

        
        
        