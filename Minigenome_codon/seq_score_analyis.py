# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 16:09:44 2018

@author: Ryan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches



def sequence_score(nucleotide_sequence):
    seq = nucleotide_sequence
    seq = seq.upper()
    seq = seq.replace(' ', '')
    seq = seq.replace('T', 'U')
    score = 0
    seq_score =[]
    for i in seq:
        if i == 'A':
            score = score+1
        elif i == 'G':
            score = score+1
        elif i == 'C':
            score = score-1
        elif i == 'U':
            score = score-1
        seq_score.append(score)
    return seq_score

normal = sequence_score('atgtgatgatgagtcaagagaatcattgacaacacagtcatagttccaaaacttcctgcaaatgaggatccagtggaatacccggcagattacttcagaaaatcaaaggagattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaacagctacttgtatggagcattaaaggacatccggggtaagttggataaagattggtcaagtttcggaataaacatcgggaaagcaggggatacaatcggaatatttgaccttgtatccttgaaagccctggacggcgtacttccagatggagtatcggatgcttccagaaccagcgcagatgacaaatggttgcctttgtatctacttggcttatacagagtgggcagaacacaaatgcctgaatacagaaaaaagctcatggatgggctgacaaatcaatgcaaaatgatcaatgaacagtttgaacctcttgtgccagaaggtcgtgacatttttgatgtgtggggaaatgacagtaattacacaaaaattgtcgctgcagtggacatgttcttccacatgttcaaaaaacatgaatgtgcctcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgcaaaataaccggaatgtctacagaagatgtaacgacctggatcttgaaccgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgacaaggccgattcatacatgccttatttgatcgactttggattgtcttctaagtctccatattcttccgtcaaaaaccctgccttccacttctgggggcaattgacagctcttctgctcagatccaccagagcaaggaatgcccgacagcctgatgacattgagtatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgccgacttggcacaacagttttgtgttggagataacaaatacactccagatgatagtaccggaggattgacgactaatgcaccgccacaaggcagagatgtggtcgaatggctcggatggtttgaagatcaaaacagaaaaccgactcctgatatgatgcagtatgcgaaaagagcagtcatgtcactgcaaggcctaagagagaagacaattggcaagtatgctaagtcagaatttgacaaatga')
high = sequence_score('atgtgatgatgagtcaaaagaatcattgataatacagtcatagttccaaaacttcctgcaaatgaagatccagtcgaatacccggcagattacttcagaaaatcaaaagaaattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaattcttacttgtatggagcattaaaagatatcagaggtaaattggataaagattggtcaagtttcggaataaatatcgggaaagcaggggatacaatcggaatatttgatcttgtatccttgaaagcattggatggcgtacttccagatggagtatcggatgcttccagaacatctgcagatgataaatggttgcctttgtatctacttggcttatacagagtcggcagaacacaaatgcctgaatacagaaaaaaactcatggatgggttgacaaatcaatgtaaaatgatcaatgaacaatttgaacctcttgtcccagaaggtcgtgatatttttgatgtctggggaaatgatagtaattacacaaaaattgtcgctgcagtcgatatgttcttccacatgttcaaaaaacatgaatgtgcatcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgtaaaataacaggaatgtctacagaagatgtaacgacatggatcttgaatcgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgataaagcagattcttacatgccttatttgatcgattttggattgtcttctaaatctccatattcttccgtcaaaaatcctgcattccacttctggggacaattgacagctcttttgctcagatccacaagagcaaggaatgcacgacaacctgatgatattgaatatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgcagatttggcacaacaattttgtgttggagataataaatacactccagatgatagtacaggaggattgacgactaatgcaccgccacaaggcagagatgtcgttgaatggctcggatggttcgaagatcaaaatagaaaaccgactcctgatatgatgcaatatgcgaaaagagcagtcatgtcattgcaaggcctaagagaaaaaacaattggcaaatatgctaaatcagaatttgataaatga')
low = sequence_score('atgtgatgatgagtcaagcgcatcattgacaacaccgtcatagttccaaagcttcccgccaacgaggacccagtggaatacccggccgactacttccgcaagtcaaaggagattcccctttacatcaacactaccaagagtctgtcagacctacgcgggtatgtctaccagggcctcaagtccgggaacgtatcaatcatacatgtcaacagctacctgtatggggccttaaaggacatccggggtaagctggacaaggactggtcaagtttcgggataaacatcgggaaggccggggacaccatcgggatatttgaccttgtatccctgaaggccctggacggcgtacttccagacggggtatcggacgcttcccgcaccagcgccgacgacaagtggttgcccctgtatctacttggcttataccgcgtgggccgcacccagatgcccgagtaccgcaagaagctcatggacgggctgaccaaccagtgcaagatgatcaacgagcagtttgagccccttgtgccagagggtcgtgacatttttgacgtgtgggggaacgacagtaactacaccaagattgtcgctgccgtggacatgttcttccacatgttcaagaagcatgagtgcgcctcgttccgctacgggactattgtttcccgcttcaaggactgcgctgccctggccacctttgggcacctctgcaagataaccgggatgtctaccgaggacgtaacgacctggatcctgaaccgagaggttgccgacgagatggtccagatgatgcttccaggccaggagattgacaaggccgactcatacatgccctatctgatcgactttgggctgtcttctaagtctccctattcttccgtcaagaaccccgccttccacttctgggggcagctgaccgctcttctgctccgctccacccgcgccaggaacgcccgacagcccgacgacattgagtatacctctcttactaccgccggtctgctgtacgcttatgccgtagggtcctctgccgacctggcccagcagttttgcgttggggacaacaagtacactcccgacgacagtaccggggggctgacgactaacgccccgccacagggccgcgacgtggttgaatggctcggatggttcgaggaccagaaccgcaagccgactcccgacatgatgcagtatgcgaagcgcgccgtcatgtcactgcagggcctacgcgagaagaccattggcaagtatgctaagtcagagtttgacaagtga')

seq_place = [] 
for i in range(1,len(normal)+1):
    seq_place.append(i)        
    
    
plt.scatter(seq_place, normal, marker = 'o', c = 'gray', alpha = 0.5)
plt.scatter(seq_place, high, marker = 'o', c = 'red', alpha = 0.5)
plt.scatter(seq_place, low, marker = 'o', c = 'blue', alpha = 0.5)
plt.ylabel('Score')
plt.xlabel('Sequence position')
plt.title('Sequence Score')

def sequence_ratio(nucleotide_seqs, window):
    seq = nucleotide_seqs
    seq = seq.upper()
    seq = seq.replace(' ', '')
    seq = seq.replace('T', 'U')
    ratio = []
    for i in range(0, len(seq)-window):
        working_seq = seq[i:i+window]
        A = working_seq.count('A')
        G = working_seq.count('G')
        C = working_seq.count('C')
        U = working_seq.count('U')
        if C+U == 0:
            calc = (A+G)/1
        else:
            calc = (A+G)/(C+U)
        ratio.append(calc)
    return ratio
        
rat_normal = sequence_ratio('atgtgatgatgagtcaagagaatcattgacaacacagtcatagttccaaaacttcctgcaaatgaggatccagtggaatacccggcagattacttcagaaaatcaaaggagattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaacagctacttgtatggagcattaaaggacatccggggtaagttggataaagattggtcaagtttcggaataaacatcgggaaagcaggggatacaatcggaatatttgaccttgtatccttgaaagccctggacggcgtacttccagatggagtatcggatgcttccagaaccagcgcagatgacaaatggttgcctttgtatctacttggcttatacagagtgggcagaacacaaatgcctgaatacagaaaaaagctcatggatgggctgacaaatcaatgcaaaatgatcaatgaacagtttgaacctcttgtgccagaaggtcgtgacatttttgatgtgtggggaaatgacagtaattacacaaaaattgtcgctgcagtggacatgttcttccacatgttcaaaaaacatgaatgtgcctcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgcaaaataaccggaatgtctacagaagatgtaacgacctggatcttgaaccgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgacaaggccgattcatacatgccttatttgatcgactttggattgtcttctaagtctccatattcttccgtcaaaaaccctgccttccacttctgggggcaattgacagctcttctgctcagatccaccagagcaaggaatgcccgacagcctgatgacattgagtatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgccgacttggcacaacagttttgtgttggagataacaaatacactccagatgatagtaccggaggattgacgactaatgcaccgccacaaggcagagatgtggtcgaatggctcggatggtttgaagatcaaaacagaaaaccgactcctgatatgatgcagtatgcgaaaagagcagtcatgtcactgcaaggcctaagagagaagacaattggcaagtatgctaagtcagaatttgacaaatga', 9)
rat_high = sequence_ratio('atgtgatgatgagtcaaaagaatcattgataatacagtcatagttccaaaacttcctgcaaatgaagatccagtcgaatacccggcagattacttcagaaaatcaaaagaaattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaattcttacttgtatggagcattaaaagatatcagaggtaaattggataaagattggtcaagtttcggaataaatatcgggaaagcaggggatacaatcggaatatttgatcttgtatccttgaaagcattggatggcgtacttccagatggagtatcggatgcttccagaacatctgcagatgataaatggttgcctttgtatctacttggcttatacagagtcggcagaacacaaatgcctgaatacagaaaaaaactcatggatgggttgacaaatcaatgtaaaatgatcaatgaacaatttgaacctcttgtcccagaaggtcgtgatatttttgatgtctggggaaatgatagtaattacacaaaaattgtcgctgcagtcgatatgttcttccacatgttcaaaaaacatgaatgtgcatcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgtaaaataacaggaatgtctacagaagatgtaacgacatggatcttgaatcgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgataaagcagattcttacatgccttatttgatcgattttggattgtcttctaaatctccatattcttccgtcaaaaatcctgcattccacttctggggacaattgacagctcttttgctcagatccacaagagcaaggaatgcacgacaacctgatgatattgaatatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgcagatttggcacaacaattttgtgttggagataataaatacactccagatgatagtacaggaggattgacgactaatgcaccgccacaaggcagagatgtcgttgaatggctcggatggttcgaagatcaaaatagaaaaccgactcctgatatgatgcaatatgcgaaaagagcagtcatgtcattgcaaggcctaagagaaaaaacaattggcaaatatgctaaatcagaatttgataaatga', 9)
rat_low = sequence_ratio('atgtgatgatgagtcaagcgcatcattgacaacaccgtcatagttccaaagcttcccgccaacgaggacccagtggaatacccggccgactacttccgcaagtcaaaggagattcccctttacatcaacactaccaagagtctgtcagacctacgcgggtatgtctaccagggcctcaagtccgggaacgtatcaatcatacatgtcaacagctacctgtatggggccttaaaggacatccggggtaagctggacaaggactggtcaagtttcgggataaacatcgggaaggccggggacaccatcgggatatttgaccttgtatccctgaaggccctggacggcgtacttccagacggggtatcggacgcttcccgcaccagcgccgacgacaagtggttgcccctgtatctacttggcttataccgcgtgggccgcacccagatgcccgagtaccgcaagaagctcatggacgggctgaccaaccagtgcaagatgatcaacgagcagtttgagccccttgtgccagagggtcgtgacatttttgacgtgtgggggaacgacagtaactacaccaagattgtcgctgccgtggacatgttcttccacatgttcaagaagcatgagtgcgcctcgttccgctacgggactattgtttcccgcttcaaggactgcgctgccctggccacctttgggcacctctgcaagataaccgggatgtctaccgaggacgtaacgacctggatcctgaaccgagaggttgccgacgagatggtccagatgatgcttccaggccaggagattgacaaggccgactcatacatgccctatctgatcgactttgggctgtcttctaagtctccctattcttccgtcaagaaccccgccttccacttctgggggcagctgaccgctcttctgctccgctccacccgcgccaggaacgcccgacagcccgacgacattgagtatacctctcttactaccgccggtctgctgtacgcttatgccgtagggtcctctgccgacctggcccagcagttttgcgttggggacaacaagtacactcccgacgacagtaccggggggctgacgactaacgccccgccacagggccgcgacgtggttgaatggctcggatggttcgaggaccagaaccgcaagccgactcccgacatgatgcagtatgcgaagcgcgccgtcatgtcactgcagggcctacgcgagaagaccattggcaagtatgctaagtcagagtttgacaagtga', 9)        
        
    
plt.plot(seq_place[0:1260], rat_normal, c = 'red', alpha = 0.5)
plt.plot(seq_place[0:1260], rat_low, c = 'blue', alpha = 0.5)
plt.plot(seq_place[0:1260], rat_normal, c = 'gray', alpha = 0.5)
plt.ylabel('Purine to Pyrimidine Ratio')
plt.xlabel('Sequence position')
plt.title('Ratio')    




