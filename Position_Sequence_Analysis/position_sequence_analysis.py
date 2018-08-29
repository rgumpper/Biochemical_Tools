# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 11:17:41 2018

@author: Ryan
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import seaborn as sns

def word_cloud_seq(sequence, window = 9, start = 0):
    seq = sequence.upper()
    seq = seq.replace(' ', '')
    seq = seq.replace('T', 'U')
    position_seq = [[] for i in range(window)]
    for i in range(start, len(seq)-window, window):
        win_seq = seq[i:i+window]
        for x in range(window):
            position_seq[x].append(win_seq[x])
    position = [[] for i in range(window)]
    for i in range(window):
        A = position_seq[i].count('A')
        G = position_seq[i].count('G')
        C = position_seq[i].count('C')
        U = position_seq[i].count('U')
        A = A/len(position_seq[i])
        G = G/len(position_seq[i])
        C = C/len(position_seq[i])
        U = U/len(position_seq[i])
        position[i].append(A)
        position[i].append(G)
        position[i].append(C)
        position[i].append(U)
    position = np.array(position)
    return position

        
seq = ['atgtctgttacagtcaagagaatcattgacaacacagtcatagttccaaaacttcctgcaaatgaggatccagtggaatacccggcagattacttcagaaaatcaaaggagattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaacagctacttgtatggagcattaaaggacatccggggtaagttggataaagattggtcaagtttcggaataaacatcgggaaagcaggggatacaatcggaatatttgaccttgtatccttgaaagccctggacggcgtacttccagatggagtatcggatgcttccagaaccagcgcagatgacaaatggttgcctttgtatctacttggcttatacagagtgggcagaacacaaatgcctgaatacagaaaaaagctcatggatgggctgacaaatcaatgcaaaatgatcaatgaacagtttgaacctcttgtgccagaaggtcgtgacatttttgatgtgtggggaaatgacagtaattacacaaaaattgtcgctgcagtggacatgttcttccacatgttcaaaaaacatgaatgtgcctcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgcaaaataaccggaatgtctacagaagatgtaacgacctggatcttgaaccgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgacaaggccgattcatacatgccttatttgatcgactttggattgtcttctaagtctccatattcttccgtcaaaaaccctgccttccacttctgggggcaattgacagctcttctgctcagatccaccagagcaaggaatgcccgacagcctgatgacattgagtatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgccgacttggcacaacagttttgtgttggagataacaaatacactccagatgatagtaccggaggattgacgactaatgcaccgccacaaggcagagatgtggtcgaatggctcggatggtttgaagatcaaaacagaaaaccgactcctgatatgatgcagtatgcgaaaagagcagtcatgtcactgcaaggcctaagagagaagacaattggcaagtatgctaagtcagaatttgacaaatga',
        'atgtgatgatgagtcaaaagaatcattgataatacagtcatagttccaaaacttcctgcaaatgaagatccagtcgaatacccggcagattacttcagaaaatcaaaagaaattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaattcttacttgtatggagcattaaaagatatcagaggtaaattggataaagattggtcaagtttcggaataaatatcgggaaagcaggggatacaatcggaatatttgatcttgtatccttgaaagcattggatggcgtacttccagatggagtatcggatgcttccagaacatctgcagatgataaatggttgcctttgtatctacttggcttatacagagtcggcagaacacaaatgcctgaatacagaaaaaaactcatggatgggttgacaaatcaatgtaaaatgatcaatgaacaatttgaacctcttgtcccagaaggtcgtgatatttttgatgtctggggaaatgatagtaattacacaaaaattgtcgctgcagtcgatatgttcttccacatgttcaaaaaacatgaatgtgcatcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgtaaaataacaggaatgtctacagaagatgtaacgacatggatcttgaatcgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgataaagcagattcttacatgccttatttgatcgattttggattgtcttctaaatctccatattcttccgtcaaaaatcctgcattccacttctggggacaattgacagctcttttgctcagatccacaagagcaaggaatgcacgacaacctgatgatattgaatatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgcagatttggcacaacaattttgtgttggagataataaatacactccagatgatagtacaggaggattgacgactaatgcaccgccacaaggcagagatgtcgttgaatggctcggatggttcgaagatcaaaatagaaaaccgactcctgatatgatgcaatatgcgaaaagagcagtcatgtcattgcaaggcctaagagaaaaaacaattggcaaatatgctaaatcagaatttgataaatga',
        'atgtgatgatgagtcaagcgcatcattgacaacaccgtcatagttccaaagcttcccgccaacgaggacccagtggaatacccggccgactacttccgcaagtcaaaggagattcccctttacatcaacactaccaagagtctgtcagacctacgcgggtatgtctaccagggcctcaagtccgggaacgtatcaatcatacatgtcaacagctacctgtatggggccttaaaggacatccggggtaagctggacaaggactggtcaagtttcgggataaacatcgggaaggccggggacaccatcgggatatttgaccttgtatccctgaaggccctggacggcgtacttccagacggggtatcggacgcttcccgcaccagcgccgacgacaagtggttgcccctgtatctacttggcttataccgcgtgggccgcacccagatgcccgagtaccgcaagaagctcatggacgggctgaccaaccagtgcaagatgatcaacgagcagtttgagccccttgtgccagagggtcgtgacatttttgacgtgtgggggaacgacagtaactacaccaagattgtcgctgccgtggacatgttcttccacatgttcaagaagcatgagtgcgcctcgttccgctacgggactattgtttcccgcttcaaggactgcgctgccctggccacctttgggcacctctgcaagataaccgggatgtctaccgaggacgtaacgacctggatcctgaaccgagaggttgccgacgagatggtccagatgatgcttccaggccaggagattgacaaggccgactcatacatgccctatctgatcgactttgggctgtcttctaagtctccctattcttccgtcaagaaccccgccttccacttctgggggcagctgaccgctcttctgctccgctccacccgcgccaggaacgcccgacagcccgacgacattgagtatacctctcttactaccgccggtctgctgtacgcttatgccgtagggtcctctgccgacctggcccagcagttttgcgttggggacaacaagtacactcccgacgacagtaccggggggctgacgactaacgccccgccacagggccgcgacgtggttgaatggctcggatggttcgaggaccagaaccgcaagccgactcccgacatgatgcagtatgcgaagcgcgccgtcatgtcactgcagggcctacgcgagaagaccattggcaagtatgctaagtcagagtttgacaagtga']

for i in seq:
    normal = word_cloud_seq(i)     
    x_axis = [1,2,3,4,5,6,7,8,9]
    y_axis = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    width = 0.4
    f = plt.figure()
    p1 = plt.bar(x_axis, normal[:, 0], width, color = 'r')
    p2 = plt.bar(x_axis, normal[:, 1], width, bottom = normal[:, 0], color = 'b')
    p3 = plt.bar(x_axis, normal[:, 2], width, bottom = normal[:, 0] + normal[:, 1], color = 'g')
    p4 = plt.bar(x_axis, normal[:, 3], width, bottom = normal[:, 0] + normal[:, 1] + normal[:, 2], color = 'c' )
    plt.ylabel('Nucleotide Percent')
    plt.xlabel('Window Position')
    plt.xticks(x_axis)
    plt.yticks(y_axis)
    plt.xlim([0, 11])
    plt.title('Nucleotide Percent vs Window Position')
    plt.legend((p1[0], p2[0], p3[0], p4[0]), ('A', 'G', 'C', 'U'))
    plt.show()
    
norm = 'atgtctgttacagtcaagagaatcattgacaacacagtcatagttccaaaacttcctgcaaatgaggatccagtggaatacccggcagattacttcagaaaatcaaaggagattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaacagctacttgtatggagcattaaaggacatccggggtaagttggataaagattggtcaagtttcggaataaacatcgggaaagcaggggatacaatcggaatatttgaccttgtatccttgaaagccctggacggcgtacttccagatggagtatcggatgcttccagaaccagcgcagatgacaaatggttgcctttgtatctacttggcttatacagagtgggcagaacacaaatgcctgaatacagaaaaaagctcatggatgggctgacaaatcaatgcaaaatgatcaatgaacagtttgaacctcttgtgccagaaggtcgtgacatttttgatgtgtggggaaatgacagtaattacacaaaaattgtcgctgcagtggacatgttcttccacatgttcaaaaaacatgaatgtgcctcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgcaaaataaccggaatgtctacagaagatgtaacgacctggatcttgaaccgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgacaaggccgattcatacatgccttatttgatcgactttggattgtcttctaagtctccatattcttccgtcaaaaaccctgccttccacttctgggggcaattgacagctcttctgctcagatccaccagagcaaggaatgcccgacagcctgatgacattgagtatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgccgacttggcacaacagttttgtgttggagataacaaatacactccagatgatagtaccggaggattgacgactaatgcaccgccacaaggcagagatgtggtcgaatggctcggatggtttgaagatcaaaacagaaaaccgactcctgatatgatgcagtatgcgaaaagagcagtcatgtcactgcaaggcctaagagagaagacaattggcaagtatgctaagtcagaatttgacaaatga'
hi = 'atgtgatgatgagtcaaaagaatcattgataatacagtcatagttccaaaacttcctgcaaatgaagatccagtcgaatacccggcagattacttcagaaaatcaaaagaaattcctctttacatcaatactacaaaaagtttgtcagatctaagaggatatgtctaccaaggcctcaaatccggaaatgtatcaatcatacatgtcaattcttacttgtatggagcattaaaagatatcagaggtaaattggataaagattggtcaagtttcggaataaatatcgggaaagcaggggatacaatcggaatatttgatcttgtatccttgaaagcattggatggcgtacttccagatggagtatcggatgcttccagaacatctgcagatgataaatggttgcctttgtatctacttggcttatacagagtcggcagaacacaaatgcctgaatacagaaaaaaactcatggatgggttgacaaatcaatgtaaaatgatcaatgaacaatttgaacctcttgtcccagaaggtcgtgatatttttgatgtctggggaaatgatagtaattacacaaaaattgtcgctgcagtcgatatgttcttccacatgttcaaaaaacatgaatgtgcatcgttcagatacggaactattgtttccagattcaaagattgtgctgcattggcaacatttggacacctctgtaaaataacaggaatgtctacagaagatgtaacgacatggatcttgaatcgagaagttgcagatgaaatggtccaaatgatgcttccaggccaagaaattgataaagcagattcttacatgccttatttgatcgattttggattgtcttctaaatctccatattcttccgtcaaaaatcctgcattccacttctggggacaattgacagctcttttgctcagatccacaagagcaaggaatgcacgacaacctgatgatattgaatatacatctcttactacagcaggtttgttgtacgcttatgcagtaggatcctctgcagatttggcacaacaattttgtgttggagataataaatacactccagatgatagtacaggaggattgacgactaatgcaccgccacaaggcagagatgtcgttgaatggctcggatggttcgaagatcaaaatagaaaaccgactcctgatatgatgcaatatgcgaaaagagcagtcatgtcattgcaaggcctaagagaaaaaacaattggcaaatatgctaaatcagaatttgataaatga'
lo = 'atgtgatgatgagtcaagcgcatcattgacaacaccgtcatagttccaaagcttcccgccaacgaggacccagtggaatacccggccgactacttccgcaagtcaaaggagattcccctttacatcaacactaccaagagtctgtcagacctacgcgggtatgtctaccagggcctcaagtccgggaacgtatcaatcatacatgtcaacagctacctgtatggggccttaaaggacatccggggtaagctggacaaggactggtcaagtttcgggataaacatcgggaaggccggggacaccatcgggatatttgaccttgtatccctgaaggccctggacggcgtacttccagacggggtatcggacgcttcccgcaccagcgccgacgacaagtggttgcccctgtatctacttggcttataccgcgtgggccgcacccagatgcccgagtaccgcaagaagctcatggacgggctgaccaaccagtgcaagatgatcaacgagcagtttgagccccttgtgccagagggtcgtgacatttttgacgtgtgggggaacgacagtaactacaccaagattgtcgctgccgtggacatgttcttccacatgttcaagaagcatgagtgcgcctcgttccgctacgggactattgtttcccgcttcaaggactgcgctgccctggccacctttgggcacctctgcaagataaccgggatgtctaccgaggacgtaacgacctggatcctgaaccgagaggttgccgacgagatggtccagatgatgcttccaggccaggagattgacaaggccgactcatacatgccctatctgatcgactttgggctgtcttctaagtctccctattcttccgtcaagaaccccgccttccacttctgggggcagctgaccgctcttctgctccgctccacccgcgccaggaacgcccgacagcccgacgacattgagtatacctctcttactaccgccggtctgctgtacgcttatgccgtagggtcctctgccgacctggcccagcagttttgcgttggggacaacaagtacactcccgacgacagtaccggggggctgacgactaacgccccgccacagggccgcgacgtggttgaatggctcggatggttcgaggaccagaaccgcaagccgactcccgacatgatgcagtatgcgaagcgcgccgtcatgtcactgcagggcctacgcgagaagaccattggcaagtatgctaagtcagagtttgacaagtga'


normal = word_cloud_seq(norm)
high = word_cloud_seq(hi)
low = word_cloud_seq(lo)


normal_df = pd.DataFrame(data = normal, index = x_axis, columns = ['A', 'G', 'C', 'U'])
high_df = pd.DataFrame(data = high, index = x_axis, columns = ['A', 'G', 'C', 'U'])
low_df = pd.DataFrame(data = low, index = x_axis, columns = ['A', 'G', 'C', 'U'])

g = plt.figure('Normal')
sns.heatmap(normal_df, annot = True, vmin = 0.1, vmax = 0.40)
plt.ylabel('Window Position')
plt.xlabel('Nucleotide')
plt.title('Nucleotide Percentage by Position')
plt.show()

k = plt.figure('High')
sns.heatmap(high_df, annot = True, vmin = 0.1, vmax = 0.40)
plt.ylabel('Window Position')
plt.xlabel('Nucleotide')
plt.title('Nucleotide Percentage by Position')
plt.show()

h = plt.figure('Low')
sns.heatmap(low_df, annot = True, vmin = 0.1, vmax = 0.40)
plt.ylabel('Window Position')
plt.xlabel('Nucleotide')
plt.title('Nucleotide Percentage by Position')
plt.show()




