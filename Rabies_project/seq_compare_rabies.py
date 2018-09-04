# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 17:34:36 2018

@author: Ryan
"""

from Bio import SeqIO
import pandas as pd

def align_to_df(file):
    df_name = {}
    with open(file) as handle:
        for record in SeqIO.parse(handle, 'clustal'):
            sequence = str(record.seq)
            sequence = sequence.replace('', ',')
            sequence = sequence.split(',')
            df_name[record.id] = sequence
    df_name = pd.DataFrame.from_dict(df_name, orient = 'index')
    return df_name

N_seq = align_to_df('N_aligned_rabies.aln')

L_seq = align_to_df('L_aligned_rabies.aln')

P_seq = align_to_df('P_aligned_rabies.aln')



df1 = N_seq[N_seq.index.isin(L_seq.index)]

import matplotlib as plt
import seaborn as sns

N_seq.sort_index(inplace = True)
L_seq.sort_index(inplace = True)
P_seq.sort_index(inplace = True)

n_ind = N_seq.index.tolist()
l_ind = L_seq.index.tolist()
p_ind = P_seq.index.tolist()


def redesign_index(index):
    empty = []
    for i in index:
        seq = i[6:11]
        seq = int(seq)
        empty.append(seq)
    return empty

new_ind_n = redesign_index(n_ind)
new_ind_l = redesign_index(l_ind)
new_ind_p = redesign_index(p_ind)

new_ind_n_compare = []
for i in new_ind_n:
    ind = i + 4
    new_ind_n_compare.append(ind)
    

new_ind_p_compare = []
for i in new_ind_n:
    ind = i +1
    new_ind_p_compare.append(ind)


N_seq_for_l = N_seq.copy(deep = True)
N_seq_for_p = N_seq.copy(deep = True)

N_seq_for_l.index = new_ind_n_compare
N_seq_for_p.index = new_ind_p_compare

L_seq.index = new_ind_l
P_seq.index = new_ind_p
N_seq_for_l = N_seq_for_l.filter(items = L_seq.index, axis = 0)
N_seq_for_p = N_seq_for_p.filter(items = P_seq.index, axis = 0)

L_seq = L_seq.filter(items = N_seq_for_l.index, axis = 0)
P_seq = P_seq.filter(items = N_seq_for_p.index, axis = 0)








