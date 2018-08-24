import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def SeqAGRatio(seq = '', window = 9):
    seq = seq.upper()
    seq = seq.replace(' ', '')
    AGRatio = []
    for i in range(0, (len(seq)-window)):
        score = 0
        nucs = seq[i:i+window]
        for j in nucs:
            if j == 'A':
                score = score +1
            elif j == 'G':
                score = score +1
        ratio = score/window
        AGRatio.append(ratio)
    sequence_place = []
    for i in range(window, (len(seq))):
        sequence_place.append(i)
    mean = []
    avg_score = np.mean(AGRatio)
    for i in range(0, len(sequence_place)):
        mean.append(avg_score)
    
    plt.plot(sequence_place, AGRatio)
    plt.plot(sequence_place, mean)
    plt.ylabel('AG Ratio')
    plt.xlabel('Sequence Number')
    plt.title('Sequence AG Ratio')
    plt.show()
