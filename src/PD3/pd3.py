from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
import os
from math import log
import numpy as np
import pandas as pd

def gap_function(x, y):
    if y==0:
        return 0
    elif y==1:
        return -2
    return -(2 + y/4.0 + log(y)/2.0)

# read sequences
dir = "data/PD3_data"
filenames = os.listdir(dir)
names = [f.replace(".fasta", "") for f in filenames]
filepaths = [os.path.join(dir, filename) for filename in filenames]
sequences = []
for f in filepaths:
    sequences.append(SeqIO.read(f, "fasta").seq)

# alignment example
# alignment = pairwise2.align.globalmc(sequences[0], sequences[1], 5, -4, gap_function, gap_function)

# get score
# print(alignment[0].score)
#
# # print info
# for a in alignment:
#     print(format_alignment(*a))


## bogdan

def distmat(sequences):
    dist = np.array([[pairwise2.align.globalxs(s1, s2, -0.5, -0.05)[0].score for s1 in sequences] for s2 in sequences])
    return -dist + dist.max()


def UPGMA(dm, names=None):
    assert dm.shape[0] == dm.shape[1]
    if not names:
        names = [i for i in range(dm.shape[0])]

    h = [([i, names[i]], 1) for i in range(dm.shape[0])]

    for it in range(dm.shape[0] - 1):
        n = dm.shape[0]
        amax = np.argmin(dm + np.diag(np.ones(n) * float('inf')))

        i = amax // n
        j = amax % n

        # Update distance matrix
        dm[i] = (h[i][1] * dm[i] + h[j][1] * dm[j]) / (h[i][1] + h[j][1])
        dm[:, i] = dm[i]

        dm = dm[(np.arange(n) != j), :][:, (np.arange(n) != j)]

        # Update structure
        h[i] = ([h[i][0], h[j][0]], h[i][1] + h[j][1])
        h.pop(j)

    return h[0][0]


distance_matrix = distmat(sequences)
tree = UPGMA(distance_matrix, names)
tree


def Fitch_Margoliash(dm, names=None):
    tree = UPGMA(dm, names)
    df = pd.DataFrame(columns=['dist', 'node1', 'node2'])

    def plot_tree(tree, dm,  node_number, df):
        A, B = tree
        node_name = str(node_number)
        if type(A[0]) is not int:
            node_number += 1
            A_indices = plot_tree(A, dm, node_number, df)
            A_name = str(node_number)
        else:
            A_indices = [A[0]]
            A_name = A[1]
        if type(B[0]) is not int:
            node_number += 1
            B_indices = plot_tree(B, dm, node_number, df)
            B_name = str(node_number)
        else:
            B_indices = [B[0]]
            B_name = B[1]

        C_indices = list(set(range(dm.shape[0])) - set(A_indices) - set(B_indices))
        AB_dist = dm[A_indices][:, B_indices].mean()
        if len(C_indices)==0:
            print(AB_dist, A_name, B_name)
            df.loc[df.shape[0]] = [AB_dist, A_name, B_name]
        else:
            BC_dist = dm[B_indices][:, C_indices].mean()
            AC_dist = dm[A_indices][:, C_indices].mean()
            a = (AB_dist + AC_dist - BC_dist)/2
            b = (AB_dist + BC_dist - AC_dist) / 2
            c = (BC_dist + AC_dist - AB_dist) / 2
            print(a, A_name, node_name)
            print(b, B_name, node_name)
            df.loc[df.shape[0]] = [a, A_name, node_name]
            df.loc[df.shape[0]] = [b, B_name, node_name]
        return A_indices + B_indices
    plot_tree(tree, dm, 0, df)

    return df


df1 = Fitch_Margoliash(distance_matrix, names)

