from Bio import pairwise2
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
import os
from math import log


def gap_function(x, y):
    if y==0:
        return 0
    elif y==1:
        return -2
    return -(2 + y/4.0 + log(y)/2.0)

dir = "data/PD3_data"
filenames = os.listdir(dir)
filepaths = [os.path.join(dir, filename) for filename in filenames]
sequences = []
for f in filepaths:
    sequences.append(SeqIO.read(f, "fasta").seq)


alignment = pairwise2.align.globalmc(sequences[0], sequences[1], 5, -4, gap_function, gap_function)

for a in alignment:
    print(format_alignment(*a))