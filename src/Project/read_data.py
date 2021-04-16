######################
# Import modules
from Bio import Phylo, SeqIO, AlignIO, Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import networkx, pylab, matplotlib
import os
import pandas as pd
from Bio.Phylo import PhyloXML

# get metadata
metadata = pd.read_csv("data/project/poland/1618582302224.metadata.tsv", delimiter="\t")
# get_data
input_file = 'data/project/1618507098376.sequences.fasta'
records = SeqIO.parse(input_file, 'fasta')
records = list(records) # make a copy, otherwise our generator
                        # is exhausted after calculating maxlen
maxlen = max(len(record.seq) for record in records)

# pad sequences so that they all have the same length
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '-')
        record.seq = Seq.Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

# write to temporary file and do alignment
output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')

aln = AlignIO.read(output_file, "fasta")

for i, a in enumerate(aln):
    a.id = metadata.iloc[i].date + " " + metadata.iloc[i].division + " " + a.id

# Read the sequences and align



# Print the alignment
print(aln)

# Calculate the distance matrix
calculator = DistanceCalculator('identity')
calculator2 = DistanceCalculator('blosum62')
dm = calculator.get_distance(aln)

# Print the distance Matrix
print('\nDistance Matrix\n===================')
print(dm)

# Construct the phylogenetic tree using UPGMA algorithm
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)
tree2 = constructor.nj(dm)
tree_xml = tree.as_phyloxml()

L = tree.get_terminals()
for clade in tree.get_terminals():
    if clade.name in L:
        clade.color = 'red'
# Draw the phylogenetic tree
Phylo.draw(tree)
Phylo.draw(tree2)
Phylo.draw_graphviz(tree)
net = Phylo.to_networkx(tree)
networkx.draw(net)

# Print the phylogenetic tree in the terminal
print('\nPhylogenetic Tree\n===================')
Phylo.draw_ascii(tree)


def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

def branch_length(t):
    t = t.clades
    if t[0] == None:
        t[0] = min(branch_length(t[1]), branch_length(t[2]))/2
        t[1][0] -= t[0]
        t[2][0] -= t[0]
        return t[0]
    else:
        return t[0]


###############
import time
aln = AlignIO.read('data/project/poland/1618582302224.sequences_padded.fasta', 'fasta')
aln2 = aln[:100]
calculator = DistanceCalculator('identity')
t1 = time.time()
dm = calculator.get_distance(aln2)
t2 = time.time()
print(t2-t1)