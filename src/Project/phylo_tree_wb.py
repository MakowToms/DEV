from Bio import Phylo, SeqIO, AlignIO, Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import networkx
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
from Bio.Phylo import PhyloXMLIO
import pickle

# get data and add padding
input_file = 'data/project/poland/1618582302224.sequences.fasta'
records = SeqIO.parse(input_file, 'fasta')
records = list(records)
maxlen = max(len(record.seq) for record in records)

# pad sequences so that they all have the same length
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '-')
        record.seq = Seq.Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

# write to temporary file and do alignment


#####
# get metadata
metadata = pd.read_csv("data/project/poland/1618582302224.metadata.tsv", delimiter="\t")
metadata = metadata[metadata.strain.str.startswith("hCoV-19")].reset_index()
metadata = metadata.replace("?", np.nan)
divisions = list(metadata.groupby("division").count().reset_index().sort_values("strain", ascending=False).head(16).division)
metadata = metadata[metadata.division.isin(divisions)]
metadata = metadata.sort_values("date").groupby("division").head(5)
records_new = [records[i] for i in np.array(metadata.index)]

output_file = '{}_new.fasta'.format(os.path.splitext(input_file)[0])
with open(output_file, 'w') as f:
    SeqIO.write(records_new, f, 'fasta')

aln = AlignIO.read('data/project/poland/1618582302224.sequences_new.fasta', 'fasta')

for i, a in enumerate(aln):
    a.id = metadata.iloc[i].date + " " + metadata.iloc[i].division + " " + a.id

print(aln)

# Calculate the distance matrix
# calculator = DistanceCalculator('identity')
# dm = calculator.get_distance(aln)
#
# with open("data/project/dist_mat.p", "wb") as f:
#     pickle.dump(dm, f)

with open("data/project/dist_mat.p", "rb") as f:
    dm = pickle.load(f)
# Construct the phylogenetic tree using UPGMA algorithm
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)
tree_xml = tree.as_phyloxml()
PhyloXMLIO.write(tree_xml, "data/project/tree1.xml")

def labels(c):
    if not c.is_terminal():
        return ""
    else:
        return " ".join(c.name.split()[:2])

fig = plt.figure(figsize=(25, 25), dpi=100)
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes, do_show=False, label_func=labels)
plt.savefig("plots/tree1")
Phylo.draw_ascii(tree)