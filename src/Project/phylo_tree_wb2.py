from Bio import Phylo, SeqIO, AlignIO, Seq
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import networkx
import matplotlib.pyplot as plt
import matplotlib
import os
import pandas as pd
import numpy as np
from Bio.Phylo import PhyloXMLIO
import pickle

clades = ["G", "GH", "GR", "GRY", "GV", "L", "O", "S", "V"]
records = []
clades_names = []
metadata = pd.DataFrame()
for clade in clades:
    input_dir = 'data/project/clade_{}'.format(clade)
    files = os.listdir(input_dir)
    for f in files:
        if f.endswith(".fasta"):
            print(os.path.exists(os.path.join(input_dir, f)))
            recs = SeqIO.parse(os.path.join(input_dir, f), 'fasta')
            r = list(recs)
            clades_names = clades_names + [clade] * len(r)
            records = records + r
        if f.endswith(".tsv"):
            mtdt = pd.read_csv(os.path.join(input_dir, f), delimiter="\t", na_values="?").loc[:, ["date", "division"]]
            metadata = metadata.append(mtdt, ignore_index=True)

metadata = metadata.fillna("")
print(records)

maxlen = max(len(record.seq) for record in records)

# pad sequences so that they all have the same length
for record in records:
    if len(record.seq) != maxlen:
        sequence = str(record.seq).ljust(maxlen, '-')
        record.seq = Seq.Seq(sequence)
assert all(len(record.seq) == maxlen for record in records)

output_file = 'data/project/all_records.fasta'
with open(output_file, 'w') as f:
    SeqIO.write(records, f, 'fasta')

aln = AlignIO.read(output_file, 'fasta')

for i, a in enumerate(aln):
    a.id = clades_names[i] + " " + str(metadata.iloc[i].date) + " " + str(metadata.iloc[i].division) + " " + a.id

print(aln)

# Calculate the distance matrix
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

calculator2 = DistanceCalculator('blastn')
dm2 = calculator2.get_distance(aln)

with open("data/project/dist_mat.p", "wb") as f:
    pickle.dump(dm, f)
with open("data/project/dist_mat_blastn.p", "wb") as f:
    pickle.dump(dm2, f)

with open("data/project/dist_mat.p", "rb") as f:
    dm = pickle.load(f)

constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)
tree_xml = tree.as_phyloxml()


# PhyloXMLIO.write(tree_xml, "data/project/tree1.xml")
from matplotlib import colors as mcolors

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
def labels(c):
    if not c.is_terminal():
        return ""
    else:
        return " ".join(c.name.split(" ")[:3])


def labels_col(c):
    if c.split(" ")[0] == "GRY":
        return colors['darkblue']
    elif c.split(" ")[0] == "GR":
        return colors['blue']
    elif c.split(" ")[0] == "G":
        return colors['seagreen']
    elif c.split(" ")[0] == "GH":
        return colors['green']
    elif c.split(" ")[0] == "S":
        return colors['lime']
    elif c.split(" ")[0] == "O":
        return colors['gold']
    elif c.split(" ")[0] == "GV":
        return colors['darkorange']
    elif c.split(" ")[0] == "V":
        return colors['tomato']
    elif c.split(" ")[0] == "L":
        return colors['red']


matplotlib.rc('font', size=30)
fig = plt.figure(figsize=(35, 25), dpi=100)
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes, do_show=False, label_func=labels, label_colors=labels_col)
plt.axis('off')
plt.savefig("plots/tree3.svg", format="svg", transparent=True)
Phylo.draw_ascii(tree)


