import numpy as np
import pandas as pd
from Bio import Phylo

RES = "/mnt/c/Bioinfo_RAM/Results/"
tab = pd.read_csv(RES + "panel30_tableau_final.tsv", sep="\t", dtype=str).fillna("")
cat = dict(zip(tab["accession"].str.strip(), tab["categorie"].str.strip()))
st = dict(zip(tab["accession"].str.strip(), tab["ST"].str.strip()))

tree = Phylo.read(RES + "ecoli35_core.treefile", "newick")
tips = [t.name for t in tree.get_terminals() if t.name in cat]   # les 30 du panel
n = len(tips)
D = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        D[i, j] = D[j, i] = tree.distance(tips[i], tips[j])
np.fill_diagonal(D, np.inf)
nn = D.argmin(axis=1)

for d, i in sorted((D[i, nn[i]], i) for i in range(n)):
    a, b = tips[i], tips[nn[i]]
    flag = "MEME ORIGINE" if cat[a] == cat[b] else ""
    print(f"{d:.5f}  {cat[a]:<14}{a[4:]:<12}ST{st[a] or '?':<6}-> "
          f"{cat[b]:<14}{b[4:]:<12}ST{st[b] or '?':<6}{flag}")
