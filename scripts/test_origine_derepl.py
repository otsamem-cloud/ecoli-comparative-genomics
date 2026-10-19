import sys
import numpy as np
import pandas as pd
from Bio import Phylo

RES = "/mnt/c/Bioinfo_RAM/Results/"
seuils = [float(x) for x in sys.argv[1:]] or [0.0001, 0.001]

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
labels_all = np.array([cat[t] for t in tips])

def groupes(t):
    """Liaison simple : deux génomes à moins de t sont dans le même groupe."""
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] < t:
                parent[find(i)] = find(j)
    g = {}
    for i in range(n):
        g.setdefault(find(i), []).append(i)
    return list(g.values())

rng = np.random.default_rng(42)

def test(idx, n_perm=500):
    d = D[np.ix_(idx, idx)].copy()
    np.fill_diagonal(d, np.inf)
    nn = d.argmin(axis=1)
    lab = labels_all[idx]
    obs = np.mean(lab[nn] == lab)
    perm = np.empty(n_perm)
    for k in range(n_perm):
        p = rng.permutation(lab)
        perm[k] = np.mean(p[nn] == p)
    return obs, perm.mean(), (np.sum(perm >= obs) + 1) / (n_perm + 1)

for t in seuils:
    g = groupes(t)
    print(f"\nSeuil {t} : {len(g)} groupes indépendants sur {n} génomes")
    for m in g:
        if len(m) > 1:
            print("  groupe redondant :", ", ".join(
                f"{tips[i][4:]} ({cat[tips[i]]}, ST{st[tips[i]] or '?'})" for i in m))
    res = []
    for _ in range(200):
        idx = np.array([rng.choice(m) for m in g])
        res.append(test(idx))
    res = np.array(res)
    print(f"  Voisin de même origine : {100*np.median(res[:, 0]):.0f} % observé (médiane), "
          f"{100*np.median(res[:, 1]):.0f} % attendu au hasard")
    print(f"  p < 0,05 dans {100*np.mean(res[:, 2] < 0.05):.0f} % des 200 tirages de représentants")
