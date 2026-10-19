import sys
import numpy as np
import pandas as pd
from Bio import Phylo

RES = "/mnt/c/Bioinfo_RAM/Results/"
seuils = [float(x) for x in sys.argv[1:]] or [0.0, 0.0001, 0.001]

tab = pd.read_csv(RES + "panel30_tableau_final.tsv", sep="\t", dtype=str).fillna("")
cat = dict(zip(tab["accession"].str.strip(), tab["categorie"].str.strip()))

tree = Phylo.read(RES + "ecoli35_core.treefile", "newick")
tips = [t.name for t in tree.get_terminals() if t.name in cat]   # les 30 du panel
n = len(tips)
D = np.zeros((n, n))     # distances sur l'arbre du core
for i in range(n):
    for j in range(i + 1, n):
        D[i, j] = D[j, i] = tree.distance(tips[i], tips[j])
labels_all = np.array([cat[t] for t in tips])

# Distance de Jaccard sur les familles accessoires (shell + cloud)
df = pd.read_csv(RES + "presence_absence_35.tsv", sep="\t")
acc = df[df["classe"].isin(["shell", "cloud"])]
M = (acc[tips].to_numpy() > 0).astype(int)
inter = M.T @ M
taille = M.sum(axis=0)
J = 1 - inter / (taille[:, None] + taille[None, :] - inter)

def groupes(t):
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

def test(idx, dist, n_perm=500):
    d = dist[np.ix_(idx, idx)].copy()
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
    res = {"core": [], "accessoire": []}
    for _ in range(200):
        idx = np.array([rng.choice(m) for m in g])
        res["core"].append(test(idx, D))
        res["accessoire"].append(test(idx, J))
    for nom, r in res.items():
        r = np.array(r)
        print(f"  {nom:<11}: {100*np.median(r[:, 0]):.0f} % observé, "
              f"{100*np.median(r[:, 1]):.0f} % attendu, "
              f"p<0,05 dans {100*np.mean(r[:, 2] < 0.05):.0f} % des tirages")
