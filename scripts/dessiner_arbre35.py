import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import Phylo

RES = "/mnt/c/Bioinfo_RAM/Results/"
refs = {
    "GCF_000005845.2": "K-12 MG1655",
    "GCF_000010385.1": "SE11",
    "GCF_000013265.1": "UTI89",
    "GCF_000019645.1": "SMS-3-5",
    "GCF_000210475.1": "ETEC H10407",
}

tab = pd.read_csv(RES + "panel30_tableau_final.tsv", sep="\t", dtype=str).fillna("")
categorie = {a.strip(): c.strip() for a, c in zip(tab["accession"], tab["categorie"])}
st = {a.strip(): s.strip() for a, s in zip(tab["accession"], tab["ST"])}
for acc in refs:
    categorie[acc] = "référence"

def etiquette(acc):
    if acc in refs:
        return refs[acc]
    s = f"ST{st[acc]}" if st.get(acc) else "ST ?"
    return f"{categorie[acc]} | {acc[4:]} | {s}"

# Couleurs par origine
cats = sorted(set(categorie.values()) - {"référence"})
palette = plt.get_cmap("tab10")
couleur = {c: palette(i) for i, c in enumerate(cats)}
couleur["référence"] = "black"

tree = Phylo.read(RES + "ecoli35_core.treefile", "newick")
tree.root_at_midpoint()
tree.ladderize()

# Couleur de chaque étiquette (clé = texte affiché)
col_etiq = {etiquette(t.name): couleur[categorie[t.name]] for t in tree.get_terminals()}

xmax = max(tree.depths().values())
fig, ax = plt.subplots(figsize=(11, 10))
Phylo.draw(
    tree, axes=ax, do_show=False,
    label_func=lambda c: etiquette(c.name) if c.is_terminal() else None,
    label_colors=lambda lab: col_etiq.get(lab, "black"),
    branch_labels=lambda c: None,
)
ax.set_xlim(-0.05 * xmax, 2.0 * xmax)
ax.set_yticks([])
ax.set_ylabel("")
ax.set_xlabel("substitutions par site")
ax.set_title("Arbre du core genome (35 génomes d'E. coli, enraciné au point médian, sans supports)")
ax.legend(handles=[Patch(color=couleur[c], label=c) for c in cats + ["référence"]],
          loc="lower left", frameon=False)
fig.savefig(RES + "arbre35_origine.png", dpi=200, bbox_inches="tight")
print("Image enregistrée : Results/arbre35_origine.png")

# --- Test : le plus proche voisin a-t-il la même origine ? (30 génomes du panel) ---
tips = [t.name for t in tree.get_terminals() if t.name not in refs]
n = len(tips)
D = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        D[i, j] = D[j, i] = tree.distance(tips[i], tips[j])
np.fill_diagonal(D, np.inf)
voisin = D.argmin(axis=1)
labels = np.array([categorie[t] for t in tips])
obs = np.mean(labels[voisin] == labels)

rng = np.random.default_rng(42)
perm = []
for _ in range(2000):
    p = rng.permutation(labels)
    perm.append(np.mean(p[voisin] == p))
perm = np.array(perm)
pval = (np.sum(perm >= obs) + 1) / (len(perm) + 1)

print(f"Plus proche voisin de même origine : {100*obs:.0f} % observé, "
      f"{100*perm.mean():.0f} % attendu au hasard (p = {pval:.3f}, {n} génomes)")
