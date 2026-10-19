import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import Phylo

RES = "/mnt/c/Bioinfo_RAM/Results/"
SEUIL = 95   # convention UFBoot : >= 95 = branche bien soutenue
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

cats = sorted(set(categorie.values()) - {"référence"})
palette = plt.get_cmap("tab10")
couleur = {c: palette(i) for i, c in enumerate(cats)}
couleur["référence"] = "black"

# 1. Supports lus sur l'arbre non enraciné, indexés par bipartition
tree = Phylo.read(RES + "ecoli35_boot.treefile", "newick")
tous = frozenset(t.name for t in tree.get_terminals())
ref = sorted(tous)[0]

def cle(noms):
    """Bipartition normalisée : le côté qui ne contient pas le génome de référence."""
    s = frozenset(noms)
    return s if ref not in s else tous - s

supports = {}
for cl in tree.get_nonterminals():
    if cl.confidence is not None:
        supports[cle(t.name for t in cl.get_terminals())] = cl.confidence

# 2. Enracinement au point médian, puis supports ré-attribués par bipartition
tree.root_at_midpoint()
tree.ladderize()
for cl in tree.get_nonterminals():
    if cl is tree.root:
        cl.confidence = None
    else:
        cl.confidence = supports.get(cle(t.name for t in cl.get_terminals()))

# 3. Dessin
col_etiq = {etiquette(t.name): couleur[categorie[t.name]] for t in tree.get_terminals()}
xmax = max(tree.depths().values())
fig, ax = plt.subplots(figsize=(11, 10))
Phylo.draw(
    tree, axes=ax, do_show=False,
    label_func=lambda c: etiquette(c.name) if c.is_terminal() else None,
    label_colors=lambda lab: col_etiq.get(lab, "black"),
    branch_labels=lambda c: (f"{c.confidence:.0f}"
        if c.confidence is not None and c.confidence < SEUIL else None),
)
ax.set_xlim(-0.05 * xmax, 2.0 * xmax)
ax.set_yticks([])
ax.set_ylabel("")
ax.set_xlabel("substitutions par site")
ax.set_title(f"Arbre du core genome (35 génomes, UFBoot 1000) : supports < {SEUIL} affichés")
ax.legend(handles=[Patch(color=couleur[c], label=c) for c in cats + ["référence"]],
          loc="lower left", frameon=False)
fig.savefig(RES + "arbre35_origine_supports.png", dpi=200, bbox_inches="tight")
print("Image enregistrée : Results/arbre35_origine_supports.png")

# 4. Résumé des supports
n = len(supports)
ok = sum(v >= SEUIL for v in supports.values())
print(f"Branches internes : {n} | support >= {SEUIL} : {ok} ({100*ok/n:.0f} %)")
faibles = sorted(((v, k) for k, v in supports.items() if v < SEUIL), key=lambda x: x[0])
for v, k in faibles:
    petit = k if len(k) <= len(tous) - len(k) else tous - k
    print(f"  {v:>4.0f}  isole {len(petit)} génome(s) : "
          + ", ".join(sorted(etiquette(x) for x in petit)))
