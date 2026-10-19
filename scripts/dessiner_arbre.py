import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

RES = "/mnt/c/Bioinfo_RAM/Results/"
noms = {
    "GCF_000005845.2": "K-12 MG1655",
    "GCF_000010385.1": "SE11",
    "GCF_000013265.1": "UTI89",
    "GCF_000019645.1": "SMS-3-5",
    "GCF_000210475.1": "ETEC H10407",
}

tree = Phylo.read(RES + "ecoli_core_tree.treefile", "newick")

# Supports lus AVANT d'enraciner : on n'affiche que s'ils sont tous identiques
supports = {int(c.confidence) for c in tree.get_nonterminals() if c.confidence is not None}
valeur = str(supports.pop()) if len(supports) == 1 else None
tree.root_at_midpoint()

fig, ax = plt.subplots(figsize=(9, 4))
Phylo.draw(
    tree, axes=ax, do_show=False,
    label_func=lambda c: noms.get(c.name) if c.is_terminal() else None,
    branch_labels=lambda c: valeur if (valeur and not c.is_terminal() and c is not tree.root) else None,
)
ax.set_yticks([])
ax.set_ylabel("")
ax.set_title("Arbre du core genome d'E. coli (enraciné au point médian)")
ax.set_xlabel("substitutions par site")
fig.savefig(RES + "arbre_core_ecoli_enracine.png", dpi=200, bbox_inches="tight")
print("Image enregistrée : Results/arbre_core_ecoli_enracine.png")
