from itertools import groupby
from collections import Counter
from pathlib import Path

CLU = Path("/mnt/c/Bioinfo_RAM/Data/comparatif/clusters_cluster.tsv")
OUT = Path("/mnt/c/Bioinfo_RAM/Results/presence_absence.tsv")

def lire_paires(chemin):
    """Génère (représentant, génome) pour chaque protéine, en streaming."""
    with open(chemin) as f:
        for ligne in f:
            rep, membre = ligne.rstrip("\n").split("\t")
            yield rep, membre.split("|")[0]

# Les lignes d'une même famille sont consécutives -> groupby suffit
familles, genomes = [], set()
for rep, groupe in groupby(lire_paires(CLU), key=lambda x: x[0]):
    copies = Counter(g for _, g in groupe)  # génome -> nb de copies
    familles.append((rep, copies))
    genomes.update(copies)
genomes = sorted(genomes)

classes, specifiques = Counter(), Counter()
with open(OUT, "w") as f:
    f.write("famille\t" + "\t".join(genomes) + "\tn_genomes\tclasse\n")
    for rep, copies in familles:
        n = len(copies)
        if n == len(genomes):
            classe = "core"
        elif n == 1:
            classe = "specifique"
            specifiques[next(iter(copies))] += 1
        else:
            classe = "accessoire"
        classes[classe] += 1
        f.write(rep + "\t" + "\t".join(str(copies.get(g, 0)) for g in genomes)
                + f"\t{n}\t{classe}\n")

print(f"Familles totales : {len(familles)}")
for c in ("core", "accessoire", "specifique"):
    print(f"  {c:<11}: {classes[c]}")
print("Familles spécifiques par génome :")
for g in genomes:
    print(f"  {g} : {specifiques[g]}")
