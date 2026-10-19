import sys
from itertools import groupby
from collections import Counter
from math import ceil
from statistics import median
from pathlib import Path

CLU, OUT = Path(sys.argv[1]), Path(sys.argv[2])
K12 = "GCF_000005845.2"

def lire_paires(chemin):
    """Génère (représentant, génome) pour chaque protéine, en streaming."""
    with open(chemin) as f:
        for ligne in f:
            rep, membre = ligne.rstrip("\n").split("\t")
            yield rep, membre.split("|")[0]

familles, genomes = [], set()
for rep, groupe in groupby(lire_paires(CLU), key=lambda x: x[0]):
    copies = Counter(g for _, g in groupe)
    familles.append((rep, copies))
    genomes.update(copies)
genomes = sorted(genomes)
N = len(genomes)
s_core, s_soft, s_shell = ceil(0.99 * N), ceil(0.95 * N), ceil(0.15 * N)

def classer(n):
    if n >= s_core:  return "core"
    if n >= s_soft:  return "soft_core"
    if n >= s_shell: return "shell"
    return "cloud"

classes, uniques = Counter(), Counter()
core_unique = 0
with open(OUT, "w") as f:
    f.write("famille\t" + "\t".join(genomes) + "\tn_genomes\tcopie_unique\tclasse\n")
    for rep, copies in familles:
        n = len(copies)
        classe = classer(n)
        unique = int(all(v == 1 for v in copies.values()))
        classes[classe] += 1
        if n == 1:
            uniques[next(iter(copies))] += 1
        if classe in ("core", "soft_core") and unique:
            core_unique += 1
        f.write(rep + "\t" + "\t".join(str(copies.get(g, 0)) for g in genomes)
                + f"\t{n}\t{unique}\t{classe}\n")

print(f"Génomes : {N} | seuils (nb de génomes) : core>={s_core}, soft>={s_soft}, shell>={s_shell}")
print(f"Familles totales : {len(familles)}")
for c in ("core", "soft_core", "shell", "cloud"):
    print(f"  {c:<10}: {classes[c]}")
print(f"Core + soft core à copie unique (utilisables pour l'arbre) : {core_unique}")
v = [uniques[g] for g in genomes]
print(f"Familles présentes dans un seul génome : min {min(v)}, médiane {median(v)}, max {max(v)}")
print(f"  K-12 MG1655 : {uniques[K12]}")
