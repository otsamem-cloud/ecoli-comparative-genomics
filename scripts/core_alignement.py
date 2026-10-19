import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO

BASE = Path("/mnt/c/Bioinfo_RAM")
CLU = BASE / "Data/comparatif/clusters_cluster.tsv"
PROT = BASE / "Data/comparatif/all_proteins.faa"
PA = BASE / "Results/presence_absence.tsv"
OUT = BASE / "Results/core_concat.faa"
WORK = Path.home() / "align_core"
WORK.mkdir(exist_ok=True)

# 1. Familles core à copie unique
with open(PA) as f:
    genomes = f.readline().rstrip("\n").split("\t")[1:-2]
    ok = set()
    for ligne in f:
        c = ligne.rstrip("\n").split("\t")
        if c[-1] == "core" and all(x == "1" for x in c[1:-2]):
            ok.add(c[0])

# 2. Pour chaque famille : génome -> identifiant de la protéine
membres = {}
with open(CLU) as f:
    for ligne in f:
        rep, m = ligne.rstrip("\n").split("\t")
        if rep in ok:
            membres.setdefault(rep, {})[m.split("|")[0]] = m
besoin = {m for d in membres.values() for m in d.values()}

# 3. Une seule lecture du gros FASTA, on garde seulement le nécessaire
seqs = {r.id: str(r.seq).replace("U", "X").replace("O", "X").replace("*", "") for r in SeqIO.parse(PROT, "fasta") if r.id in besoin}

# 4. Un alignement MAFFT par famille, 4 en parallèle
familles = sorted(membres)

def aligner(i):
    fam = familles[i]
    entree, sortie = WORK / f"fam{i}.faa", WORK / f"fam{i}.aln"
    with open(entree, "w") as f:
        for g in genomes:
            f.write(f">{g}\n{seqs[membres[fam][g]]}\n")
    with open(sortie, "w") as out:
        subprocess.run(["mafft", "--auto", "--quiet", "--thread", "1",
                        str(entree)], stdout=out, check=True)
    return sortie

with ThreadPoolExecutor(4) as ex:
    sorties = list(ex.map(aligner, range(len(familles))))

# 5. Concaténation : une longue séquence par génome
concat = {g: [] for g in genomes}
for s in sorties:
    for rec in SeqIO.parse(s, "fasta"):
        concat[rec.id].append(str(rec.seq).upper())

with open(OUT, "w") as f:
    for g in genomes:
        f.write(f">{g}\n{''.join(concat[g])}\n")

print(f"Familles alignées : {len(familles)}")
print(f"Longueur de la supermatrice : {len(''.join(concat[genomes[0]]))} positions")
