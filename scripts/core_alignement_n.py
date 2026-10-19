import sys, subprocess, itertools
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO

PA, CLU, PROT, OUT = map(Path, sys.argv[1:5])
WORK = Path.home() / "align_core35"
WORK.mkdir(exist_ok=True)

# 1. Familles core + soft core à copie unique, et liste des génomes
with open(PA) as f:
    genomes = f.readline().rstrip("\n").split("\t")[1:-3]
    ok = set()
    for ligne in f:
        c = ligne.rstrip("\n").split("\t")
        if c[-1] in ("core", "soft_core") and c[-2] == "1":
            ok.add(c[0])

# 2. Pour chaque famille : génome -> identifiant de la protéine
membres = {}
with open(CLU) as f:
    for ligne in f:
        rep, m = ligne.rstrip("\n").split("\t")
        if rep in ok:
            membres.setdefault(rep, {})[m.split("|")[0]] = m
besoin = {m for d in membres.values() for m in d.values()}

def nettoyer(s):
    return s.replace("U", "X").replace("O", "X").replace("*", "")

# 3. Une seule lecture du gros FASTA, on garde seulement le nécessaire
seqs = {r.id: nettoyer(str(r.seq)) for r in SeqIO.parse(PROT, "fasta") if r.id in besoin}

# 4. Un alignement MAFFT par famille, 4 en parallèle
familles = sorted(membres)
compteur = itertools.count(1)

def aligner(i):
    presents = membres[familles[i]]
    entree, sortie = WORK / f"fam{i}.faa", WORK / f"fam{i}.aln"
    with open(entree, "w") as f:
        for g in genomes:
            if g in presents:
                f.write(f">{g}\n{seqs[presents[g]]}\n")
    with open(sortie, "w") as out:
        r = subprocess.run(["mafft", "--auto", "--quiet", "--thread", "1", str(entree)],
                           stdout=out, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"MAFFT a échoué sur {entree} : {r.stderr[-300:]}")
    n = next(compteur)
    if n % 500 == 0:
        print(f"  {n}/{len(familles)} alignements", flush=True)
    return sortie

with ThreadPoolExecutor(4) as ex:
    sorties = list(ex.map(aligner, range(len(familles))))

# 5. Concaténation : une longue séquence par génome, tirets si génome absent
concat = {g: [] for g in genomes}
for s in sorties:
    aln = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(s, "fasta")}
    L = len(next(iter(aln.values())))
    for g in genomes:
        concat[g].append(aln.get(g, "-" * L))

with open(OUT, "w") as f:
    for g in genomes:
        f.write(f">{g}\n{''.join(concat[g])}\n")

total = sum(len("".join(v)) for v in concat.values())
tirets = sum("".join(v).count("-") for v in concat.values())
print(f"Génomes : {len(genomes)} | familles alignées : {len(familles)}")
print(f"Longueur de la supermatrice : {len(''.join(concat[genomes[0]]))} positions")
print(f"Part de tirets : {100 * tirets / total:.2f} %")
