from pathlib import Path
from Bio import SeqIO

DATA = Path("/mnt/c/Bioinfo_RAM/Data")
OUT = Path("/mnt/c/Bioinfo_RAM/Results/comparative_stats.tsv")

# Tous les .fna sauf ceux du panel30
fichiers = sorted(p for p in DATA.rglob("*.fna") if "panel30" not in p.parts)

lignes = []
for fna in fichiers:
    accession = fna.name.split("_ASM")[0]  # ex: GCF_000005845.2
    longueurs, total_gc, total = [], 0, 0
    for rec in SeqIO.parse(fna, "fasta"):
        s = str(rec.seq).upper()
        longueurs.append(len(s))
        total_gc += s.count("G") + s.count("C")
        total += len(s)
    lignes.append((accession, total, round(100 * total_gc / total, 2),
                   len(longueurs), max(longueurs),
                   total - max(longueurs)))

entete = ("accession", "taille_totale_pb", "GC_pourcent",
          "nb_replicons", "chromosome_pb", "plasmides_pb")
OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as f:
    f.write("\t".join(entete) + "\n")
    for l in lignes:
        f.write("\t".join(map(str, l)) + "\n")
        print(*l, sep="\t")
