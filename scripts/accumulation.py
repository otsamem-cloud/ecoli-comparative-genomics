import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PA, PNG = sys.argv[1], sys.argv[2]
K12, R = "GCF_000005845.2", 200
rng = np.random.default_rng(42)  # graine fixe : résultat reproductible

df = pd.read_csv(PA, sep="\t")
genomes = [c for c in df.columns if c.startswith("GCF_")]
M = df[genomes].to_numpy() > 0   # familles x génomes (présence/absence)
N = M.shape[1]

pan = np.empty((R, N), dtype=int)
core = np.empty((R, N), dtype=int)
for r in range(R):
    P = M[:, rng.permutation(N)]
    pan[r] = np.logical_or.accumulate(P, axis=1).sum(axis=0)    # union cumulée
    core[r] = np.logical_and.accumulate(P, axis=1).sum(axis=0)  # intersection cumulée

k = np.arange(1, N + 1)
pan_m, core_m = pan.mean(axis=0), core.mean(axis=0)
nouveaux = np.diff(pan_m)  # familles nouvelles apportées par le k-ième génome (k>=2)

# Loi de Heaps : nouveaux(k) = kappa * k^(-alpha)  ->  droite en log-log
pente, _ = np.polyfit(np.log(k[1:]), np.log(nouveaux), 1)
alpha = -pente

print(f"Pangénome : {pan_m[4]:.0f} familles à 5 génomes -> {pan_m[-1]:.0f} à {N}")
print(f"Core strict : {core_m[4]:.0f} à 5 génomes -> {core_m[-1]:.0f} à {N}")
print(f"Nouvelles familles apportées : {nouveaux[0]:.0f} par le 2e génome, {nouveaux[-1]:.0f} par le {N}e")
print(f"Loi de Heaps : alpha = {alpha:.2f} ({'pangénome ouvert' if alpha < 1 else 'pangénome fermé'})")

uniques = M[(M.sum(axis=1) == 1)].sum(axis=0)
i = genomes.index(K12)
rang = int((uniques > uniques[i]).sum()) + 1
print(f"K-12 : {uniques[i]} familles uniques, rang {rang}/{N} (1 = le plus de familles uniques)")

fig, ax = plt.subplots(figsize=(8, 5))
for arr, label in ((pan, "Pangénome"), (core, "Core strict")):
    lo, hi = np.percentile(arr, [5, 95], axis=0)
    ax.plot(k, arr.mean(axis=0), marker="o", ms=3, label=label)
    ax.fill_between(k, lo, hi, alpha=0.25)
ax.set_xlabel("Nombre de génomes")
ax.set_ylabel("Nombre de familles de gènes")
ax.set_title(f"Courbes d'accumulation ({R} ordres aléatoires)")
ax.legend()
fig.savefig(PNG, dpi=200, bbox_inches="tight")
print("Image enregistrée :", PNG)
