import pandas as pd, numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = "C:/work/alsu/nygc_site_2026-10-04/"; K = D + "kaz/"; I = "C:/work/ALSU-analysis/images/"
fam = pd.read_csv(K + "panel_kaz.fam", sep=r"\s+", header=None, usecols=[0, 1], names=["FID", "IID"], dtype=str)
lab = pd.read_csv(K + "labels.txt", sep=r"\s+", header=None, names=["FID", "IID", "grp"], dtype=str)
fam = fam.merge(lab[["IID", "grp"]], on="IID", how="left")
assert fam.grp.notna().all() and len(fam) == 4192
out = {"n": fam.grp.value_counts().to_dict()}
# CV errors
cv = {}
for k in range(2, 8):
    for l in open(K + f"adm_K{k}.log"):
        if l.startswith("CV error"):
            cv[k] = float(l.split(":")[1])
out["cv"] = cv
order = ["UZB", "KAZ", "EUR", "SAS", "EAS", "AFR"]
comp = {}
Qs = {}
for k in range(2, 8):
    Q = pd.read_csv(K + f"panel_kaz.{k}.Q", sep=r"\s+", header=None).values
    Qs[k] = Q
    m = pd.DataFrame(Q).groupby(fam.grp.values).mean().reindex(order)
    comp[k] = m.round(3).values.tolist()
out["comp_means"] = comp
# K=4: name components by the four reference groups
Q4 = Qs[4]; m4 = pd.DataFrame(Q4).groupby(fam.grp.values).mean()
names4 = {c: m4[c].loc[["EUR", "SAS", "EAS", "AFR"]].idxmax() for c in range(4)}
assert len(set(names4.values())) == 4, names4
out["k4_names"] = {int(c): n for c, n in names4.items()}
out["k4_means"] = {g: {names4[c]: float(m4.loc[g, c]) for c in range(4)} for g in order}
# K=5: which component is specific for UZB/KAZ?
Q5 = Qs[5]; m5 = pd.DataFrame(Q5).groupby(fam.grp.values).mean()
out["k5_means"] = {g: [float(x) for x in m5.loc[g]] for g in order}
ref_max = {c: m5[c].loc[["EUR", "SAS", "EAS", "AFR"]].idxmax() for c in range(5)}
out["k5_ref_dominant"] = {int(c): n for c, n in ref_max.items()}
# spread of KAZ vs UZB at K=4
for g in ("UZB", "KAZ"):
    s = pd.DataFrame(Q4[fam.grp.values == g])
    out[f"k4_sd_{g}"] = {names4[c]: float(s[c].std()) for c in range(4)}
# PCA
ev = [float(x) for x in open(K + "pca_kaz.eigenval").read().split()]
e = pd.read_csv(K + "pca_kaz.eigenvec", sep=r"\s+", header=None, names=["FID", "IID"] + [f"PC{i}" for i in range(1, 21)], dtype={"FID": str, "IID": str})
e = e.merge(fam[["IID", "grp"]], on="IID")
out["pca_means"] = {g: {f"PC{i}": [float(e[e.grp == g][f"PC{i}"].mean()), float(e[e.grp == g][f"PC{i}"].std())] for i in range(1, 5)} for g in order}
out["pca_pct10"] = [100 * x / sum(ev[:10]) for x in ev[:6]]
# nearest-reference label for Kazakh and Uzbek samples (k=25, PC1-4)
ref = e[e.grp.isin(["EUR", "SAS", "EAS", "AFR"])]
X = ref[["PC1", "PC2", "PC3", "PC4"]].values; L = ref.grp.values
def knn(g):
    Y = e[e.grp == g][["PC1", "PC2", "PC3", "PC4"]].values
    import collections
    r = collections.Counter()
    for y in Y:
        idx = np.argsort(((X - y) ** 2).sum(1))[:25]; r[collections.Counter(L[idx]).most_common(1)[0][0]] += 1
    return dict(r)
out["knn_KAZ"] = knn("KAZ"); out["knn_UZB"] = knn("UZB")
# FST table
f = pd.read_csv(K + "fst_pairs.tsv", sep="\t")
out["fst"] = {r.pair: dict(w=float(r.weighted_fst), m=float(r.mean_fst)) for r in f.itertuples()}
# is KAZ samples' PC4 offset like the Uzbek one?
json.dump(out, open(D + "kaz_results.json", "w"), indent=1)
R = json.load(open(D + "results.json")); R["kaz"] = out; json.dump(R, open(D + "results.json", "w"), indent=1)
# figures
col = {"UZB": "#2e7d32", "KAZ": "#6a1b9a", "EUR": "#7570b3", "SAS": "#e7298a", "EAS": "#1b9e77", "AFR": "#d95f02"}
fig, axs = plt.subplots(1, 3, figsize=(13, 4), dpi=130)
for ax, (a, b) in zip(axs, [("PC1", "PC2"), ("PC2", "PC3"), ("PC3", "PC4")]):
    for g in ["EUR", "SAS", "EAS", "AFR", "UZB", "KAZ"]:
        s = e[e.grp == g]; ax.scatter(s[a], s[b], s=5 if g in ("EUR", "SAS", "EAS", "AFR") else 7, c=col[g], alpha=.45 if g in ("EUR", "SAS", "EAS", "AFR", "UZB") else .9, label=f"{g} (n={len(s)})")
    ax.set_xlabel(a); ax.set_ylabel(b)
axs[0].legend(fontsize=7, markerscale=2)
fig.suptitle("PCA of Uzbek, Kazakh (GSA) and NYGC reference samples, 34,619 SNPs", fontsize=10); fig.tight_layout(); fig.savefig(I + "nygc_kaz_pca.png"); plt.close()
# bars K=4 and K=5
fig, axs = plt.subplots(2, 6, figsize=(13, 5), dpi=130, sharey=True, gridspec_kw={"width_ratios": [1, 1, 1, 1, 1, 1]})
cmap = ["#7570b3", "#e7298a", "#1b9e77", "#d95f02", "#9e9e9e", "#26a69a", "#ffa000"]
for row, k in enumerate((4, 5)):
    Q = Qs[k]
    for ax, g in zip(axs[row], order):
        s = Q[fam.grp.values == g]
        o = np.argsort(-s[:, int(np.argmax(s.mean(0)))]) if g in ("UZB", "KAZ") else np.arange(len(s))
        s = s[o]; bot = np.zeros(len(s))
        for c in range(k):
            ax.bar(range(len(s)), s[:, c], bottom=bot, width=1, color=cmap[c], linewidth=0); bot += s[:, c]
        ax.set_title(f"{g} (n={len(s)})", fontsize=8); ax.set_xticks([]); ax.set_xlim(-.5, len(s) - .5)
    axs[row][0].set_ylabel(f"K={k}")
fig.suptitle("ADMIXTURE with Kazakh samples added (unsupervised; component colours are arbitrary per K)", fontsize=10); fig.tight_layout(); fig.savefig(I + "nygc_kaz_admixture.png"); plt.close()
print(json.dumps({k: out[k] for k in ["n", "cv", "k4_means", "k5_means", "k5_ref_dominant", "knn_KAZ", "knn_UZB", "pca_means"]}, indent=1)[:4500])
print({g: out["fst"][g] for g in out["fst"] if "KAZ" in g})
