import pandas as pd, numpy as np, json, collections, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = "C:/work/alsu/nygc_site_2026-10-04/"; I = "C:/work/ALSU-analysis/images/"; S = "C:/work/ALSU-analysis/"
# DRAGEN pca.eigenvec (plink 1.9, no header): FID IID PC1..PC20; same order as adm_in.fam
cols = ["FID", "IID"] + [f"PC{i}" for i in range(1, 21)]
e = pd.read_csv(D + "dragen_pca.eigenvec", sep=r"\s+", header=None, names=cols)
ev = [float(x) for x in open(D + "pca.eigenval").read().split()]
pm = pd.read_csv(S + "data/pop_mapping.txt", sep=r"\s+", header=None, names=["IID", "lab"], engine="python")
e["IID2"] = e.IID.replace({"Rustamov_Aziz": "GWAS2026_0035"})
e = e.merge(pm.rename(columns={"IID": "IID2"}), on="IID2", how="left")
e["cohort"] = e.FID.astype(str) == "0"
e["grp"] = np.where(e.cohort, "ALSU", e.lab)
ref = e[~e.cohort]; coh = e[e.cohort].copy()
R = json.load(open(D + "results.json"))
out = {"eigenval": ev, "n_total": int(len(e)), "n_cohort": int(len(coh)), "n_ref": int(len(ref))}
tot10 = sum(ev[:10])
out["pct10"] = [100 * x / tot10 for x in ev[:10]]
# group stats on PC1..PC4
g = {}
for name, d in list(ref.groupby("grp")) + [("ALSU", coh)]:
    g[name] = {f"PC{i}": [float(d[f"PC{i}"].mean()), float(d[f"PC{i}"].std())] for i in range(1, 5)}
    g[name]["n"] = int(len(d))
out["groups"] = g
# kNN assignment of cohort samples (k=25, PC1-4, Euclidean) against the four references
X = ref[[f"PC{i}" for i in range(1, 5)]].values; lab = ref.grp.values
Y = coh[[f"PC{i}" for i in range(1, 5)]].values
k = 25
votes = []
for y in Y:
    d = ((X - y) ** 2).sum(1); idx = np.argsort(d)[:k]
    c = collections.Counter(lab[idx]); votes.append(c)
coh["knn"] = [c.most_common(1)[0][0] for c in votes]
coh["knn_frac"] = [c.most_common(1)[0][1] / k for c in votes]
cnt = coh.knn.value_counts().to_dict(); out["knn_counts"] = {a: int(b) for a, b in cnt.items()}
out["knn_pure"] = {a: int(((coh.knn == a) & (coh.knn_frac == 1.0)).sum()) for a in cnt}
# mixed neighbours: samples whose k nearest refs are from >= 2 groups
out["knn_mixed"] = int((coh.knn_frac < 1.0).sum())
# distance of each cohort sample from the cohort median in PC1-4 (robust outliers)
med = coh[[f"PC{i}" for i in range(1, 5)]].median().values
dist = np.sqrt(((Y - med) ** 2).sum(1)); coh["dmed"] = dist
mad = np.median(np.abs(dist - np.median(dist))) * 1.4826
coh["z"] = (dist - np.median(dist)) / mad
out["outlier_z5"] = int((coh.z > 5).sum())
o = coh.sort_values("dmed", ascending=False).head(12)
out["far_cohort"] = [dict(id=("GWAS2026_0035" if r.IID2 == "GWAS2026_0035" else f"ALSU_{i:04d}"), pc1=float(r.PC1), pc2=float(r.PC2), knn=r.knn, frac=float(r.knn_frac)) for i, (_, r) in enumerate(o.iterrows())]
# AFR-like cohort members and EAS-like by self-reported label counts
tab = pd.crosstab(coh.lab.fillna("unlabelled"), coh.knn)
out["label_by_knn"] = {a: {b: int(tab.loc[a, b]) for b in tab.columns} for a in tab.index}
# pc1/pc2 spread for the cohort
out["coh_pc12_sd"] = [float(coh.PC1.std()), float(coh.PC2.std())]
out["afr_like_cohort"] = int((coh.PC1 > 0.01).sum())
out["afr_like_ids_pc"] = [[float(r.PC1), float(r.PC2)] for _, r in coh[coh.PC1 > 0.01].iterrows()]
json.dump(out, open(D + "pca_results.json", "w"), indent=1)
R["pca"] = out; json.dump(R, open(D + "results.json", "w"), indent=1)
# figures
colr = {"AFR": "#d95f02", "EAS": "#1b9e77", "EUR": "#7570b3", "SAS": "#e7298a"}
def panel(ax, a, b, zoom=False):
    for p, c in colr.items():
        s = ref[ref.grp == p]; ax.scatter(s[a], s[b], s=5, c=c, alpha=.5, label=p)
    ax.scatter(coh[a], coh[b], s=6, c="k", alpha=.6, label="ALSU")
    if zoom: ax.set_xlim(-.0115, -.004)
    ax.set_xlabel(f"{a} ({100 * ev[int(a[2:]) - 1] / tot10:.1f}% of first 10)"); ax.set_ylabel(f"{b} ({100 * ev[int(b[2:]) - 1] / tot10:.1f}%)")
fig, axs = plt.subplots(2, 2, figsize=(10, 8), dpi=130)
panel(axs[0, 0], "PC1", "PC2"); panel(axs[0, 1], "PC1", "PC2", zoom=True); panel(axs[1, 0], "PC3", "PC4"); panel(axs[1, 1], "PC2", "PC3", zoom=False)
axs[0, 0].set_title("PC1 vs PC2", fontsize=10); axs[0, 1].set_title("PC1 vs PC2, Eurasian part enlarged", fontsize=10)
axs[1, 0].set_title("PC3 vs PC4", fontsize=10); axs[1, 1].set_title("PC2 vs PC3", fontsize=10)
axs[0, 0].legend(fontsize=7, markerscale=2)
fig.suptitle("Global PCA, 1,256 Uzbek samples + 2,712 NYGC reference samples, 64,959 LD-pruned SNPs", fontsize=10)
fig.tight_layout(); fig.savefig(I + "nygc_global_pca_panels.png"); plt.close()
fig, ax = plt.subplots(figsize=(6, 3), dpi=130)
ax.bar(range(1, 21), ev, color="#1565c0"); ax.set_yscale("log"); ax.set_xlabel("PC"); ax.set_ylabel("eigenvalue (log)"); ax.set_title("PCA eigenvalues: three large components, then a plateau", fontsize=9)
fig.tight_layout(); fig.savefig(I + "nygc_global_pca_scree.png"); plt.close()
print(json.dumps({k2: out[k2] for k2 in ["pct10", "knn_counts", "knn_pure", "knn_mixed", "outlier_z5", "coh_pc12_sd", "afr_like_cohort", "afr_like_ids_pc"]}, indent=1))
print(out["far_cohort"][:6]); print(json.dumps(out["label_by_knn"], indent=0)[:1500])
