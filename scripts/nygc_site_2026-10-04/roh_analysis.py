import pandas as pd, numpy as np, json, warnings
from scipy.stats import spearmanr, mannwhitneyu
warnings.filterwarnings("ignore")
D = "C:/work/alsu/nygc_site_2026-10-04/"
h = pd.read_csv(D + "UZB_expanded_ROH.hom.indiv", sep=r"\s+")
fam = pd.read_csv("C:/work/alsu/figures_2026-10-04/adm_in.fam", sep=r"\s+", header=None, usecols=[0, 1], names=["FID", "IID"])
Q = pd.read_csv("C:/work/alsu/figures_2026-10-04/adm_in.4.Q", sep=r"\s+", header=None); Q.columns = ["SAS", "EAS", "AFR", "EUR"]
d = pd.concat([fam, Q], axis=1); d = d[d.FID.astype(str) == "0"]
m = d.merge(h[["IID", "NSEG", "KB"]], on="IID", how="inner")
m["FROH"] = m.KB / 2881033
out = dict(n=int(len(m)), froh_mean=float(m.FROH.mean()), froh_median=float(m.FROH.median()), n_gt_0625=int((m.FROH > 0.0625).sum()))
out["spearman"] = {c: dict(rho=float(spearmanr(m[c], m.FROH)[0]), p=float(spearmanr(m[c], m.FROH)[1])) for c in ["EUR", "EAS", "SAS"]}
tert = {}
for c in ["EUR", "EAS", "SAS"]:
    q = pd.qcut(m[c], 3, labels=["low", "mid", "high"])
    g = m.groupby(q).FROH.agg(["count", "median", "mean"])
    tert[c] = {k: dict(n=int(g.loc[k, "count"]), median=float(g.loc[k, "median"]), mean=float(g.loc[k, "mean"])) for k in g.index}
    tert[c]["bounds"] = [float(m[c].quantile(1 / 3)), float(m[c].quantile(2 / 3))]
out["tertiles"] = tert
# is the consanguineous tail enriched for an ancestry? (FROH > 0.0625 vs rest)
hi = m[m.FROH > 0.0625]; lo = m[m.FROH <= 0.0625]
out["tail"] = {c: dict(hi=float(hi[c].mean()), lo=float(lo[c].mean()), p=float(mannwhitneyu(hi[c], lo[c]).pvalue)) for c in ["EUR", "EAS", "SAS"]}
# U-shape test on each component: quadratic term
quad = {}
for c in ["EUR", "EAS", "SAS"]:
    x = m[c].values; y = m.FROH.values
    A = np.vstack([np.ones_like(x), x - x.mean(), (x - x.mean()) ** 2]).T
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ b; s2 = (res ** 2).sum() / (len(y) - 3); cov = s2 * np.linalg.inv(A.T @ A)
    quad[c] = dict(lin=float(b[1]), quad=float(b[2]), t_quad=float(b[2] / np.sqrt(cov[2, 2])))
out["quad"] = quad
# within self-reported Uzbeks only (phenotype sheet nationality code 1), to avoid mixing nationality groups
import csv
rows = list(csv.reader(open("C:/work/alsu/GWAS от 27.08 - For_Plink.csv", encoding="utf-8-sig")))
sh = pd.DataFrame([r[:12] if len(r) >= 12 else r[:12] + [""] * (12 - len(r)) for r in rows[1:]]).rename(columns={0: "IID", 8: "nat"})[["IID", "nat"]]
sh["IID"] = sh.IID.str.strip(); sh["nat"] = sh.nat.astype(str).str.strip()
mu = m.merge(sh, on="IID", how="left"); mu = mu[mu.nat == "1"]
out["uzbek_only"] = dict(n=int(len(mu)), **{c: dict(rho=float(spearmanr(mu[c], mu.FROH)[0]), p=float(spearmanr(mu[c], mu.FROH)[1])) for c in ["EUR", "EAS", "SAS"]})
R = json.load(open(D + "results.json")); R["roh"] = out; json.dump(R, open(D + "results.json", "w"), indent=1)
print(json.dumps(out, indent=1)[:3500])
