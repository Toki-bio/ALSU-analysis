import pandas as pd, numpy as np, json, csv, io
from scipy.stats import kruskal
D = "C:/work/alsu/nygc_site_2026-10-04/"
rows = list(csv.reader(open("C:/work/alsu/GWAS от 27.08 - For_Plink.csv", encoding="utf-8-sig")))
hdr = rows[0]
sheet = pd.DataFrame([r[:12] + [""] * (12 - len(r)) if len(r) < 12 else r[:12] for r in rows[1:]])
sheet = sheet.rename(columns={0: "IID", 3: "birth", 8: "nat"})[["IID", "birth", "nat"]]
sheet["IID"] = sheet.IID.str.strip()
fam = pd.read_csv("C:/work/alsu/figures_2026-10-04/adm_in.fam", sep=r"\s+", header=None, usecols=[0, 1], names=["FID", "IID"])
Q = pd.read_csv("C:/work/alsu/figures_2026-10-04/adm_in.4.Q", sep=r"\s+", header=None)
Q.columns = ["SAS", "EAS", "AFR", "EUR"]
d = pd.concat([fam, Q], axis=1)
d["IID2"] = d.IID.replace({"Rustamov_Aziz": "GWAS2026_0035"})
coh = d[d.FID.astype(str) == "0"].copy()
ref = d[d.FID.astype(str) != "0"]
m = coh.merge(sheet, on="IID", how="left", indicator=True)
print("cohort", len(coh), "matched to phenotype sheet by ID:", int((m._merge == "both").sum()))
names = {"1": "Uzbek", "2": "Kazakh", "3": "Tajik", "4": "Russian", "5": "Tatar", "6": "Korean", "7": "Karakalpak", "8": "Other"}
m = m[m._merge == "both"].copy()
m["nat"] = m.nat.astype(str).str.strip()
m = m[m.nat.isin(names)]
m["natn"] = m.nat.map(names)
comps = ["EUR", "EAS", "SAS"]
tab = []
for g, s in m.groupby("natn"):
    tab.append(dict(group=g, n=int(len(s)), **{c: [float(s[c].mean()), float(s[c].std())] for c in comps}))
tab.sort(key=lambda r: -r["n"])
# Kruskal-Wallis by nationality, groups with n>=5
big = [g for g, s in m.groupby("natn") if len(s) >= 5]
def kw(df, col, by):
    gs = [s[col].values for g, s in df.groupby(by) if len(s) >= 5]
    return float(kruskal(*gs).pvalue)
res = {c: kw(m[m.natn.isin(big)], c, "natn") for c in comps}
# control 1: label scramble (same group sizes), 200 permutations
rng = np.random.default_rng(1)
scr = {c: [] for c in comps}
mb = m[m.natn.isin(big)].copy()
for _ in range(200):
    mb["perm"] = rng.permutation(mb.natn.values)
    for c in comps:
        scr[c].append(kw(mb, c, "perm"))
scr_summary = {c: dict(median_p=float(np.median(v)), frac_below_0p05=float(np.mean(np.array(v) < .05))) for c, v in scr.items()}
# control 2: positive control: NYGC reference group labels (from mapping) vs components
pm = pd.read_csv("C:/work/ALSU-analysis/data/pop_mapping.txt", sep=r"\s+", header=None, names=["IID", "lab"], engine="python")
rr = ref.merge(pm, on="IID", how="left"); rr = rr[rr.lab.isin(["EUR", "EAS", "SAS", "AFR"])]
pos = {c: float(kruskal(*[s[c].values for g, s in rr.groupby("lab")]).pvalue) for c in comps}
# birthplace within self-reported Uzbeks
u = m[m.natn == "Uzbek"].copy()
u["birth"] = u.birth.astype(str).str.strip()
bp = u.groupby("birth").size(); keep = bp[bp >= 10].index.tolist()
ub = u[u.birth.isin(keep) & (u.birth != "")]
btab = []
for g, s in ub.groupby("birth"):
    btab.append(dict(group=g, n=int(len(s)), **{c: [float(s[c].mean()), float(s[c].std())] for c in comps}))
btab.sort(key=lambda r: -r["n"])
bres = {c: kw(ub, c, "birth") for c in comps}
out = dict(n_cohort=int(len(coh)), n_matched=int(len(m)), nat_table=tab, kw_nat=res, scramble=scr_summary, positive_control=pos, birth_table=btab, kw_birth=bres, n_birth_groups=len(keep), n_birth_samples=int(len(ub)))
json.dump(out, open(D + "eth_results.json", "w"), indent=1, ensure_ascii=False)
R = json.load(open(D + "results.json")); R["eth"] = out; json.dump(R, open(D + "results.json", "w"), indent=1)
print(json.dumps({k: out[k] for k in ["n_matched", "kw_nat", "scramble", "positive_control", "kw_birth", "n_birth_groups", "n_birth_samples"]}, indent=1))
for r in tab: print(r["group"], r["n"], [round(r[c][0], 3) for c in comps])
for r in btab[:15]: print("birth", r["group"], r["n"], [round(r[c][0], 3) for c in comps])
