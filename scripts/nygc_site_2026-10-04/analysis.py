import pandas as pd, numpy as np, json, os, urllib.request, time, itertools
D = "C:/work/alsu/nygc_site_2026-10-04/"
R = {}
pops = ["UZB", "SAS", "EUR", "EAS", "AFR"]
wt = {"all": {"AFR_EAS":.166827,"AFR_EUR":.141407,"AFR_SAS":.129399,"AFR_UZB":.131005,"EAS_EUR":.0872392,"EAS_SAS":.057043,"EAS_UZB":.0400238,"EUR_SAS":.0313661,"EUR_UZB":.0152432,"SAS_UZB":.0142574},
      "filt": {"AFR_EAS":.166538,"AFR_EUR":.141102,"AFR_SAS":.129073,"AFR_UZB":.130663,"EAS_EUR":.0874646,"EAS_SAS":.0571663,"EAS_UZB":.0400666,"EUR_SAS":.031257,"EUR_UZB":.0151516,"SAS_UZB":.0141336}}
pairs = {}
for k in wt["all"]:
    f = pd.read_csv(D + f"filt_{k}.fst", sep="\t").dropna(subset=["FST"])
    f["FST"] = f.FST.clip(lower=0)  # not used for mean; keep raw below
    raw = pd.read_csv(D + f"filt_{k}.fst", sep="\t").dropna(subset=["FST"])
    chrs = sorted(raw.CHR.unique())
    m_all = raw.FST.mean()
    # leave-one-chromosome-out jackknife for the mean per-SNP FST (SNP-count weighted)
    loo = []
    for c in chrs:
        s = raw[raw.CHR != c].FST
        loo.append(s.mean())
    loo = np.array(loo); g = len(chrs)
    se = np.sqrt((g - 1) / g * ((loo - loo.mean()) ** 2).sum())
    pairs[k] = dict(weighted_all=wt["all"][k], weighted_filt=wt["filt"][k], mean_filt=float(m_all), mean_se=float(se), n=int(len(raw)))
R["pairs"] = pairs
# matrix (weighted, all 82,744 sites; headline)
M = pd.DataFrame(0.0, index=pops, columns=pops)
for k, v in wt["all"].items():
    a, b = k.split("_"); M.loc[a, b] = M.loc[b, a] = v
R["matrix"] = M.round(4).to_dict()
# classical MDS on sqrt? use FST directly as distance (as page step14)
def cmds(M):
    D2 = M.values ** 2
    n = len(D2); J = np.eye(n) - 1 / n
    B = -0.5 * J @ D2 @ J
    w, v = np.linalg.eigh(B); i = np.argsort(w)[::-1]; w, v = w[i], v[:, i]
    X = v[:, :2] * np.sqrt(np.maximum(w[:2], 0))
    return X, w
X, w = cmds(M)
R["mds"] = {p: [float(X[i, 0]), float(X[i, 1])] for i, p in enumerate(pops)}
R["mds_eig"] = [float(x) for x in w]
R["mds_var_pct"] = [float(100 * w[0] / w[w > 0].sum()), float(100 * w[1] / w[w > 0].sum())]
# per-SNP UZB-EUR
ue = pd.read_csv(D + "filt_EUR_UZB.fst", sep="\t").dropna(subset=["FST"])
pb = pd.read_csv(D + "pbs_all.tsv", sep="\t")
ue = ue.merge(pb[["SNP", "AF_UZB", "AF_EUR", "AF_EAS", "AF_SAS", "AF_AFR"]], on="SNP", how="left")
q = ue.FST.quantile([.5, .75, .9, .95, .99, .999]).to_dict()
R["ue_dist"] = dict(n=int(len(ue)), mean=float(ue.FST.mean()), median=float(ue.FST.median()), max=float(ue.FST.max()),
                    q={str(k): float(v) for k, v in q.items()},
                    gt005=int((ue.FST > .05).sum()), gt01=int((ue.FST > .1).sum()), gt02=int((ue.FST > .2).sum()), gt03=int((ue.FST > .3).sum()))
bins = [0, .005, .01, .02, .03, .05, .1, .2, 1.01]
h = pd.cut(ue.FST.clip(lower=0), bins, right=False).value_counts().sort_index()
R["ue_hist"] = [[str(i), int(c)] for i, c in h.items()]
top = ue.sort_values("FST", ascending=False).head(30).copy()
def genes(ch, pos):
    for _ in range(3):
        try:
            u = f"https://rest.ensembl.org/overlap/region/human/{ch}:{pos}-{pos}?feature=gene;content-type=application/json"
            r = json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "alsu"}), timeout=30))
            g = [x.get("external_name") or x["gene_id"] for x in r if x.get("biotype") in ("protein_coding", "lncRNA", "miRNA", "processed_transcript") or True]
            return sorted(set(g))
        except Exception as e:
            time.sleep(1.5)
    return None
cache_p = D + "gene_cache.json"
cache = json.load(open(cache_p)) if os.path.exists(cache_p) else {}
def gene_of(snp):
    if snp not in cache:
        ch, pos = snp.split(":")[:2]; cache[snp] = genes(ch, pos); time.sleep(.15)
    return cache[snp]
rows = []
for _, r in top.iterrows():
    rows.append(dict(snp=r.SNP, chr=int(r.CHR), pos=int(r.POS), fst=float(r.FST), uzb=float(r.AF_UZB), eur=float(r.AF_EUR), eas=float(r.AF_EAS), sas=float(r.AF_SAS), afr=float(r.AF_AFR), genes=gene_of(r.SNP)))
R["ue_top30"] = rows
# per chromosome
pc = ue.groupby("CHR").FST.agg(["count", "mean", "median", "max"]).reset_index()
R["ue_chr"] = pc.round(5).to_dict("records")
# PBS
P = pb.copy()
st = dict(n=int(len(P)), mean=float(P.PBS_UZB.mean()), median=float(P.PBS_UZB.median()), max=float(P.PBS_UZB.max()), sd=float(P.PBS_UZB.std()),
          p95=float(P.PBS_UZB.quantile(.95)), p99=float(P.PBS_UZB.quantile(.99)), p999=float(P.PBS_UZB.quantile(.999)),
          t1=int(P.TIER1.sum()), t2=int(P.TIER2.sum()), t3=int(P.TIER3.sum()), np=int(P.NEAR_PRIVATE.sum()),
          gt01=int((P.PBS_UZB >= .1).sum()), gt015=int((P.PBS_UZB >= .15).sum()), gt03=int((P.PBS_UZB >= .3).sum()),
          maxdaf=float(P.DELTA_AF.max()), gt03daf=int((P.DELTA_AF >= .3).sum()))
R["pbs"] = st
tp = P.sort_values("PBS_UZB", ascending=False).head(15)
R["pbs_top15"] = [dict(snp=r.SNP, pbs=float(r.PBS_UZB), uzb=float(r.AF_UZB), eur=float(r.AF_EUR), eas=float(r.AF_EAS), sas=float(r.AF_SAS), afr=float(r.AF_AFR), daf=float(r.DELTA_AF), genes=gene_of(r.SNP)) for _, r in tp.iterrows()]
# old Tier-1 lookups
old = ["11:207698", "5:53879140", "3:133749168", "10:8512594", "11:20665570", "12:22967890", "12:125520190", "12:5664803"]
flag = set(open(D + "flagged_overlap.txt").read().split())
bimf = pd.read_csv(D + "merged_filt.bim", sep="\t", header=None)
inbim = set(bimf[1])
ol = []
for s in old:
    row = P[P.SNP == s]
    if len(row):
        r = row.iloc[0]
        ol.append(dict(snp=s, status="tested", pbs=float(r.PBS_UZB), uzb=float(r.AF_UZB), eur=float(r.AF_EUR), eas=float(r.AF_EAS), sas=float(r.AF_SAS), afr=float(r.AF_AFR)))
    else:
        ol.append(dict(snp=s, status="flagged (overlapping indel record)" if s in flag else "not in the NYGC-matched panel"))
R["old_tier1"] = ol
# correlation of per-SNP UZB-EUR FST with EUR-SAS etc. as sanity: UZB AF vs weighted mixture of refs
Pn = P.dropna()
R["af_corr"] = {k: float(np.corrcoef(Pn.AF_UZB, Pn[k])[0, 1]) for k in ["AF_EUR", "AF_EAS", "AF_SAS", "AF_AFR"]}
# best-fit mixture of EUR/EAS/SAS (non-negative least squares, AF space)
from scipy.optimize import nnls
A = np.vstack([Pn.AF_EUR, Pn.AF_EAS, Pn.AF_SAS, Pn.AF_AFR]).T; y = Pn.AF_UZB.values
c, rn = nnls(A, y); R["nnls_mix"] = dict(zip(["EUR", "EAS", "SAS", "AFR"], (c / c.sum()).round(3).tolist()), resid=float(rn / np.sqrt(len(y))))
json.dump(cache, open(cache_p, "w"))
json.dump(R, open(D + "results.json", "w"), indent=1)
print(json.dumps({k: R[k] for k in ["mds", "mds_var_pct", "ue_dist", "pbs", "old_tier1", "af_corr", "nnls_mix"]}, indent=1)[:5000])
