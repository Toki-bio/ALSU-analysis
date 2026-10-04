import pandas as pd, numpy as np, json
from scipy.stats import mannwhitneyu
D = "C:/work/alsu/nygc_site_2026-10-04/"
pairs = pd.read_csv(D + "pcrelate_vs_king_pihat_admix.tsv", sep="\t")
sm = pd.read_csv(D + "sample_map.tsv", sep="\t", dtype=str)
pf = pd.read_csv(D + "proj_sub.fam", sep=r"\s+", header=None, usecols=[0, 1], names=["FID", "IID"], dtype=str)
PQ = pd.read_csv(D + "proj_sub.4.Q", sep=r"\s+", header=None); PQ.columns = ["SAS", "EAS", "AFR", "EUR"]
qv = pd.concat([pf, PQ], axis=1).set_index("IID")[["EUR", "EAS", "SAS", "AFR"]]
n_have = int(qv.EUR.notna().sum())
# validation against the full NYGC run for samples present in both
nf = pd.read_csv("C:/work/alsu/figures_2026-10-04/adm_in.fam", sep=r"\s+", header=None, usecols=[0, 1], names=["FID", "IID"], dtype=str)
NQ = pd.read_csv("C:/work/alsu/figures_2026-10-04/adm_in.4.Q", sep=r"\s+", header=None); NQ.columns = ["SAS", "EAS", "AFR", "EUR"]
N = pd.concat([nf, NQ], axis=1); N = N[N.FID == "0"].drop_duplicates("IID")
mm = qv.reset_index().merge(sm, left_on="IID", right_on="safe_iid").merge(N, left_on="original_iid", right_on="IID", suffixes=("_p", "_n"))
val = {c: dict(r=float(np.corrcoef(mm[c + "_p"], mm[c + "_n"])[0, 1]), rmse=float(np.sqrt(((mm[c + "_p"] - mm[c + "_n"]) ** 2).mean()))) for c in ["EUR", "EAS", "SAS", "AFR"]}
a = qv.reindex(pairs.iid1).values; b = qv.reindex(pairs.iid2).values
pairs["q4_dist"] = np.sqrt(((a - b) ** 2).sum(1))
pairs["eas_abs"] = np.abs(a[:, 1] - b[:, 1])
pairs["eur_abs"] = np.abs(a[:, 0] - b[:, 0])
pairs["sas_abs"] = np.abs(a[:, 2] - b[:, 2])
have = pairs.q4_dist.notna()
out = dict(validation=val, n_validation=int(len(mm)), n_snps_projection=17484, n_samples_total=int(len(sm)), n_samples_with_q=n_have, n_pairs=int(len(pairs)), n_pairs_with_q=int(have.sum()))
conc = pairs[have & (pairs.status == "concordant_king_pihat")]; disc = pairs[have & (pairs.status == "discordant_pihat_only")]
out["conc"] = dict(n=int(len(conc)), q=float(conc.q4_dist.mean()), eas=float(conc.eas_abs.mean()), eur=float(conc.eur_abs.mean()), sas=float(conc.sas_abs.mean()))
out["disc"] = dict(n=int(len(disc)), q=float(disc.q4_dist.mean()), eas=float(disc.eas_abs.mean()), eur=float(disc.eur_abs.mean()), sas=float(disc.sas_abs.mean()))
out["mw_p"] = float(mannwhitneyu(disc.q4_dist, conc.q4_dist, alternative="greater").pvalue)
# tertiles over pairs with q
sub = pairs[have].copy()
sub["bin"] = pd.qcut(sub.q4_dist, 3, labels=["low", "mid", "high"])
tt = []
for k in ["low", "mid", "high"]:
    s = sub[sub.bin == k]
    tt.append(dict(bin=k, n=int(len(s)), conc=int((s.status == "concordant_king_pihat").sum()), disc=int((s.status == "discordant_pihat_only").sum()),
                   rate=float((s.status == "discordant_pihat_only").mean()), lo=float(s.q4_dist.min()), hi=float(s.q4_dist.max())))
out["tertiles"] = tt
un = pairs[~have]
out["unmatched"] = dict(n=int(len(un)), conc=int((un.status == "concordant_king_pihat").sum()), disc=int((un.status == "discordant_pihat_only").sum()), rate=float((un.status == "discordant_pihat_only").mean()))
# by pair type coverage
out["by_type"] = {t: dict(n=int((pairs.pair_type == t).sum()), with_q=int(((pairs.pair_type == t) & have).sum())) for t in pairs.pair_type.unique()}
# the same test on the pairs that also had old K=5 distances (so the two are comparable)
old = pairs.q5_dist.notna() & have
c2 = pairs[old & (pairs.status == "concordant_king_pihat")]; d2 = pairs[old & (pairs.status == "discordant_pihat_only")]
out["both"] = dict(n=int(old.sum()), conc_new=float(c2.q4_dist.mean()), disc_new=float(d2.q4_dist.mean()), conc_old=float(c2.q5_dist.mean()), disc_old=float(d2.q5_dist.mean()),
                   r=float(np.corrcoef(pairs[old].q4_dist, pairs[old].q5_dist)[0, 1]))
# control: scramble the Q vectors across samples (distance should no longer separate groups)
rng = np.random.default_rng(3)
sc = []
vals = qv.dropna().values
ids = qv.dropna().index.tolist()
pos = {i: k for k, i in enumerate(ids)}
for _ in range(2000):
    perm = rng.permutation(len(ids)); v2 = vals[perm]
    ok = pairs.iid1.isin(pos) & pairs.iid2.isin(pos)
    aa = v2[pairs[ok].iid1.map(pos).values]; bb = v2[pairs[ok].iid2.map(pos).values]
    dist = np.sqrt(((aa - bb) ** 2).sum(1)); st = pairs[ok].status.values
    sc.append(float(dist[st == "discordant_pihat_only"].mean() - dist[st == "concordant_king_pihat"].mean()))
obs = out["disc"]["q"] - out["conc"]["q"]
out["scramble_diff"] = dict(mean=float(np.mean(sc)), sd=float(np.std(sc)), observed=obs, p_perm=float((np.sum(np.array(sc) >= obs) + 1) / (len(sc) + 1)), n_perm=len(sc))
cnt = pd.concat([pairs.iid1, pairs.iid2]).value_counts()
out["reuse"] = dict(n_unique_samples=int(len(cnt)), max_pairs_per_sample=int(cnt.max()), top10_share=float(cnt.head(10).sum() / cnt.sum()), n_samples_in_100plus=int((cnt >= 100).sum()))
json.dump(out, open(D + "rel_q_results.json", "w"), indent=1)
R = json.load(open(D + "results.json")); R["rel_q"] = out; json.dump(R, open(D + "results.json", "w"), indent=1)
print(json.dumps(out, indent=1))
