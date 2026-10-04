import pandas as pd, numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = "C:/work/alsu/nygc_site_2026-10-04/"; I = "C:/work/ALSU-analysis/images/"
R = json.load(open(D + "results.json"))
pb = pd.read_csv(D + "pbs_all.tsv", sep="\t")
ue = pd.read_csv(D + "filt_EUR_UZB.fst", sep="\t").dropna(subset=["FST"]).merge(pb[["SNP","AF_UZB","AF_EUR","AF_EAS","AF_SAS","AF_AFR"]], on="SNP")
w = R["nnls_mix"]
pred = lambda d: w["EUR"]*d.AF_EUR + w["EAS"]*d.AF_EAS + w["SAS"]*d.AF_SAS + w["AFR"]*d.AF_AFR
ue["pred"] = pred(ue); pb["pred"] = pred(pb)
# mixture test
res = (ue.AF_UZB - ue.pred)
top = ue.sort_values("FST", ascending=False).head(30)
between = ((top.AF_UZB >= np.minimum(top.AF_EUR, top.AF_EAS)) & (top.AF_UZB <= np.maximum(top.AF_EUR, top.AF_EAS))).sum()
R["mix_test"] = dict(rmse_all=float(np.sqrt((res**2).mean())), r_all=float(np.corrcoef(ue.AF_UZB, ue.pred)[0,1]),
  top30_between_eur_eas=int(between), top30_mean_abs_resid=float((top.AF_UZB-top.pred).abs().mean()), top30_max_abs_resid=float((top.AF_UZB-top.pred).abs().max()),
  frac_abs_resid_gt_0p1=float((res.abs()>.1).mean()), n_abs_resid_gt_0p15=int((res.abs()>.15).sum()),
  max_abs_resid=float(res.abs().max()))
# the same for PBS top 15
tp = pb.sort_values("PBS_UZB", ascending=False).head(15)
R["mix_test"]["pbs_top15_mean_abs_resid"] = float((tp.AF_UZB - tp.pred).abs().mean())
# SNPs with UZB outside [min,max] of the four references by >0.05
lo = pb[["AF_EUR","AF_EAS","AF_SAS","AF_AFR"]].min(axis=1); hi = pb[["AF_EUR","AF_EAS","AF_SAS","AF_AFR"]].max(axis=1)
out = pb[(pb.AF_UZB < lo - .05) | (pb.AF_UZB > hi + .05)]
R["mix_test"]["uzb_outside_ref_range_gt0p05"] = int(len(out))
json.dump(R, open(D + "results.json", "w"), indent=1)
print(R["mix_test"])
col = {"UZB":"#2e7d32"}
# 1 FST histogram UZB-EUR
fig, ax = plt.subplots(figsize=(7,3.6), dpi=130)
ax.hist(ue.FST.clip(lower=0), bins=80, color="#1565c0"); ax.set_yscale("log")
ax.axvline(R["ue_dist"]["mean"], color="#c62828", ls="--", lw=1); ax.text(R["ue_dist"]["mean"]+.003, ax.get_ylim()[1]*.3, f"mean {R['ue_dist']['mean']:.4f}", color="#c62828", fontsize=8)
ax.set_xlabel("per-SNP Weir-Cockerham FST, UZB vs EUR (negative values set to 0)"); ax.set_ylabel("SNPs (log scale)")
ax.set_title(f"UZB vs EUR, {len(ue):,} SNPs (NYGC 30x reference)", fontsize=10); fig.tight_layout(); fig.savefig(I + "nygc_fst_uzb_eur_hist.png"); plt.close()
# 2 Manhattan-style scatter helper
def manh(df, ycol, ylabel, title, fn, thr=None, label=None):
    d = df.copy(); d["CHR"] = d.SNP.str.split(":").str[0].astype(int); d["BP"] = d.SNP.str.split(":").str[1].astype(int)
    d = d.sort_values(["CHR","BP"]); off = 0; xs=[]; tick=[]; lab=[]
    for c, g in d.groupby("CHR"):
        xs.append(g.BP + off); tick.append(off + g.BP.max()/2); lab.append(c); off += g.BP.max() + 5e6
    d["x"] = pd.concat(xs)
    fig, ax = plt.subplots(figsize=(10,3.4), dpi=130)
    for i,(c,g) in enumerate(d.groupby("CHR")): ax.scatter(g.x, g[ycol], s=2, c="#1565c0" if i%2==0 else "#78909c")
    if thr: ax.axhline(thr, color="#c62828", ls="--", lw=1); ax.text(ax.get_xlim()[0], thr, f" {label}", color="#c62828", fontsize=8, va="bottom")
    ax.set_xticks(tick); ax.set_xticklabels(lab, fontsize=7); ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10); fig.tight_layout(); fig.savefig(I + fn); plt.close()
manh(ue.rename(columns={"SNP":"SNP"}), "FST", "FST (UZB vs EUR)", "Per-SNP FST, UZB vs EUR (NYGC reference)", "nygc_fst_uzb_eur_manhattan.png")
manh(pb, "PBS_UZB", "PBS (UZB branch)", "Per-SNP PBS, Uzbek branch of UZB-EUR-EAS (NYGC reference); Tier 1 threshold 0.3 lies far above the data", "nygc_pbs_manhattan.png", thr=.3, label="Tier 1 threshold 0.3")
# 3 mixture scatter
fig, axs = plt.subplots(1,2, figsize=(9,4), dpi=130)
axs[0].hexbin(ue.pred, ue.AF_UZB, gridsize=60, bins="log", cmap="viridis", mincnt=1); axs[0].plot([0,1],[0,1],"r--",lw=.8)
axs[0].set_xlabel("predicted Uzbek allele freq. (non-negative fit to EUR/EAS/SAS/AFR)"); axs[0].set_ylabel("observed Uzbek allele freq."); axs[0].set_title(f"All {len(ue):,} SNPs, r = {R['mix_test']['r_all']:.3f}", fontsize=9)
axs[0].xaxis.label.set_size(7)
t = top; axs[1].scatter(t.pred, t.AF_UZB, c="#c62828", s=16); axs[1].plot([0,1],[0,1],"k--",lw=.8)
axs[1].set_xlim(0,1); axs[1].set_ylim(0,1); axs[1].set_xlabel("predicted"); axs[1].set_ylabel("observed"); axs[1].set_title("30 highest-FST SNPs (UZB vs EUR)", fontsize=9)
fig.tight_layout(); fig.savefig(I + "nygc_mixture_fit.png"); plt.close()
# 4 FST heatmap & MDS (static) for step14
pops = ["UZB","SAS","EUR","EAS","AFR"]; M = np.array([[0 if a==b else R["matrix"][a][b] for b in pops] for a in pops])
fig, axs = plt.subplots(1,2, figsize=(10,4.2), dpi=130)
im = axs[0].imshow(M, cmap="YlOrRd"); axs[0].set_xticks(range(5)); axs[0].set_xticklabels(pops); axs[0].set_yticks(range(5)); axs[0].set_yticklabels(pops)
for i in range(5):
    for j in range(5):
        if i!=j: axs[0].text(j,i,f"{M[i,j]:.4f}",ha="center",va="center",fontsize=8)
axs[0].set_title("Weighted FST (NYGC)", fontsize=10)
cc = {"UZB":"#2e7d32","SAS":"#e7298a","EUR":"#7570b3","EAS":"#1b9e77","AFR":"#d95f02"}
for p,(x,y) in R["mds"].items(): axs[1].scatter(x,y,s=70,c=cc[p]); axs[1].annotate(p,(x,y),textcoords="offset points",xytext=(6,6),fontsize=9)
axs[1].set_xlabel(f"MDS 1 ({R['mds_var_pct'][0]:.1f}% of positive eigenvalues)"); axs[1].set_ylabel(f"MDS 2 ({R['mds_var_pct'][1]:.1f}%)"); axs[1].set_title("Classical MDS of the FST matrix", fontsize=10); axs[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(I + "nygc_fst_heatmap_mds.png"); plt.close()
# slide MDS (deck) same map
fig, ax = plt.subplots(figsize=(5.6,4), dpi=200)
for p,(x,y) in R["mds"].items(): ax.scatter(x,y,s=90,c=cc[p],zorder=3); ax.annotate(p,(x,y),textcoords="offset points",xytext=(7,7),fontsize=10,weight="bold")
ax.set_xlabel("MDS dimension 1"); ax.set_ylabel("MDS dimension 2"); ax.set_title("Classical MDS of the FST distance matrix (NYGC)", fontsize=10); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig("C:/work/alsu/figures_2026-10-04/fst_mds.png"); plt.close()
