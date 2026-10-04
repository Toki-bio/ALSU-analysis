import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D = "C:/work/alsu/nygc_site_2026-10-04/"; I = "C:/work/ALSU-analysis/images/"
L = json.load(open(D + "ld_bins.json"))
N = {"UZB": 1256, "EUR": 633, "EAS": 585, "SAS": 601, "AFR": 893}
col = {"UZB": "#2e7d32", "EUR": "#7570b3", "EAS": "#1b9e77", "SAS": "#e7298a", "AFR": "#d95f02"}
mid = [(b[0] + b[1]) / 2 for b in L["UZB"]]
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), dpi=130, sharex=True)
for k, rows in L.items():
    raw = [r[2] for r in rows]
    adj = [r[2] - 1 / N[k] for r in rows]
    axs[0].plot(mid, raw, "o-", c=col[k], label=k, lw=2 if k == "UZB" else 1.2)
    axs[1].plot(mid, adj, "o-", c=col[k], label=k, lw=2 if k == "UZB" else 1.2)
for a, t in zip(axs, ["mean r\u00b2 (raw)", "mean r\u00b2 minus 1/n (sample-size baseline removed)"]):
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("distance between SNPs (kb, bin midpoint)"); a.set_ylabel(t, fontsize=8); a.grid(alpha=.3, which="both")
axs[0].legend(fontsize=8)
fig.suptitle("LD decay on the same 30,000 random SNPs (MAF \u2265 0.05 in each group), NYGC 30x and Uzbek imputed data", fontsize=9)
fig.tight_layout(); fig.savefig(I + "nygc_ld_decay.png"); plt.close()
R = json.load(open(D + "results.json")); R["ld"] = L; R["ld_n"] = N
json.dump(R, open(D + "results.json", "w"), indent=1)
print("ok")
