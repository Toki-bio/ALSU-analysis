from sg import *
import pandas as pd, numpy as np

t = read("steps/admixture_aware_relatedness.html")
RQ = R["rel_q"]; V = RQ["validation"]; TT = RQ["tertiles"]; SD = RQ["scramble_diff"]; RU = RQ["reuse"]
pairs = pd.read_csv(D + "pcrelate_vs_king_pihat_admix.tsv", sep="\t")
im = pd.read_csv(D + "hub_q.imiss", sep=r"\s+"); he = pd.read_csv(D + "hub_q.het", sep=r"\s+")
qc = im.merge(he, on=["FID", "IID"])
cnt = pd.concat([pairs.iid1, pairs.iid2]).value_counts()
hubs = list(cnt[cnt >= 100].index)
hub_rows = []
for h in hubs:
    s = pairs[(pairs.iid1 == h) | (pairs.iid2 == h)]
    r = qc[qc.IID == h].iloc[0]
    hub_rows.append([h, len(s), f"{s.pihat.mean():.3f}", f"{100 * (s.status == 'discordant_pihat_only').mean():.0f}%", f"{100 * r.F_MISS:.1f}%", f"{r.F:+.2f}", "no"])
inv = pairs.iid1.isin(hubs) | pairs.iid2.isin(hubs)
nonhub = pairs[~inv]
f_med = qc.F.median(); m_med = qc.F_MISS.median(); m_mean = qc.F_MISS.mean()
n_hub5 = int((cnt >= 1000).sum())

val_rows = [[c, f"{V[c]['r']:.3f}", f"{V[c]['rmse']:.3f}"] for c in ("EUR", "EAS", "SAS", "AFR")]
ter_rows = [[b["bin"], f"{b['lo']:.2f}&ndash;{b['hi']:.2f}", f"{b['n']:,}", f"{b['conc']:,}", f"{b['disc']:,}", f"{100 * b['rate']:.0f}%"] for b in TT]

blk = f"""
<h2 style="border-bottom:3px solid #16a34a;color:#1b5e20">UPDATE 2026-10-04: what the PI_HAT&ndash;KING discordance is made of</h2>
{box("warn", "<strong>Two corrections to the analysis below.</strong> (1) The ancestry distances were taken from the deprecated global K=5 run; they were recomputed on the NYGC-reference K=4 ancestry for all 1,193 samples. (2) Doing so exposed that the discordance is dominated by a handful of low-quality samples, not by ancestry. The sections below are kept for the record; where they conflict with this update, this update applies.")}

<h3>1. Five samples explain almost all flagged pairs</h3>
<p>Of the {len(pairs):,} pairs with PI_HAT &ge; 0.185, {int(inv.sum()):,} ({100 * inv.mean():.1f}%) contain at least one of {len(hubs)} samples, and {n_hub5} of them (ALSU_0242, ALSU_0402, ALSU_0499, ALSU_1030, GWAS2026_0016) appear in <strong>{RU['max_pairs_per_sample']:,} pairs each</strong>, i.e. they exceed PI_HAT 0.185 with nearly every other sample in the merged data set ({RU['n_unique_samples']:,} samples are involved in total). A person is not a close relative of 1,184 others; this is a data-quality signature. Only {len(nonhub)} pairs involve none of these samples.</p>
{table(["Sample (alias)", "Pairs flagged", "Mean PI_HAT in those pairs", "Share discordant (PI_HAT-only)", "Missing genotypes", "Inbreeding coefficient F", "In final 1,256 cohort"], hub_rows)}
<p>For comparison, over all {len(qc):,} samples the median missingness is {100 * m_med:.2f}% (mean {100 * m_mean:.1f}%) and the median F is {f_med:+.2f}. The six samples have 11&ndash;18% missing genotypes and a strongly negative F, i.e. far more heterozygous genotypes than expected. High missingness together with excess heterozygosity is typical of poor-quality or mixed DNA, and it inflates PLINK&rsquo;s PI_HAT against everyone. This is consistent with, but does not prove, the cause; the original array QC used a call-rate cut-off of 20% missing, which these samples pass. (QC measured on the 104,411 LD-pruned SNPs of the merged set; script <code>hub.sh</code>, DRAGEN <code>/staging/tmp/scratch/hub_check/</code>.)</p>
{box("key", "<strong>None of the six samples is in the final cohort of 1,256</strong> (the final QC of step 7 applies a 5% missingness limit, which samples with 11&ndash;18% missing genotypes would not pass), and only 4 of the 6,228 flagged pairs have both members in the final cohort. The discordance study therefore describes a pre-QC array data set. For the final cohort, relatedness is minimal (11 pairs at second degree or closer; <a href='step15.html'>step 15</a>).")}
<p><strong>Consequence for the pair supported by both PC-Relate and RelateAdmix:</strong> ALSU_0218 / ALSU_0242 (PC-Relate 0.153, RelateAdmix 0.233) involves ALSU_0242, one of the five poor-quality samples. It should not be treated as a trustworthy relative pair without re-genotyping or a check of the sample.</p>

<h3>2. Ancestry distances recomputed on the NYGC reference (K=4)</h3>
<p>The 1,193 samples of the relatedness set were projected onto the fixed allele frequencies of the NYGC K=4 ADMIXTURE run (<code>admixture -P</code> with <code>adm_in.4.P</code>) on the 17,484 SNPs present in both data sets with identical alleles; script <code>/staging/tmp/proj_k4.sh</code>, output <code>/staging/tmp/scratch/proj_k4/</code>. This gives ancestry for all 1,193 samples; the earlier K=5 distances existed for only 1,181 of the 6,228 pairs. Validation: for the {RQ['n_validation']:,} samples that are also in the full K=4 run the projected and the full-run components agree as follows.</p>
{table(["Component", "Correlation", "RMSE"], val_rows)}
<p>Mean Euclidean K=4 distance between the two members of a pair: <strong>{RQ['conc']['q']:.3f}</strong> for KING-and-PI_HAT concordant pairs (n = {RQ['conc']['n']:,}) and <strong>{RQ['disc']['q']:.3f}</strong> for PI_HAT-only discordant pairs (n = {RQ['disc']['n']:,}). The direction matches the earlier K=5 result (0.130 vs 0.160 on the 1,181 pairs that had it; on those same pairs the new distances are 0.140 vs 0.165, and the two distance sets correlate at r = {RQ['both']['r']:.2f}). By tertile of distance:</p>
{table(["Distance tertile", "K=4 distance", "Pairs", "Concordant", "Discordant", "Discordance rate"], ter_rows)}
{box("info", f"<strong>The ancestry effect is much weaker than the earlier text claimed.</strong> The naive Mann&ndash;Whitney test (p = {RQ['mw_p']:.0e}) treats 6,228 pairs as independent, but they are built from {RU['n_unique_samples']:,} samples of which ten account for {100 * RU['top10_share']:.0f}% of all pair slots. A permutation test that reassigns whole ancestry vectors among samples gives a difference of {SD['observed']:.3f} against a null of {SD['mean']:.3f} &plusmn; {SD['sd']:.3f} (one-sided p = {SD['p_perm']:.3f}, {SD['n_perm']:,} permutations): borderline at best. The discordance rate is not monotone across tertiles ({100 * TT[0]['rate']:.0f}%, {100 * TT[1]['rate']:.0f}%, {100 * TT[2]['rate']:.0f}%; the earlier K=5 table read 83%, 90%, 93%). The statement &ldquo;direct K=5 distance did explain discordance&rdquo; is therefore withdrawn; the poor-quality samples are a far stronger explanation.")}

<h3>3. RelateAdmix and PC-Relate do not need re-running on the NYGC reference</h3>
<p>RelateAdmix received its ancestry input from an ADMIXTURE K=5 run on the merged 1,193-sample LD-pruned file itself (<code>merged_pruned.5.Q/.P</code>; script <code>scripts/remote_relateadmix_run.sh</code>), and PC-Relate uses PCs of the same file. Neither used the deprecated global 1000G run, so their results (27 pairs with RelateAdmix kinship-equivalent &ge; 0.0884; 48 with PC-Relate &ge; 0.0884; 1 resp. 5 of the PI_HAT-only pairs) are unaffected by the reference change. Only the &ldquo;Step 11 K=5 distance&rdquo; analysis and the birthplace map depended on the deprecated run, and both are replaced here.</p>

<h3>4. Revised recommendation</h3>
<ul>
<li>Do not use PI_HAT &ge; 0.185 lists from the merged pre-QC array data to remove samples. Check per-sample missingness and heterozygosity first: a sample with PI_HAT above threshold against hundreds of others is a QC failure.</li>
<li>Re-check the QC cut-offs of the original array pipeline (20% missingness) in the light of these samples; the final cohort already excludes them.</li>
<li>Keep KING as the operational screen and PC-Relate/RelateAdmix as validation, on samples that passed QC.</li>
</ul>
<h3>5. Birthplace map on the NYGC reference</h3>
{img("alsu_uzbekistan_k4_nygc_birthplace_map.png", "Mean K=4 ADMIXTURE proportions (expanded cohort + NYGC 30x reference) by Uzbekistan birthplace region; pies sized by &radic;N; 1,021 mapped samples. Component names come from the reference means (EUR 0.97, EAS 1.00, SAS 0.89, AFR 0.96). Values: data/alsu_uzbekistan_k4_nygc_birthplace_summary.tsv; script scripts/build_uzbekistan_admixture_map_nygc.py.")}
<p>The regional pattern is flat: mean EUR-like shares range from 0.44 to 0.54 and EAS-like shares from 0.20 to 0.38 across regions, with small-sample regions (n &lt; 30) noisy. The statistical test is on <a href="step11.html">step 11</a>.</p>
"""
# replace the old banner, mark old sections
t = banner(t, "<strong>Updated 2026-10-04:</strong> the ancestry distances were recomputed on the NYGC-reference K=4 ancestry and the main finding was revised (the PI_HAT-only discordance is dominated by five poor-quality samples that are not in the final cohort). See the update below the summary cards; text further down is kept for the record. <a href=\"deprecated_old_runs.html\">List of current and deprecated runs</a>.")
t = put(t, "relate", blk, before='<div class="figure">')
# old figure -> keep in a collapsed block
i = t.index('<div class="figure">');
# find the end of the first figure div (the one that directly follows our block)
i = t.index('<div class="figure">', t.index("<!--NYGC:relate:END-->"))
j = t.index("</div>\n</div>", i) + len("</div>\n</div>") if "</div>\n</div>" in t[i:i + 3000] else None
if j is None:
    j = t.index("<details", i)
old_fig = t[i:j]
if "NYGC:relate:OLDFIG" not in t:
    t = t[:i] + '<!--NYGC:relate:OLDFIG--><details style="margin:14px 0"><summary style="cursor:pointer"><strong>Superseded: August 2026 map on the deprecated global K=5 run (kept for the record)</strong></summary>' + old_fig + "</details>" + t[j:]
a = '<h3 id="direct-step-11-k-5-distance-did-explain-discordance">Direct Step 11 K=5 distance did explain discordance</h3>'
t = t.replace(a, '<h3 id="direct-step-11-k-5-distance-did-explain-discordance">Direct Step 11 K=5 distance did explain discordance <span style="color:#991b1b;font-size:.8em">[WITHDRAWN 2026-10-04: deprecated K=5 run; effect weak and confounded by five poor-quality samples, see the update at the top]</span></h3>')
b = '<h2 id="draft-discussion">Draft Discussion</h2>'
t = t.replace(b, b + '<div class="callout" style="background:#fff4e5;border-left:4px solid #d97706"><strong>Correction 2026-10-04:</strong> statements below that ancestry (Step 11 K=5 distance) explains the PI_HAT&ndash;KING discordance are withdrawn; see the update at the top of the page. The discordance is dominated by five poor-quality samples absent from the final cohort.</div>', 1)
write("steps/admixture_aware_relatedness.html", t)
print("relate ok")
