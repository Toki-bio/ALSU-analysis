from sg import *

t = read("steps/step9.html")
P = R["pairs"]; U = R["ue_dist"]; MT = R["mix_test"]
order = ["SAS_UZB", "EUR_UZB", "EAS_UZB", "AFR_UZB", "EUR_SAS", "EAS_SAS", "EAS_EUR", "AFR_SAS", "AFR_EUR", "AFR_EAS"]


def nm(k):
    a, b = k.split("_")
    return f"UZB vs {a}" if b == "UZB" else f"{a} vs {b}" if "UZB" not in k else f"UZB vs {b}"


def nm2(k):
    a, b = k.split("_")
    if "UZB" in k:
        return "UZB vs " + (a if b == "UZB" else b)
    return f"{b} vs {a}" if k in ("EAS_EUR",) else f"{a} vs {b}"


rows = [[nm2(k), f(P[k]["weighted_all"]), f(P[k]["weighted_filt"]),
         f"{P[k]['mean_filt']:.4f} &plusmn; {P[k]['mean_se']:.4f}", f"{P[k]['n']:,}"] for k in order]


def gene(g):
    named = [x for x in (g or []) if not x.startswith("ENSG")]
    if named:
        return ", ".join(named)
    return "Ensembl gene without symbol" if g else "no gene overlap"


top = [[i + 1, r["snp"], f(r["fst"], 3), f(r["uzb"], 2), f(r["eur"], 2), f(r["eas"], 2), f(r["sas"], 2), f(r["afr"], 2),
        esc(gene(r["genes"]))] for i, r in enumerate(R["ue_top30"])]
chr_rows = [[c["CHR"], f"{c['count']:,}", f(c["mean"], 4), f(c["median"], 4), f(c["max"], 3)] for c in R["ue_chr"]]
hist = [[esc(h[0]), f"{h[1]:,}", f"{100 * h[1] / U['n']:.1f}%"] for h in R["ue_hist"]]
se_diff = (P['EUR_UZB']['mean_filt'] - P['SAS_UZB']['mean_filt']) / (P['EUR_UZB']['mean_se'] ** 2 + P['SAS_UZB']['mean_se'] ** 2) ** .5
uz = [r['uzb'] for r in R['ue_top30']]
chm = [c['mean'] for c in R['ue_chr']]

blk = f"""
<h2>1. Overview</h2>
<p><strong>F<sub>ST</sub></strong> measures how much allele frequencies differ between two populations: 0 means identical frequencies, 1 means fixed different alleles. This step estimates it for the Uzbek cohort against four continental reference groups and uses the per-SNP values to ask two questions: how far is the Uzbek cohort from each reference, and does any single locus differ much more than the genome-wide level would predict.</p>
{box("warn", "<strong>This page was rebuilt on 2026-10-04.</strong> All numbers below come from the NYGC 30x high-coverage GRCh38 reference (md5-verified) and the expanded cohort of 1,256 Uzbek samples. The earlier hg19 pipeline (2025), the August 2026 run on the old GRCh38 reference file, and the &ldquo;256 loci with F<sub>ST</sub> &gt; 0.5&rdquo; are withdrawn; the hg19 pipeline is kept collapsed at the bottom for provenance only.")}
{cards([("1,256", "Uzbek samples"), ("2,712", "NYGC reference (EUR 633, EAS 585, SAS 601, AFR 893)"), ("82,744", "SNPs, all ten pairs"), (f"{U['n']:,}", "SNPs after indel filter"), (f(P['EUR_UZB']['weighted_all']), "UZB vs EUR weighted F<sub>ST</sub>"), (f(U['max'], 3), "highest single-SNP F<sub>ST</sub>")])}
<h3>How to read an F<sub>ST</sub> value</h3>
{table(["F<sub>ST</sub>", "Reading", "This analysis"], [["0.00&ndash;0.05", "little differentiation (same continental group)", "UZB vs SAS, UZB vs EUR, UZB vs EAS (0.014&ndash;0.040); EUR vs SAS (0.031)"], ["0.05&ndash;0.15", "moderate", "EAS vs SAS (0.057), EUR vs EAS (0.087), UZB vs AFR and SAS vs AFR (about 0.13)"], ["0.15&ndash;0.25", "great", "EAS vs AFR (0.167)"], ["&gt;0.25", "very great", "no pair; no single SNP reaches 0.25"]])}

<h2>2. Data and method</h2>
<ul>
<li><strong>Cohort:</strong> expanded ALSU Uzbek cohort, 1,256 samples (array, imputed; see <a href="cohorts_and_sample_sets.html">cohorts and sample sets</a>).</li>
<li><strong>Reference:</strong> 1000 Genomes high-coverage (NYGC) GRCh38, 3,202 samples, all 22 chromosome files verified against the official md5 manifest. Superpopulations used: EUR, EAS, SAS, AFR. AMR is not used. Reference samples were not filtered for relatedness.</li>
<li><strong>Markers:</strong> 82,744 LD-pruned SNPs present in both data sets, merged by position and alleles (<code>/staging/ALSU-analysis/spring2026/full_expanded_cohort/nygc30x_1256/merged_nygc.bim</code> on DRAGEN). A second, stricter set drops the 8,573 sites where the NYGC call set has an overlapping indel, multiallelic or structural-variant record within 1 bp, because single-base comparisons are not reliable there (list: <code>nygc30x_1256_filt/flagged_overlap.txt</code>). This leaves {U['n']:,} sites for the UZB pairs. The number of SNPs per pair in the table can be smaller because PLINK drops sites that are monomorphic in both groups.</li>
<li><strong>Estimator:</strong> PLINK 1.9 <code>--fst --within</code>, Weir &amp; Cockerham. &ldquo;Weighted&rdquo; is the ratio of summed numerators to summed denominators over all SNPs (the standard genome-wide value); &ldquo;mean&rdquo; is the plain average of per-SNP values. The &plusmn; is a leave-one-chromosome-out jackknife standard error of the mean per-SNP value (22 blocks), so it reflects variation between chromosomes.</li>
</ul>

<h2>3. Pairwise F<sub>ST</sub>, all ten pairs</h2>
{table(["Pair", "Weighted F<sub>ST</sub>, 82,744 SNPs", "Weighted F<sub>ST</sub>, indel-filtered", "Mean per-SNP F<sub>ST</sub> &plusmn; jackknife SE", "SNPs (filtered)"], rows)}
{box("key", f"<strong>Result.</strong> The Uzbek cohort is about equally close to South Asians ({f(P['SAS_UZB']['weighted_all'])}) and Europeans ({f(P['EUR_UZB']['weighted_all'])}), clearly farther from East Asians ({f(P['EAS_UZB']['weighted_all'])}) and farthest from Africans ({f(P['AFR_UZB']['weighted_all'])}). The ranking is the same with or without the 8,573 flagged sites; no weighted value moves by more than 0.0005. The UZB&ndash;AFR distance ({f(P['AFR_UZB']['weighted_all'])}) is almost identical to SAS&ndash;AFR ({f(P['AFR_SAS']['weighted_all'])}). The matrix, heatmap and MDS are on <a href='step14.html'>step 14</a>.")}
<p>The UZB&ndash;SAS and UZB&ndash;EUR values differ by only {f(P['EUR_UZB']['weighted_all'] - P['SAS_UZB']['weighted_all'])}, which is about {se_diff:.0f} jackknife standard errors of the mean per-SNP value, so the difference is measurable, but both distances are small. &ldquo;Closer to SAS or to EUR&rdquo; should not be over-read: it depends on the reference panel, which has no Central Asian populations.</p>

<h2>4. Distribution of per-SNP F<sub>ST</sub>, UZB vs EUR</h2>
{img("nygc_fst_uzb_eur_hist.png", f"Per-SNP F<sub>ST</sub>, UZB vs EUR, {U['n']:,} SNPs after the indel filter. Log scale on the y-axis; the dashed line is the mean.")}
{table(["F<sub>ST</sub> bin", "SNPs", "Share"], hist)}
<p>Median {f(U['median'], 4)}, mean {f(U['mean'], 4)}, 95th percentile {f(U['q']['0.95'], 3)}, 99th {f(U['q']['0.99'], 3)}, 99.9th {f(U['q']['0.999'], 3)}, maximum {f(U['max'], 3)}. {U['gt01']:,} SNPs exceed 0.10, {U['gt02']} exceed 0.20 and none exceeds 0.30. A tail towards 0.5&ndash;1 would indicate fixed or near-fixed differences; here the tail ends at 0.24.</p>
{img("nygc_fst_uzb_eur_manhattan.png", "Per-SNP F<sub>ST</sub> along the genome, UZB vs EUR. No chromosome or region stands out.")}

<h2>5. The 30 most differentiated SNPs, and why they are not selection candidates</h2>
{table(["#", "SNP (GRCh38)", "F<sub>ST</sub>", "UZB", "EUR", "EAS", "SAS", "AFR", "Gene overlap (Ensembl)"], top)}
<p>Frequencies are for the allele that is the Uzbek minor allele. Gene overlap was looked up by position in Ensembl on 2026-10-04; &ldquo;no gene overlap&rdquo; means no annotated gene covers the position.</p>
{box("info", f"<strong>What the pattern says.</strong> In all 30 SNPs the Uzbek frequency lies between the EUR and EAS frequencies, near 0.5 ({min(uz):.2f}&ndash;{max(uz):.2f}), while EUR and EAS sit at opposite ends. That is what a population made of roughly equal parts of two differentiated sources looks like, and it is why these SNPs rank highest for UZB vs EUR. As a direct test, a non-negative least-squares fit of the Uzbek frequency as a mixture of the four references gives EUR {R['nnls_mix']['EUR'] * 100:.0f}%, EAS {R['nnls_mix']['EAS'] * 100:.0f}%, SAS {R['nnls_mix']['SAS'] * 100:.0f}%, AFR {R['nnls_mix']['AFR'] * 100:.0f}%, close to the ADMIXTURE K=4 result (about 50% / 30% / 19% / 0%, <a href='step11.html'>step 11</a>), and predicts the observed Uzbek frequency at all {U['n']:,} SNPs with r = {MT['r_all']:.3f} (RMSE {MT['rmse_all']:.3f}). The 30 top SNPs are predicted to within {MT['top30_max_abs_resid']:.2f} at most (mean {MT['top30_mean_abs_resid']:.3f}). Only {MT['frac_abs_resid_gt_0p1'] * 100:.2f}% of all SNPs deviate from the mixture prediction by more than 0.10, and the largest deviation is {MT['max_abs_resid']:.2f}.")}
{img("nygc_mixture_fit.png", "Observed Uzbek allele frequency against the frequency predicted from a fitted mixture of the references. Alleles are oriented to the Uzbek minor allele, so observed values do not exceed 0.5.")}
<p>The list contains known pigmentation loci (for example <em>OCA2</em>, <em>BNC2</em>), which differ strongly between Europeans and East Asians and are therefore expected near the top of the F<sub>ST</sub> ranking of any population with ancestry from both. Nothing in this list needs a selection explanation, and none was tested for one.</p>

<h2>6. F<sub>ST</sub> by chromosome (UZB vs EUR)</h2>
{table(["Chr", "SNPs", "Mean", "Median", "Max"], chr_rows)}
<p>Chromosome means range from {min(chm):.4f} to {max(chm):.4f}, around the genome-wide mean of {f(U['mean'], 4)}; no chromosome is an outlier.</p>

<h2>7. What changed relative to earlier versions of this page</h2>
{table(["Earlier statement", "Status now"], [
    ["UZB vs EUR weighted F<sub>ST</sub> 0.0150; UZB vs EAS 0.0396; UZB vs SAS 0.0145; UZB vs AFR 0.1283 (August 2026, old reference)", f"Replicated on NYGC: {f(P['EUR_UZB']['weighted_all'])}, {f(P['EAS_UZB']['weighted_all'])}, {f(P['SAS_UZB']['weighted_all'])}, {f(P['AFR_UZB']['weighted_all'])}. Same ranking, differences up to 0.0027."],
    ["256 loci with F<sub>ST</sub> &gt; 0.5; loci with F<sub>ST</sub> = 0.998 (chr9 14.8 Mb, chr1 40.9 Mb) (hg19 run, 2025)", "Not reproducible: the highest per-SNP F<sub>ST</sub> on NYGC is " + f(U['max'], 3) + ". The earlier values came from the old reference file and the 2025 pipeline; do not cite."],
    ["&ldquo;Annotate the 256 extreme loci&rdquo;, &ldquo;highly differentiated loci: potential biology&rdquo;", "Withdrawn. Replaced by section 5 (the top SNPs are admixture-informative, not selection signals)."],
    ["Reference sizes EUR 503 / EAS 504 / SAS 489 / AFR 660 (1000G phase 3)", "NYGC high-coverage sizes: EUR 633 / EAS 585 / SAS 601 / AFR 893."]])}

<h2>8. Files</h2>
<ul>
<li>DRAGEN, all pairs and logs: <code>/staging/ALSU-analysis/spring2026/full_expanded_cohort/pbs_refqc_2026-10/fst_all_pairs_2026-10-04/</code> (<code>all_*.fst</code> for 82,744 SNPs, <code>filt_*.fst</code> for the indel-filtered set); script <code>/staging/tmp/scratch/alsu_fst_all.sh</code>.</li>
<li>Tables and figures on this page were produced from those files by <code>analysis.py</code> and <code>figs.py</code> (project folder <code>C:\\work\\alsu\\nygc_site_2026-10-04\\</code>, copied into the repository under <code>scripts/nygc_site_2026-10-04/</code>).</li>
</ul>
<h2>9. Limitations</h2>
<ul>
<li>The SNP set is an LD-pruned array-based panel; it has no power to find narrow selection signals and is not a catalogue of all differences.</li>
<li>The reference panel has no Central Asian populations, so &ldquo;closest&rdquo; means closest among EUR, SAS, EAS and AFR.</li>
<li>Reference individuals were not filtered for relatedness; the effect on the weighted values was not tested (the values changed by at most 0.0027 between two quite different reference files, so a large effect is unlikely).</li>
<li>The jackknife standard error is for the mean per-SNP value; PLINK does not give one for the weighted value.</li>
</ul>
"""
t = put(t, "step9", blk, replace_range=("<!-- ==================== OVERVIEW ==================== -->", "<!-- ==================== HISTORICAL HG19 PIPELINE (SECTIONS 2-10) ==================== -->"))
t = old_wrap(t, "step9", "<!-- ==================== HISTORICAL HG19 PIPELINE (SECTIONS 2-10) ==================== -->",
             "            </div>\n        </div>\n    </div>\n    \n    <script>",
             "ARCHIVED, not current: original hg19 pipeline of 2025 (1,047/1,199 samples, 1000G phase 3). Superseded; kept for provenance only. Do not cite its numbers.", end_regex=r"</div>\s*</div>\s*</div>\s*<script>")
t = t.replace("✓ Expanded cohort, GRCh38 canonical — August 2026 (1,256 samples)", "✓ NYGC 30x reference, expanded cohort (1,256 samples) — rebuilt 2026-10-04")
t = t.replace('<span class="step-badge step-badge-old">Spring 2026 (1,047 samples) — April 11, 2026</span>', '')
if '<div class="info-box" style="margin:20px 30px 0;">' in t:
    i = t.index('<div class="info-box" style="margin:20px 30px 0;">')
    j = t.index("</div>", i) + 6
    t = t[:i] + t[j:]
t = banner(t, "<strong>Updated 2026-10-04:</strong> this page was rebuilt on the NYGC 30x reference. The previous banner (tables not yet re-verified) no longer applies. See the <a href=\"deprecated_old_runs.html\">list of current and deprecated runs</a>.")
write("steps/step9.html", t)
print("step9 written", len(t))
