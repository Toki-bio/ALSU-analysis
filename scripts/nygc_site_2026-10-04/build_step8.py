from sg import *
import re

t = read("steps/step8.html")
PC = R["pca"]; G = PC["groups"]; ev = PC["eigenval"]; pct = PC["pct10"]
kc = PC["knn_counts"]; lk = PC["label_by_knn"]
order = ["EUR", "SAS", "EAS", "AFR", "ALSU"]


def m(g, pc):
    mu, sd = G[g][pc]
    return f"{mu:+.4f} &plusmn; {sd:.4f}"


grp_rows = [[g if g != "ALSU" else "ALSU cohort", f"{G[g]['n']:,}"] + [m(g, f"PC{i}") for i in range(1, 5)] for g in order]
tot20 = sum(ev[:20])
ev_rows = [[i + 1, f"{ev[i]:.2f}", f"{pct[i]:.1f}%" if i < 10 else "&mdash;", f"{100 * ev[i] / tot20:.1f}%"] for i in range(12)]
lab_names = {"Uzbek": "self-reported Uzbek", "UZB_(mixed/unknown)": "Uzbek, mixed/unknown ethnicity", "UZB_(no_phenotype)": "Uzbek cohort, no phenotype record", "Karakalpak": "Karakalpak", "Kazakh": "Kazakh", "Korean": "Korean", "Russian": "Russian", "Tajik": "Tajik", "Tatar": "Tatar"}
lab_order = ["Uzbek", "UZB_(mixed/unknown)", "UZB_(no_phenotype)", "Tajik", "Tatar", "Russian", "Kazakh", "Karakalpak", "Korean"]
lab_rows = []
for k in lab_order:
    d = lk.get(k, {})
    tot = sum(d.values())
    lab_rows.append([lab_names[k], f"{tot:,}"] + [f"{d.get(p, 0):,} ({100 * d.get(p, 0) / tot:.0f}%)" for p in ("EUR", "SAS", "EAS", "AFR")])
tot_c = PC["n_cohort"]

blk = f"""
<h2>1. Overview</h2>
<p>Principal component analysis (PCA) of genotypes places each person on axes that summarise the main differences in ancestry. Here the 1,256 Uzbek samples are analysed together with 2,712 reference samples from the NYGC 30x 1000 Genomes data (EUR, EAS, SAS, AFR; AMR is not used), so that the position of the Uzbek cohort can be read against known continental groups.</p>
{box("warn", "<strong>This page was rebuilt on 2026-10-04.</strong> The earlier version used the 1000G phase 3 GRCh38 shapeit2 file (2,548 reference samples including AMR, 83,091 SNPs, August 2026). It is superseded: the new run uses the NYGC 30x reference. The first two components look the same, but all numbers on this page are new. The old pipeline text is kept collapsed at the bottom.")}
{cards([(f"{PC['n_total']:,}", "samples in the PCA"), (f"{tot_c:,}", "Uzbek cohort samples"), (f"{PC['n_ref']:,}", "NYGC reference samples"), ("64,959", "LD-pruned SNPs"), (f"{pct[0]:.1f}% / {pct[1]:.1f}%", "PC1 / PC2 (of the first 10)")])}

<h2>2. Data and method</h2>
<ul>
<li><strong>Samples:</strong> UZB 1,256; EUR 633; EAS 585; SAS 601; AFR 893 (NYGC superpopulation codes). Reference samples were not filtered for relatedness.</li>
<li><strong>SNPs:</strong> 64,959 SNPs from the 82,744-SNP merged panel after LD pruning (the same input as the ADMIXTURE runs, <a href="step11.html">step 11</a>): <code>/staging/ALSU-analysis/spring2026/full_expanded_cohort/nygc30x_1256/adm/adm_in.*</code> on DRAGEN.</li>
<li><strong>PCA:</strong> PLINK 1.9 <code>--pca 20</code> on that file (<code>pca.eigenvec</code>, <code>pca.eigenval</code> in the same directory). Percentages below are shares of the sum of the first 10 eigenvalues (PLINK&rsquo;s <code>--pca</code> does not report the total variance).</li>
<li><strong>Nearest-reference label:</strong> for each Uzbek sample, the 25 nearest reference samples in PC1&ndash;PC4 (Euclidean) were found and the most frequent reference group was recorded. This is a descriptive label of where a sample lies; it is <em>not</em> an ancestry proportion (see <a href="step11.html">step 11</a> for ADMIXTURE).</li>
<li>Figures and tables on this page were produced by <code>pca_analysis.py</code> from the DRAGEN files above (copy in <code>scripts/nygc_site_2026-10-04/</code>).</li>
</ul>

<h2>3. Results</h2>
{img("nygc_global_pca_panels.png", "Global PCA. Black: Uzbek cohort. Coloured: NYGC reference groups. Upper right: the Eurasian part of PC1 vs PC2 enlarged.")}
<h3>Eigenvalues</h3>
{img("nygc_global_pca_scree.png", "Eigenvalues of the first 20 components (log scale).")}
{table(["PC", "Eigenvalue", "Share of the first 10", "Share of the first 20"], ev_rows)}
<p>Two bases are used for &ldquo;percent of variance&rdquo; in this project. The report and this page quote the share of the first 10 eigenvalues (PC1 {pct[0]:.1f}%, PC2 {pct[1]:.1f}%); <a href="step11.html">step 11</a> quotes the share of the first 20 (PC1 {100*ev[0]/tot20:.1f}%, PC2 {100*ev[1]/tot20:.1f}%). They describe the same PCA. Neither is a share of the total genetic variance, which PLINK&rsquo;s <code>--pca</code> does not report.</p>
<p>Three components stand out (214, 91 and 29), PC4&ndash;PC10 are small (8 down to 2.7) and from PC11 the eigenvalues form a plateau near 2.6. Only the first few components therefore carry population structure in this panel.</p>
<h3>Group means and spread (PC1&ndash;PC4)</h3>
{table(["Group", "n", "PC1", "PC2", "PC3", "PC4"], grp_rows)}
{box("key", f"<strong>What the components mean.</strong> PC1 separates AFR from all Eurasians (AFR mean {G['AFR']['PC1'][0]:+.3f}; every other group near {G['EUR']['PC1'][0]:+.3f}). PC2 separates EUR ({G['EUR']['PC2'][0]:+.3f}) from EAS ({G['EAS']['PC2'][0]:+.3f}), with SAS in between ({G['SAS']['PC2'][0]:+.3f}). PC3 separates SAS from the rest. The Uzbek cohort has almost no spread on PC1 (SD {PC['coh_pc12_sd'][0]:.4f}, i.e. almost no African-like ancestry) and a large spread on PC2 (SD {PC['coh_pc12_sd'][1]:.4f}, mean {G['ALSU']['PC2'][0]:+.4f}): it forms a band along the EUR&ndash;SAS&ndash;EAS axis, centred near SAS, and reaching both the European cluster and the East Asian cluster.")}

{box("info", f"<strong>PC4 separates the whole Uzbek cohort from every reference group.</strong> The cohort mean on PC4 is {G['ALSU']['PC4'][0]:+.4f} (SD {G['ALSU']['PC4'][1]:.4f}), while the reference groups lie between {min(G[g]['PC4'][0] for g in ('EUR','SAS','EAS','AFR')):+.4f} and {max(G[g]['PC4'][0] for g in ('EUR','SAS','EAS','AFR')):+.4f}. Two explanations are possible and this analysis cannot separate them: (a) an ancestry component common in the Uzbek cohort and absent from the reference panel (ADMIXTURE at K=7 also finds an Uzbek-specific component, <a href='step11.html'>step 11</a>), or (b) a technical difference between the Uzbek data (array genotypes imputed to sequence level) and the sequencing-based reference. A test would be to project Uzbek-like samples genotyped on the same platform, or to add Central Asian reference genomes; neither is available in this project.")}
<h2>4. Where the Uzbek samples fall</h2>
{table(["Nearest-reference group", "Uzbek samples", "Share"], [[p, f"{kc.get(p, 0):,}", f"{100 * kc.get(p, 0) / tot_c:.1f}%"] for p in ("SAS", "EUR", "EAS", "AFR")])}
<p>The label SAS is the largest class because the centre of the Uzbek band lies next to the South Asian cluster, not because most Uzbek people are of South Asian origin. {PC['knn_mixed']:,} samples have neighbours from more than one reference group, which is the usual picture for an admixed sample.</p>
<h3>Nearest-reference label by self-reported group</h3>
{table(["Self-reported group (phenotype sheet / mapping)", "n", "EUR", "SAS", "EAS", "AFR"], lab_rows)}
<p>Self-reported ethnicity and the PCA position agree in the expected direction: Russian and Tatar participants lie mostly in the European part, Tajik participants mostly in the South Asian part, Korean participants in the East Asian cluster, and the Uzbek groups spread across all three. The association between self-report, birthplace and ancestry components is tested formally on <a href="step11.html">step 11</a>.</p>
{box("warn", "<strong>One sample lies inside the African reference cluster.</strong> One Uzbek-cohort sample sits at PC1 = " + f"{PC['afr_like_ids_pc'][0][0]:+.3f}" + " among the AFR reference samples and has an ADMIXTURE K=4 African-like component of 0.90. No other cohort sample has PC1 &gt; 0.01. Possible explanations are a person of recent African ancestry, a sample mix-up, or a labelling error; the data on this page cannot tell which. It has not been excluded from any analysis; provenance should be checked against the sample sheet. (Its internal code is not shown on this public page.)")}

<h2>5. What this means for the rest of the project</h2>
<ul>
<li>The cohort is structured mainly along one axis (EUR&ndash;SAS&ndash;EAS), so population stratification is a real concern in association tests. GWAS covariates are the PCs from the <em>Uzbek-only</em> PCA (step 7), which do not depend on the reference panel (<a href="step16.html">step 16</a>).</li>
<li>For ancestry proportions use ADMIXTURE (<a href="step11.html">step 11</a>); the PCA labels above are only a map.</li>
</ul>
<h2>6. What changed</h2>
{table(["Earlier version (August 2026)", "Now"], [["3,804 samples (1,256 UZB + 2,548 1000G phase 3 including AMR); 83,091 SNPs; eigenvalues 173, 75, 25, 20, 7.9 &hellip;", f"{PC['n_total']:,} samples (1,256 UZB + 2,712 NYGC without AMR); 64,959 SNPs; eigenvalues 214, 91, 29, 7.8, 5.3 &hellip;"], ["Reference: ALL.chr*.shapeit2_integrated_v1a.GRCh38 (batch artefact at a few dozen sites)", "Reference: NYGC 30x high-coverage GRCh38, md5-verified"], ["No statement on individual outliers", "One sample inside the African cluster flagged for provenance check"]])}
<h2>7. Limitations</h2>
<ul>
<li>PCA is a projection; two or three components do not capture all structure, and the share of variance is relative to the first 10 components only.</li>
<li>The reference panel has no Central Asian groups, so the position of the Uzbek band is described relative to EUR, SAS and EAS only.</li>
<li>The reference group sizes differ and relatives are included in the reference; PCA axes are driven by the largest groups.</li>
</ul>
"""
i = t.index('<h2 class="section-header">1. Overview</h2>')
j = t.rfind('<div class="section">', 0, i)
start_lit = t[j:i]
t = old_wrap(t, "step8", start_lit, None,
             "ARCHIVED, not current: previous version of this page (August 2026; 1000G phase 3 file, 2,548 reference samples incl. AMR, 83,091 SNPs). Superseded; do not cite its numbers.",
             end_regex=r"</div>\s*</div>\s*</div>\s*<script>")
t = put(t, "step8", blk, before="<!--NYGC:step8:OLDBEGIN-->")
t = t.replace("✓ Expanded cohort — August 2026 (1,256 Uzbek samples)", "✓ NYGC 30x reference, expanded cohort (1,256 Uzbek samples) — rebuilt 2026-10-04")
t = t.replace('<span class="step-badge step-badge-old">Spring 2026 (1,047 samples) — April 11, 2026</span>', '')
if '<div class="info-box" style="margin:20px 30px 0;">' in t:
    a = t.index('<div class="info-box" style="margin:20px 30px 0;">'); b = t.index("</div>", a) + 6
    t = t[:a] + t[b:]
t = banner(t, "<strong>Updated 2026-10-04:</strong> rebuilt on the NYGC 30x reference. See the <a href=\"deprecated_old_runs.html\">list of current and deprecated runs</a>.")
write("steps/step8.html", t)
print("step8 ok")
