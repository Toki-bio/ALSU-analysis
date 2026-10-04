from sg import *
import re

t = read("steps/step14.html")
P = R["pairs"]; M = R["matrix"]; C22 = R["chr22"]
pops = ["UZB", "SAS", "EUR", "EAS", "AFR"]


def w(a, b):
    return M[a][b]


def mat_row(a):
    cells = ""
    for b in pops:
        cells += '<td style="background:#e8e8e8;text-align:center">&mdash;</td>' if a == b else f"<td>{w(a, b):.4f}</td>"
    return f"<tr><td><strong>{a}</strong></td>{cells}</tr>"


matrix = "<table id=\"fstTable\"><thead><tr><th></th>" + "".join(f"<th>{p}</th>" for p in pops) + "</tr></thead><tbody>" + "".join(mat_row(a) for a in pops) + "</tbody></table>"
close = sorted([(w("UZB", b), b) for b in pops if b != "UZB"])
cat = lambda v: ("Very low", "#2e7d32") if v < 0.02 else ("Moderate", "#e65100") if v < 0.1 else ("High", "#c62828")
rank_rows = "".join(f'<tr><td>{i + 1}</td><td><strong>{b}</strong></td><td>{v:.4f}</td><td style="color:{cat(v)[1]}">{cat(v)[0]}</td></tr>' for i, (v, b) in enumerate(close))
allp = sorted([(w(a, b), a, b) for i, a in enumerate(pops) for b in pops[i + 1:]], reverse=True)[:5]
div_rows = "".join(f'<tr><td>{i + 1}</td><td><strong>{a} &ndash; {b}</strong></td><td>{v:.4f}</td></tr>' for i, (v, a, b) in enumerate(allp))
mds = R["mds"]; mv = R["mds_var_pct"]

# chr22 check table
H = C22["hudson"]; WC = C22["wc"]
def hv(a, b):
    k = f"{a}_{b}"
    return H[k]["hudson"] if k in H else H[f"{b}_{a}"]["hudson"]
def nhv(a, b):
    k = f"{a}_{b}"
    return H[k]["n"] if k in H else H[f"{b}_{a}"]["n"]
def wcv(a, b):
    for k in (f"{a}_{b}", f"{b}_{a}"):
        if k in WC:
            return WC[k]["weighted"]
    return None
chk_pairs = [("UZB", "SAS"), ("UZB", "EUR"), ("UZB", "EAS"), ("UZB", "AFR"), ("EUR", "SAS"), ("EAS", "SAS"), ("EUR", "EAS"), ("AFR", "SAS"), ("AFR", "EUR"), ("AFR", "EAS")]
chk_rows = []
for a, b in chk_pairs:
    wcx = wcv(a, b)
    chk_rows.append([f"{a} vs {b}", f"{w(a, b):.4f}", "&mdash;" if wcx is None else f"{wcx:.4f}", f"{hv(a, b):.4f}", f"{nhv(a, b):,}", f"{(hv(a, b) / w(a, b)):.2f}"])

blk = f"""
<h2 id="current" style="border-bottom:3px solid #16a34a;">CURRENT RESULTS (NYGC 30x reference; rebuilt 2026-10-04)</h2>
{box("warn", "<strong>Everything on this page now refers to the NYGC 30x high-coverage reference.</strong> The August 2026 matrix (old 1000G file), the &ldquo;published F<sub>ST</sub> values&rdquo; table (never sourced) and the old sample sizes are withdrawn; the interactive heatmap and MDS below are drawn from the new matrix.")}

<h2 id="overview">1. Overview</h2>
<p>F<sub>ST</sub> between each of the ten pairs of five populations: <strong>UZB</strong> (Uzbek, n = 1,256), <strong>EUR</strong> (633), <strong>SAS</strong> (601), <strong>EAS</strong> (585) and <strong>AFR</strong> (893), computed with PLINK 1.9 (Weir &amp; Cockerham, <code>--fst --within</code>) on 82,744 LD-pruned SNPs shared between the Uzbek cohort and the NYGC 30x GRCh38 reference (AMR not used). The per-SNP analysis, top loci and the genome-wide distribution are on <a href="step9.html">step 9</a>; this page presents the matrix, its geometry and a check against unascertained SNPs.</p>
{cards([("10", "population pairs"), ("82,744", "LD-pruned SNPs"), ("3,968", "samples"), (f"{w('UZB', 'SAS'):.4f}", "min F<sub>ST</sub> (UZB&ndash;SAS)"), (f"{w('EAS', 'AFR'):.4f}", "max F<sub>ST</sub> (EAS&ndash;AFR)")])}
{box("key", f"<strong>Key finding.</strong> The Uzbek cohort is about equally close to South Asians ({w('UZB', 'SAS'):.4f}) and Europeans ({w('UZB', 'EUR'):.4f}); the distance to East Asians ({w('UZB', 'EAS'):.4f}) is {100 * w('UZB', 'EAS') / w('EUR', 'EAS'):.0f}% of EUR&ndash;EAS ({w('EUR', 'EAS'):.4f}) and {100 * w('UZB', 'EAS') / w('SAS', 'EAS'):.0f}% of SAS&ndash;EAS ({w('SAS', 'EAS'):.4f}), so Uzbeks lie between the Western Eurasian groups and East Asians. The distance to Africans ({w('UZB', 'AFR'):.4f}) is almost the same as SAS&ndash;AFR ({w('SAS', 'AFR'):.4f}). This agrees with ADMIXTURE K=4 (about 50% EUR-like, 30% EAS-like, 19% SAS-like; <a href='step11.html'>step 11</a>) and with a non-negative fit of Uzbek allele frequencies as a mixture of the references (EUR {R['nnls_mix']['EUR'] * 100:.0f}%, EAS {R['nnls_mix']['EAS'] * 100:.0f}%, SAS {R['nnls_mix']['SAS'] * 100:.0f}%; <a href='step9.html'>step 9</a>).")}

<!-- ════════════════════ 2. HEATMAP ════════════════════ -->
<h2 id="heatmap">2. Interactive F<sub>ST</sub> Heatmap</h2>
<p>Hover over cells to see the population pair and exact weighted F<sub>ST</sub>. Color scale: <strong style="color:#2e7d32">green</strong> (F<sub>ST</sub> &lt; 0.05) &rarr; <strong style="color:#e65100">orange</strong> (0.05&ndash;0.10) &rarr; <strong style="color:#c62828">red</strong> (&gt; 0.10).</p>
<div class="chart-wrap" style="position:relative">
  <canvas id="heatmapCanvas" width="620" height="560"></canvas>
  <div class="tooltip" id="heatmapTip"></div>
</div>

<!-- ════════════════════ 3. FULL MATRIX ════════════════════ -->
<h2 id="matrix">3. Full F<sub>ST</sub> Matrix</h2>
{matrix}
{box("info", "<strong>Reading the matrix.</strong> Weighted F<sub>ST</sub> (ratio of summed numerators to summed denominators over all SNPs). Symmetric. Populations are ordered by distance from UZB. The matrix uses all 82,744 SNPs; on the stricter set without the 8,573 sites with overlapping indel records no value differs by more than 0.0005 (<a href='step9.html'>step 9</a>, section 3).")}

<!-- ════════════════════ 4. MDS PLOT ════════════════════ -->
<h2 id="mds">4. Classical MDS from F<sub>ST</sub> Distance Matrix</h2>
<p>Classical (metric) multidimensional scaling projects the 5 &times; 5 distance matrix into two dimensions (eigendecomposition of the double-centred squared distances; F<sub>ST</sub> itself is used as the distance). The first axis carries {mv[0]:.1f}% and the second {mv[1]:.1f}% of the positive eigenvalues; the first separates AFR from everyone else, the second EUR from EAS, and UZB sits near SAS, slightly toward EAS relative to EUR. Because two dimensions cannot hold a 5-point distance matrix exactly, small distances such as UZB&ndash;SAS&ndash;EUR are only approximate in the plot; the table in section 3 is the exact result.</p>
<div class="two-col">
  <div class="chart-panel">
    <h3>MDS Projection (Dimension 1 vs 2)</h3>
    <canvas id="mdsCanvas" width="480" height="420"></canvas>
  </div>
  <div class="chart-panel">
    <h3>F<sub>ST</sub> Bar Chart &mdash; Distance from UZB</h3>
    <canvas id="barCanvas" width="480" height="420"></canvas>
  </div>
</div>
<p style="font-size:.85em;color:#666">Static version of the heatmap and MDS (computed in Python from the same matrix): <a href="../images/nygc_fst_heatmap_mds.png">nygc_fst_heatmap_mds.png</a>. MDS coordinates (dimension 1, dimension 2): {", ".join(f"{p} ({mds[p][0]:.4f}, {mds[p][1]:.4f})" for p in pops)}.</p>

<!-- ════════════════════ 5. POPULATION RANKING ════════════════════ -->
<h2 id="ranking">5. Population Proximity Ranking</h2>
<div class="two-col">
  <div>
    <h3>Closest to UZB</h3>
    <table><thead><tr><th>#</th><th>Population</th><th>Weighted F<sub>ST</sub></th><th>Category</th></tr></thead><tbody>{rank_rows}</tbody></table>
  </div>
  <div>
    <h3>Largest divergences (all pairs)</h3>
    <table><thead><tr><th>#</th><th>Pair</th><th>Weighted F<sub>ST</sub></th></tr></thead><tbody>{div_rows}</tbody></table>
  </div>
</div>

<!-- ════════════════════ 6. INTERPRETATION ════════════════════ -->
<h2 id="interpretation">6. Interpretation</h2>
<h3>6.1. Where the Uzbek cohort sits</h3>
<ul>
<li><strong>UZB&ndash;SAS ({w('UZB', 'SAS'):.4f}) and UZB&ndash;EUR ({w('UZB', 'EUR'):.4f}) are nearly equal.</strong> The difference (0.0010) is measurable but small (<a href="step9.html">step 9</a>); which of the two is &ldquo;closer&rdquo; depends on the reference panel, which contains no Central Asian groups.</li>
<li><strong>UZB&ndash;EAS ({w('UZB', 'EAS'):.4f}) is less than half of EUR&ndash;EAS ({w('EUR', 'EAS'):.4f}) and smaller than SAS&ndash;EAS ({w('SAS', 'EAS'):.4f}).</strong> That is the F<sub>ST</sub> expression of an East-Asian-like share of ancestry that Europeans and South Asians lack.</li>
<li><strong>UZB&ndash;AFR ({w('UZB', 'AFR'):.4f}) equals SAS&ndash;AFR ({w('SAS', 'AFR'):.4f}) and is below EUR&ndash;AFR ({w('EUR', 'AFR'):.4f}),</strong> as expected when there is no African ancestry component (ADMIXTURE K=4: about 0%).</li>
</ul>
<p>Historical explanations for these patterns (Indo-Iranian, Steppe, Turkic or Mongol contributions) are plausible context but were not tested here: the panel has no ancient or Central Asian reference samples.</p>
<h3>6.2. Continental structure</h3>
<ul>
<li>AFR is farthest from every other group (0.129&ndash;0.167), EAS is farther from SAS (0.057) and EUR (0.087) than UZB is from either.</li>
<li>EUR&ndash;SAS ({w('EUR', 'SAS'):.4f}) is the closest pair among the three reference Eurasian groups, and both are close to UZB.</li>
</ul>
<h3>6.3. Check against unascertained SNPs (replaces the earlier &ldquo;published values&rdquo; table)</h3>
<p>The earlier version of this page compared reference-pair values with &ldquo;published 1000 Genomes F<sub>ST</sub> ranges&rdquo;. Those ranges have no source in the project records and are removed. In their place, F<sub>ST</sub> was recomputed on chromosome 22 using <em>all</em> common biallelic SNVs (MAF &ge; 0.01; about 190,000 in the NYGC data, 60,979 in the imputed Uzbek data, 53,765 shared) instead of the LD-pruned array-based panel. Reference pairs: PLINK Weir &amp; Cockerham; Hudson (Bhatia et al. 2013) estimator computed from allele frequencies for all pairs, used for the Uzbek pairs because the Uzbek and NYGC data are not in one PLINK file. The Hudson estimator reproduces PLINK on the reference pairs to within 0.005.</p>
{table(["Pair", "Panel (82,744 SNPs, W&amp;C)", "chr22, all SNVs (W&amp;C)", "chr22, all SNVs (Hudson)", "SNPs used (Hudson)", "Hudson / panel"], chk_rows)}
{box("info", "<strong>Two conclusions.</strong> (1) <em>The ordering is the same</em> in the panel and in the unascertained chr22 set: UZB is about equally close to SAS and EUR, closer to EAS than SAS and EUR are, and equal to SAS in distance from AFR. (2) <em>The magnitudes are not the same.</em> The array-based, LD-pruned panel gives smaller F<sub>ST</sub> between Eurasian groups than all common SNVs (for example EUR&ndash;EAS 0.087 vs 0.104) and slightly larger F<sub>ST</sub> to AFR for EUR (0.141 vs 0.126). This is SNP ascertainment: array SNPs are chosen to be common in the populations the array was designed for. Absolute F<sub>ST</sub> values from this panel should therefore not be compared with published genome-wide values without that caveat. Chromosome 22 is a single chromosome, so these numbers carry sampling noise of the order of a few thousandths.")}

<!-- ════════════════════ 7. METHODS ════════════════════ -->
<h2 id="methods">7. Methods</h2>
{box("meth", "<strong>F<sub>ST</sub> estimation:</strong> PLINK v1.9 <code>--fst --within</code> (Weir &amp; Cockerham 1984) on the merged cohort + NYGC panel (<code>merged_refqc</code>, 82,744 SNPs, GRCh38). &ldquo;Weighted&rdquo; is the ratio-of-averages estimator.")}
{box("meth", "<strong>MDS:</strong> classical metric MDS (double-centring of the squared distance matrix, two leading eigenvectors), implemented in JavaScript on this page and re-computed in Python for the static figure; both use the F<sub>ST</sub> values directly as distances.")}
<h3>7.1. Sample sizes</h3>
{table(["Population", "N", "Source"], [["UZB (Uzbek)", "1,256", "Expanded ALSU cohort (post-QC)"], ["EUR", "633", "NYGC 1000G 30x (CEU, GBR, FIN, IBS, TSI)"], ["EAS", "585", "NYGC 1000G 30x (CHB, JPT, CHS, CDX, KHV)"], ["SAS", "601", "NYGC 1000G 30x (GIH, PJL, BEB, STU, ITU)"], ["AFR", "893", "NYGC 1000G 30x (YRI, LWK, GWD, MSL, ESN, ACB, ASW)"]])}
<p style="font-size:.88em;color:#666">Reference samples were not filtered for relatedness. Population labels are the superpopulation codes of the NYGC sample table, as used in the keep-lists on DRAGEN.</p>
<h3>7.2. Commands and files</h3>
<pre style="background:#1e1e1e;color:#d4d4d4;padding:14px 16px;border-radius:8px;overflow-x:auto;font-size:.82em">for each pair A_B (UZB EUR EAS SAS AFR):
  plink --bfile merged_refqc --keep keep_A_B.txt --within within_A_B.txt --fst --out all_A_B
  (second set: --exclude nygc30x_1256_filt/flagged_overlap.txt  -> filt_A_B)</pre>
<ul>
<li>Script <code>/staging/tmp/scratch/alsu_fst_all.sh</code>; outputs <code>/staging/ALSU-analysis/spring2026/full_expanded_cohort/pbs_refqc_2026-10/fst_all_pairs_2026-10-04/</code> (DRAGEN).</li>
<li>chr22 check: <code>/staging/tmp/fst22.sh</code> and <code>/staging/tmp/fst22b.py</code>; outputs <code>/staging/tmp/scratch/fst22_check/</code> (DRAGEN).</li>
<li>Run date 2026-10-04 (matrix first computed 2026-10-02/04).</li>
</ul>

"""
a = '<h2 style="border-bottom:3px solid #16a34a;">CURRENT RESULTS (NYGC 30x reference, 2026-10-04)</h2>'
b = "</div><!-- /content -->"
t = put(t, "step14", blk, replace_range=(a, b)) if "NYGC:step14:BEGIN" not in t else put(t, "step14", blk, replace_range=("<!--NYGC:step14:BEGIN-->", b))
# JS matrix
js_old = re.search(r"const fst = \[.*?\];", t, flags=re.S).group(0)
js_new = "const fst = [\n" + ",\n".join("    [" + ", ".join(f"{(0 if a == b2 else w(a, b2)):.4f}" for b2 in pops) + "]" for a in pops) + "\n];"
t = t.replace(js_old, js_new)
t = t.replace("// Expanded cohort (1,256 UZB samples), verified directly against DRAGEN\n// step10_pbs/fst_*_fixed.log and step14_fst_matrix/fst_*.log, Aug 13 2026.", "// NYGC 30x reference, expanded cohort (1,256 UZB); weighted Fst, 82,744 SNPs.\n// Source: DRAGEN pbs_refqc_2026-10/fst_all_pairs_2026-10-04/all_*.log (2026-10-04).")
t = t.replace("🧬 Expanded cohort &bull; 5 populations &bull; 83,091 SNPs &bull; August 2026", "🧬 NYGC 30x reference &bull; 5 populations &bull; 82,744 SNPs &bull; rebuilt 2026-10-04")
t = t.replace("<title>Step 14: Fst Heatmap &amp; MDS — ALSU Pipeline</title>", "<title>Step 14: Fst Heatmap &amp; MDS — ALSU Pipeline</title>")

# stale hard-coded values in the page script
t = t.replace("const vals = [0.0144, 0.0145, 0.0393, 0.1293];", "const vals = [" + ", ".join(f"{w('UZB', b2):.4f}" for b2 in ['SAS', 'EUR', 'EAS', 'AFR']) + "];")
t = t.replace("const pctDim1 = (vals[0] / vals.reduce((a,b)=>a+b,0) * 100).toFixed(1);", f"const pctDim1 = '{mv[0]:.1f}';")
t = t.replace("const pctDim2 = (vals[1] / vals.reduce((a,b)=>a+b,0) * 100).toFixed(1);", f"const pctDim2 = '{mv[1]:.1f}';")
t = t.replace("Dimension 1 (${pctDim1}%)", "Dimension 1 (${pctDim1}% of positive eigenvalues)").replace("Dimension 2 (${pctDim2}%)", "Dimension 2 (${pctDim2}%)")
write("steps/step14.html", t)
print("step14 ok")
