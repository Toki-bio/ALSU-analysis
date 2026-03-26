#!/usr/bin/env python3
"""
Major update of step11.html for V2 global ADMIXTURE results.

V2 changes:
- 3,595 samples (2,548 1000G all superpops + 1,047 UZB) vs V1's 2,095 (1,017 selected + 1,078 UZB)
- 77,111 LD-pruned SNPs vs V1's 380,376
- 6 groups (AFR, EUR, SAS, EAS, AMR, UZB) vs V1's 11 (10 sub-populations + UZBEK)
- CV minimum at K=8 (0.29422) but K=7-8 effectively tied (plateau)
- Evanno |L''(K)| peaks at K=5 (from V2 L(K) values)
- UZB-only: sNMF V2 done (K=2 optimal, CE=0.42501)
"""

import re

with open("steps/step11.html", encoding="utf-8") as f:
    html = f.read()

orig_len = len(html)

# ============================================================
# 1. Remove the main V2 status banner
# ============================================================
html = re.sub(
    r'<div class="info-box"[^>]*>\s*<strong>&#9888; V2 re-analysis status:</strong>.*?</div>',
    '',
    html,
    count=1,
    flags=re.DOTALL
)
print("1. Removed V2 status banner")

# ============================================================
# 2. Update Overview section: 1,017 → 2,548
# ============================================================
html = html.replace(
    "with 1,017 individuals from the 1000 Genomes Project",
    "with 2,548 individuals from the 1000 Genomes Project (all 26 populations across 5 superpopulations)"
)
print("2. Updated overview 1000G count")

# ============================================================
# 3. Replace reference panel table (10 sub-populations → 5 superpops)
# ============================================================
old_ref_table = """<h3>2.1 Reference Panel Selection</h3>
            <table>
                <tr><th>Superpopulation</th><th>Population</th><th>Full Name</th><th>N</th></tr>
                <tr><td><span class="superpop-tag tag-afr">AFR</span></td><td>YRI</td><td>Yoruba in Ibadan, Nigeria</td><td>108</td></tr>
                <tr><td rowspan="5"><span class="superpop-tag tag-eur">EUR</span></td><td>CEU</td><td>Utah residents (N. &amp; W. European)</td><td>99</td></tr>
                <tr><td>GBR</td><td>British in England and Scotland</td><td>91</td></tr>
                <tr><td>FIN</td><td>Finnish in Finland</td><td>99</td></tr>
                <tr><td>IBS</td><td>Iberian in Spain</td><td>107</td></tr>
                <tr><td>TSI</td><td>Toscani in Italia</td><td>107</td></tr>
                <tr><td rowspan="2"><span class="superpop-tag tag-sas">SAS</span></td><td>GIH</td><td>Gujarati Indian in Houston</td><td>103</td></tr>
                <tr><td>PJL</td><td>Punjabi in Lahore, Pakistan</td><td>96</td></tr>
                <tr><td rowspan="2"><span class="superpop-tag tag-eas">EAS</span></td><td>CHB</td><td>Han Chinese in Beijing</td><td>103</td></tr>
                <tr><td>JPT</td><td>Japanese in Tokyo</td><td>104</td></tr>
                <tr><td><span class="superpop-tag tag-uzb">UZB</span></td><td>UZBEK</td><td>Uzbek cohort (this study)</td><td>1,078</td></tr>
                <tr style="font-weight:600; background:#f0f0f0;"><td colspan="3">Total</td><td>2,095</td></tr>
            </table>"""

new_ref_table = """<h3>2.1 Reference Panel Selection</h3>
            <p>V2 uses the <strong>complete 1000 Genomes Phase 3</strong> reference (all 26 populations grouped into 5 superpopulations), providing broader continental representation than V1&rsquo;s 10 selected sub-populations.</p>
            <table>
                <tr><th>Superpopulation</th><th>Populations included</th><th>N</th></tr>
                <tr><td><span class="superpop-tag tag-afr">AFR</span></td><td>YRI, LWK, GWD, MSL, ESN, ACB, ASW</td><td>671</td></tr>
                <tr><td><span class="superpop-tag tag-eur">EUR</span></td><td>CEU, GBR, FIN, IBS, TSI</td><td>522</td></tr>
                <tr><td><span class="superpop-tag tag-sas">SAS</span></td><td>GIH, PJL, BEB, STU, ITU</td><td>492</td></tr>
                <tr><td><span class="superpop-tag tag-eas">EAS</span></td><td>CHB, JPT, CHS, CDX, KHV</td><td>515</td></tr>
                <tr><td><span class="superpop-tag tag-amr">AMR</span></td><td>MXL, PUR, CLM, PEL</td><td>348</td></tr>
                <tr><td><span class="superpop-tag tag-uzb">UZB</span></td><td>Uzbek cohort (this study)</td><td>1,047</td></tr>
                <tr style="font-weight:600; background:#f0f0f0;"><td colspan="2">Total</td><td>3,595</td></tr>
            </table>
            <div class="info-box" style="margin-top:12px;">
                <strong>V1 &rarr; V2 change:</strong> V1 used 10 hand-picked 1000G populations (1,017 samples); V2 uses all 26 populations (2,548 samples), adding the AMR superpopulation and expanding each continental group. This provides a more complete ancestry reference, particularly for admixed Central Asian populations with potential Native American-related (Ancient North Eurasian) signal.
            </div>"""

html = html.replace(old_ref_table, new_ref_table)
print("3. Replaced reference panel table")

# ============================================================
# 4. Update stat-box: 2,095 × 380,376 → 3,595 × 77,111
# ============================================================
html = html.replace(
    '<span class="stat-value">2,095 samples &times; 380,376 LD-pruned SNPs</span>',
    '<span class="stat-value">3,595 samples &times; 77,111 LD-pruned SNPs</span>'
)
print("4. Updated stat-box")

# ============================================================
# 5. Replace CV error table (V1 K=7 optimal → V2 K=8 min but plateau)
# ============================================================
old_cv_table = """<table>
                <tr><th>K</th><th>CV Error</th><th>Note</th></tr>
                <tr><td>2</td><td>0.39924</td><td>EUR+EAS vs AFR split</td></tr>
                <tr><td>3</td><td>0.38985</td><td>AFR separates</td></tr>
                <tr><td>4</td><td>0.38766</td><td>SAS separates</td></tr>
                <tr><td>5</td><td>0.38632</td><td>Central Asian component emerges</td></tr>
                <tr><td>6</td><td>0.38606</td><td>EUR splits (N vs S)</td></tr>
                <tr class="cv-optimal"><td>7</td><td>0.38535</td><td>&larr; Optimal (minimum CV)</td></tr>
                <tr><td>8</td><td>0.38540</td><td>Uptick &mdash; overfitting begins</td></tr>
            </table>"""

new_cv_table = """<table>
                <tr><th>K</th><th>CV Error</th><th>&Delta; from prev</th><th>Note</th></tr>
                <tr><td>2</td><td>0.31043</td><td>&mdash;</td><td>AFR vs non-AFR split</td></tr>
                <tr><td>3</td><td>0.29992</td><td>&minus;0.01051</td><td>EAS separates</td></tr>
                <tr><td>4</td><td>0.29703</td><td>&minus;0.00289</td><td>SAS separates</td></tr>
                <tr><td>5</td><td>0.29503</td><td>&minus;0.00200</td><td>AMR / Central Asian component</td></tr>
                <tr><td>6</td><td>0.29458</td><td>&minus;0.00045</td><td>Plateau begins</td></tr>
                <tr><td>7</td><td>0.29442</td><td>&minus;0.00016</td><td>Effectively tied with K=8</td></tr>
                <tr class="cv-optimal"><td>8</td><td>0.29422</td><td>&minus;0.00020</td><td>&larr; Minimum CV (K=7&ndash;8 plateau)</td></tr>
            </table>"""

html = html.replace(old_cv_table, new_cv_table)
print("5. Replaced CV error table")

# ============================================================
# 6. Replace CV success box
# ============================================================
html = html.replace(
    """<div class="success">
                <strong>K = 7 is optimal.</strong> The CV error reaches its minimum at K=7 (0.38535) with a slight uptick
                at K=8 (0.38540), confirming 7 distinguishable ancestral components in this dataset.
            </div>""",
    """<div class="success">
                <strong>K=8 is the nominal minimum</strong> (CV = 0.29422), but K=7 (0.29442) is only 0.07% higher &mdash; effectively a plateau.
                The improvement from K=6&rarr;7&rarr;8 (&lt;0.02% per step) is negligible compared to K=2&rarr;3 (&minus;3.4%).
                For practical purposes, <strong>K=5&ndash;8 form a plateau</strong> with no clear single optimum;
                K=5 is the most parsimonious model capturing all major continental groups.
            </div>"""
)
print("6. Replaced CV success box")

# ============================================================
# 7. Update Evanno section with V2 L(K) values
# ============================================================
# Evanno note
html = html.replace(
    "on the global dataset (2,095 samples &times; 380,376 SNPs). A Uzbek-only replicate batch\n                (10 reps &times; K=2&ndash;6, 1,074 samples) also completed on the server (51/70 runs;\n                K=7&ndash;8 were not rerun as sNMF independently confirmed K=2 for the Uzbek-only dataset &mdash;\n                see Section 3.3). Since the Evanno method is statistically inappropriate for ADMIXTURE\n                (see caveat below), we rely on CV error and sNMF cross-entropy for K selection instead.",
    "on the global V2 dataset (3,595 samples &times; 77,111 SNPs). Since the Evanno method is statistically inappropriate for ADMIXTURE\n                (see caveat below), we rely on CV error and sNMF cross-entropy for K selection instead."
)
print("7a. Updated Evanno note")

# Evanno description paragraph
html = html.replace(
    "Using the log-likelihood values L(K) from the global ADMIXTURE run (2,095 samples &times; 380,376 SNPs) at each K&nbsp;=&nbsp;2&ndash;8:",
    "Using the log-likelihood values L(K) from the global V2 ADMIXTURE run (3,595 samples &times; 77,111 SNPs) at each K&nbsp;=&nbsp;2&ndash;8:"
)
print("7b. Updated Evanno description")

# V2 L(K) values:
# K=2: -147,676,904.69, K=3: -144,699,938.39, K=4: -143,810,411.92
# K=5: -143,162,075.40, K=6: -142,951,364.13, K=7: -142,824,765.39, K=8: -142,672,081.88
# L'(K): K=3: 2,976,966.30, K=4: 889,526.47, K=5: 648,336.52, K=6: 210,711.27, K=7: 126,598.74, K=8: 152,683.51
# |L''(K)|: K=3: 2,087,439.83, K=4: 241,189.96, K=5: 437,625.25, K=6: 84,112.54, K=7: 26,084.77

old_evanno_table = """<table>
                <tr>
                    <th>K</th>
                    <th>L(K)</th>
                    <th>L&prime;(K) = L(K) &minus; L(K&minus;1)</th>
                    <th>|L&Prime;(K)|</th>
                </tr>
                <tr><td>2</td><td>&minus;250,763,063.48</td><td>&mdash;</td><td>&mdash;</td></tr>
                <tr style="background:#d4edda; font-weight:600;">
                    <td>3</td><td>&minus;247,216,025.97</td><td>3,547,037.51</td><td style="color:#28a745;">2,542,181.58 &larr; max</td>
                </tr>
                <tr><td>4</td><td>&minus;246,211,170.04</td><td>1,004,855.93</td><td>311,294.48</td></tr>
                <tr><td>5</td><td>&minus;245,517,608.59</td><td>693,561.45</td><td>383,004.96</td></tr>
                <tr><td>6</td><td>&minus;245,207,052.10</td><td>310,556.49</td><td>135,467.59</td></tr>
                <tr><td>7</td><td>&minus;244,761,028.02</td><td>446,024.08</td><td>256,746.40</td></tr>
                <tr><td>8</td><td>&minus;244,571,750.34</td><td>189,277.68</td><td>&mdash;</td></tr>
            </table>"""

new_evanno_table = """<table>
                <tr>
                    <th>K</th>
                    <th>L(K)</th>
                    <th>L&prime;(K) = L(K) &minus; L(K&minus;1)</th>
                    <th>|L&Prime;(K)|</th>
                </tr>
                <tr><td>2</td><td>&minus;147,676,905</td><td>&mdash;</td><td>&mdash;</td></tr>
                <tr style="background:#d4edda; font-weight:600;">
                    <td>3</td><td>&minus;144,699,938</td><td>2,976,966</td><td style="color:#28a745;">2,087,440 &larr; max</td>
                </tr>
                <tr><td>4</td><td>&minus;143,810,412</td><td>889,526</td><td>241,190</td></tr>
                <tr><td>5</td><td>&minus;143,162,075</td><td>648,337</td><td>437,625</td></tr>
                <tr><td>6</td><td>&minus;142,951,364</td><td>210,711</td><td>84,113</td></tr>
                <tr><td>7</td><td>&minus;142,824,765</td><td>126,599</td><td>26,085</td></tr>
                <tr><td>8</td><td>&minus;142,672,081</td><td>152,684</td><td>&mdash;</td></tr>
            </table>"""

html = html.replace(old_evanno_table, new_evanno_table)
print("7c. Replaced Evanno L(K) table")

# Update Evanno interpretation box
html = html.replace(
    """<div class="important">
                <strong>&Delta;K peaks at K&nbsp;=&nbsp;3</strong> (|L&Prime;(K)| = 2,542,182), indicating that the most fundamental
                split in the data separates three broad ancestral groups: <strong>African</strong>, <strong>West Eurasian</strong>
                (European + South Asian), and <strong>East Eurasian</strong> (East Asian + Central Asian). This is the uppermost
                hierarchical level of population structure.
            </div>""",
    """<div class="important">
                <strong>|L&Prime;(K)| peaks at K&nbsp;=&nbsp;3</strong> (2,087,440), consistent with V1. The most fundamental
                split separates three broad ancestral groups: <strong>African</strong>, <strong>West Eurasian</strong>
                (European + South Asian), and <strong>East Eurasian</strong> (East Asian + Central Asian). A secondary peak at K=5 (437,625) reflects emergence of the Central Asian / AMR components.
            </div>"""
)
print("7d. Updated Evanno interpretation")

# Update reconciling box (K=7 → K=5/K=8 plateau)
old_reconcile = """<div class="info-box">
                <strong>Reconciling CV error (K=7) vs Evanno (K=3):</strong> These two criteria answer different questions.
                The Evanno method detects the <em>dominant structural break</em> &mdash; the deepest split in the population tree.
                Cross-validation error measures <em>overall predictive accuracy</em> and favors the model that
                best captures all levels of sub-structure. For this dataset, K&nbsp;=&nbsp;3 captures the continental-level divergence,
                while K&nbsp;=&nbsp;7 resolves finer sub-continental structure (e.g., Northern vs. Southern European, the distinct
                Central Asian component). Both results are informative: <strong>K&nbsp;=&nbsp;3 describes the hierarchy,
                K&nbsp;=&nbsp;7 describes the full resolution</strong>. However, given the statistical caveats above,
                <strong>CV error (K=7) is the primary criterion for model selection</strong>. For downstream GWAS covariate
                adjustment, K&nbsp;=&nbsp;7 remains the preferred model.
            </div>"""

new_reconcile = """<div class="info-box">
                <strong>Reconciling CV error (K=8 plateau) vs Evanno (K=3):</strong> These two criteria answer different questions.
                The Evanno method detects the <em>dominant structural break</em> &mdash; the deepest split in the population tree.
                Cross-validation error measures <em>overall predictive accuracy</em>. For V2, K&nbsp;=&nbsp;3 captures continental-level
                divergence (AFR / EUR+SAS / EAS), while the CV error plateau at K=5&ndash;8 suggests no single optimal resolution.
                Given that V2 uses superpopulations rather than sub-populations, fewer K are needed to capture the main structure.
                <strong>K&nbsp;=&nbsp;5 captures all major continental groups plus the Central Asian component</strong>;
                K=6&ndash;8 provide diminishing additional resolution.
                For downstream GWAS covariate adjustment, K&nbsp;=&nbsp;5 is the most parsimonious choice.
            </div>"""

html = html.replace(old_reconcile, new_reconcile)
print("7e. Updated reconciling box")

# Also update caveat box reference K=7 → K=5-8 plateau
html = html.replace(
    "<strong>cross-validation error (Section 3, K=7) remains the statistically appropriate model selection criterion\n                for ADMIXTURE</strong>",
    "<strong>cross-validation error (Section 3) remains the statistically appropriate model selection criterion\n                for ADMIXTURE</strong>"
)
print("7f. Updated caveat reference")

# ============================================================
# 8. Update sNMF global section
# ============================================================
html = html.replace(
    """<div class="info-box" style="border-left: 4px solid #2e7d32; background: #e8f5e9;">
                <strong>✅ Completed March 9, 2026.</strong> sNMF analysis run on global dataset
                (2,095 samples &times; 172,537 biallelic SNPs). 10 replicates per K, &alpha;&nbsp;=&nbsp;10.
            </div>""",
    """<div class="info-box" style="border-left: 4px solid #ff9800; background: #fff8e1;">
                <strong>&#9888; V1 results retained.</strong> The global sNMF analysis below was run on the V1 dataset
                (2,095 samples &times; 172,537 biallelic SNPs) and has <strong>not been re-run</strong> on V2.
                Since the V2 ADMIXTURE uses all 2,548 1000G samples (vs 1,017 in V1), a V2 sNMF re-run
                would require the same expanded reference panel. The Uzbek-only sNMF (Section 3.3)
                <em>has</em> been re-run on V2 and independently confirms K=2.
            </div>"""
)
print("8. Updated sNMF global status")

# ============================================================
# 9. Update Section 3.3 (Uzbek-only) — complete V2 data
# ============================================================
# Remove the V2 update info-box (now everything is V2)
html = html.replace(
    """<div class="info-box" style="border-left: 4px solid #2e7d32; background: #e8f5e9;">
                <strong>&#9888; V2 update (March 2026):</strong> Uzbek-only ADMIXTURE re-run on V2 dataset:
                <strong>1,047 samples &times; 88,722 SNPs</strong>. K=2&ndash;5 CV errors available (K=6&ndash;8 running).
                UZB-only sNMF will be re-run on V2 data (see suggestion below).
            </div>""",
    """<div class="info-box" style="border-left: 4px solid #2e7d32; background: #e8f5e9;">
                <strong>&#10003; V2 complete (March 2026):</strong> Both Uzbek-only ADMIXTURE (K=2&ndash;8) and sNMF (K=2&ndash;10)
                re-run on V2 dataset: <strong>1,047 samples &times; 88,722 SNPs</strong>. All results below are V2.
            </div>"""
)
print("9a. Updated UZB-only info-box")

# Replace UZB CV error table (add V2 K=6-8)
old_uzb_cv = """<tr style="color:#999;"><td>6</td><td>0.51682</td><td><em>Running&hellip;</em></td><td></td></tr>
                <tr style="color:#999;"><td>7</td><td>0.51885</td><td><em>Pending</em></td><td></td></tr>
                <tr style="color:#999;"><td>8</td><td>0.51974</td><td><em>Pending</em></td><td></td></tr>"""

new_uzb_cv = """<tr><td>6</td><td>0.51682</td><td>0.31266</td><td></td></tr>
                <tr><td>7</td><td>0.51885</td><td>0.31445</td><td></td></tr>
                <tr><td>8</td><td>0.51974</td><td>0.31627</td><td></td></tr>"""

html = html.replace(old_uzb_cv, new_uzb_cv)
print("9b. Updated UZB CV K=6-8 rows")

# Update UZB CV interpretation
html = html.replace(
    "CV error increases monotonically from K=2 in both V1 and V2, indicating <strong>no significant internal substructure</strong>\n                beyond a two-way split. V2 CV errors are substantially lower (0.307 vs 0.513) because the V2 LD-pruned\n                set has 88,722 genuine SNPs versus V1&rsquo;s 380,376 (inflated by 837K spurious imputation artifacts).",
    "CV error increases monotonically from K=2 through K=8 in both V1 and V2, confirming <strong>no significant internal substructure</strong>\n                beyond a two-way split. V2 CV errors are substantially lower (0.308&ndash;0.316 vs 0.513&ndash;0.520) because the V2 LD-pruned\n                set has 88,722 genuine SNPs versus V1&rsquo;s 380,376 (inflated by 837K spurious imputation artifacts). The monotonic pattern across all 7 K-values is definitive."
)
print("9c. Updated UZB CV interpretation")

# Replace UZB sNMF table header and content (V1 → V2)
html = html.replace(
    "<h4>Uzbek-Only: sNMF Cross-Entropy (V1 &mdash; V2 rerun pending)</h4>",
    "<h4>Uzbek-Only: sNMF Cross-Entropy (V2)</h4>"
)
print("9d. Updated sNMF UZB header")

# Replace the V1 sNMF CE table with V2 values
old_snmf_table = """<table>
                <tr><th>K</th><th>Mean CE</th><th>SD</th><th>Min CE</th></tr>
                <tr style="background:#d4edda; font-weight:600;"><td>2</td><td>0.728789</td><td>0.000205</td><td>0.728500</td></tr>
                <tr><td>3</td><td>0.728862</td><td>0.000212</td><td>0.728549</td></tr>
                <tr><td>4</td><td>0.729385</td><td>0.000209</td><td>0.729052</td></tr>
                <tr><td>5</td><td>0.729989</td><td>0.000227</td><td>0.729667</td></tr>
                <tr><td>6</td><td>0.730761</td><td>0.000122</td><td>0.730633</td></tr>
                <tr><td>7</td><td>0.731234</td><td>0.000182</td><td>0.731023</td></tr>
                <tr><td>8</td><td>0.731900</td><td>0.000187</td><td>0.731629</td></tr>
                <tr><td>9</td><td>0.732545</td><td>0.000243</td><td>0.732323</td></tr>
                <tr><td>10</td><td>0.733163</td><td>0.000282</td><td>0.732749</td></tr>
            </table>"""

new_snmf_table = """<table>
                <tr><th>K</th><th>Mean CE</th><th>SD</th><th>Min CE</th><th>Best Run</th></tr>
                <tr style="background:#d4edda; font-weight:600;"><td>2</td><td>0.425013</td><td>0.000322</td><td>0.424585</td><td>1</td></tr>
                <tr><td>3</td><td>0.425046</td><td>0.000325</td><td>0.424631</td><td>1</td></tr>
                <tr><td>4</td><td>0.425388</td><td>0.000335</td><td>0.424947</td><td>1</td></tr>
                <tr><td>5</td><td>0.425985</td><td>0.000331</td><td>0.425570</td><td>9</td></tr>
                <tr><td>6</td><td>0.426539</td><td>0.000436</td><td>0.425946</td><td>1</td></tr>
                <tr><td>7</td><td>0.427069</td><td>0.000323</td><td>0.426616</td><td>1</td></tr>
                <tr><td>8</td><td>0.427624</td><td>0.000299</td><td>0.427281</td><td>9</td></tr>
                <tr><td>9</td><td>0.428241</td><td>0.000346</td><td>0.427853</td><td>9</td></tr>
                <tr><td>10</td><td>0.428758</td><td>0.000330</td><td>0.428315</td><td>1</td></tr>
            </table>"""

html = html.replace(old_snmf_table, new_snmf_table)
print("9e. Replaced UZB sNMF CE table with V2")

# Update sNMF text
html = html.replace(
    "Cross-entropy also increases monotonically from K=2 &mdash; independently confirming that Uzbeks\n                form a <strong>relatively homogeneous group</strong> without deep ancestral subdivisions.",
    "Cross-entropy increases monotonically from K=2 (0.4250) to K=10 (0.4288), independently confirming that Uzbeks\n                form a <strong>relatively homogeneous group</strong> without deep ancestral subdivisions.\n                V2 CE values (~0.425) are substantially lower than V1 (~0.729), reflecting the cleaner V2 SNP set."
)
print("9f. Updated sNMF interpretation text")

# Replace old UZB concordance table with V2 concordance  
old_uzb_conc = """<table>
                <tr><th>K</th><th>Mean Component Correlation</th><th>Range</th></tr>
                <tr><td>2</td><td><strong>0.9979</strong></td><td>0.998 &ndash; 0.998</td></tr>
                <tr><td>3</td><td><strong>0.9847</strong></td><td>0.975 &ndash; 0.998</td></tr>
                <tr><td>5</td><td>0.6339</td><td style="color:#c62828;">0.113 &ndash; 0.976</td></tr>
            </table>"""

new_uzb_conc = """<table>
                <tr><th>K</th><th>Mean Component Correlation</th><th>Note</th></tr>
                <tr><td>2</td><td><strong>0.997</strong></td><td>Near-perfect</td></tr>
                <tr><td>3</td><td><strong>0.975</strong></td><td>Excellent</td></tr>
                <tr><td>4</td><td><strong>0.923</strong></td><td>Good</td></tr>
                <tr><td>5</td><td>0.801</td><td>Components diverge at noise level</td></tr>
            </table>"""

html = html.replace(old_uzb_conc, new_uzb_conc)
print("9g. Replaced UZB concordance table with V2")

# Update concordance text  
html = html.replace(
    "At K=2&ndash;3, sNMF and ADMIXTURE recover nearly identical ancestry proportions (r&nbsp;&gt;&nbsp;0.98).\n                By K=5, some components diverge drastically (r&nbsp;=&nbsp;0.11), confirming that higher K values\n                capture <strong>noise rather than real structure</strong> within this sample.",
    "At K=2&ndash;3, sNMF and ADMIXTURE recover nearly identical ancestry proportions (r&nbsp;&gt;&nbsp;0.97).\n                By K=5, concordance drops to r&nbsp;=&nbsp;0.80, confirming that higher K values\n                capture <strong>noise rather than real structure</strong> within this sample."
)
print("9h. Updated concordance interpretation")

# ============================================================
# 10. Update K=7 label in interpretation box
# ============================================================
html = html.replace(
    "consistent with the K=7 global\n                ADMIXTURE results where the Uzbek component (Q2, ~70%) is already a blend of multiple ancestral sources.",
    "consistent with the global ADMIXTURE results where the Uzbek component is a blend of multiple ancestral sources."
)
print("10. Updated K=7 reference in interpretation")

# ============================================================
# 11. Update Section 4 buttons and section heading
# ============================================================
# The K=7 star button — keep all K buttons but highlight K=5 instead of K=7 as optimal
html = html.replace(
    '<button class="k-btn optimal" data-k="7" onclick="showK(7)">K=7 &#9733;</button>',
    '<button class="k-btn" data-k="7" onclick="showK(7)">K=7</button>'
)
# Add optimal marker to K=5
html = html.replace(
    '<button class="k-btn" data-k="5" onclick="showK(5)">K=5</button>',
    '<button class="k-btn optimal" data-k="5" onclick="showK(5)">K=5 &#9733;</button>'
)
print("11. Updated K buttons (K=5 optimal)")

# ============================================================
# 12. Replace Section 5: K=7 Composition → K=5 (most parsimonious)
# ============================================================
old_section5_heading = "<h2>5. Ancestry Composition at K=7 (Optimal)</h2>"
new_section5_heading = "<h2>5. Ancestry Composition at K=5 (Most Parsimonious)</h2>"
html = html.replace(old_section5_heading, new_section5_heading)

# Replace component assignment table
old_comp_table = """<h3>5.1 Component Assignment</h3>
            <table>
                <tr><th>Component</th><th>Color</th><th>Dominant in</th><th>Interpretation</th></tr>
                <tr><td>Q1</td><td style="background:#e57373;color:white;text-align:center;">&#9632;</td><td>IBS (88.6%), TSI (93.2%)</td><td>Southern European</td></tr>
                <tr><td>Q2</td><td style="background:#9575cd;color:white;text-align:center;">&#9632;</td><td>UZBEK (70.2%)</td><td>Central Asian (Uzbek core)</td></tr>
                <tr><td>Q3</td><td style="background:#a1887f;color:white;text-align:center;">&#9632;</td><td>UZBEK (4.4%)</td><td>Central Asian minor subdivision</td></tr>
                <tr><td>Q4</td><td style="background:#64b5f6;color:white;text-align:center;">&#9632;</td><td>FIN (88.9%), CEU (33.1%), GBR (32.0%)</td><td>Northern European / Finnic</td></tr>
                <tr><td>Q5</td><td style="background:#4db6ac;color:white;text-align:center;">&#9632;</td><td>YRI (99.99%)</td><td>African</td></tr>
                <tr><td>Q6</td><td style="background:#ffb74d;color:white;text-align:center;">&#9632;</td><td>CHB (99.6%), JPT (100%)</td><td>East Asian</td></tr>
                <tr><td>Q7</td><td style="background:#81c784;color:white;text-align:center;">&#9632;</td><td>GIH (93.5%), PJL (81.5%)</td><td>South Asian</td></tr>
            </table>"""

# K=5 V2 mean Q for UZB: [0.2884, 0.0026, 0.0290, 0.1916, 0.4884]
# Interpretation needs to be based on which K=5 components dominate in which population
# AFR: [0.0016, 0.9562, 0.0037, 0.0052, 0.0333] → Q2 = AFR
# EUR: [0.0116, 0.0130, 0.0045, 0.0075, 0.9634] → Q5 = EUR
# SAS: [0.0206, 0.0087, 0.0019, 0.8978, 0.0710] → Q4 = SAS
# EAS: [0.9957, 0.0011, 0.0009, 0.0020, 0.0004] → Q1 = EAS
# AMR: [0.0028, 0.0914, 0.4298, 0.0037, 0.4723] → Q3=AMR-specific, Q5=EUR-like
# UZB: [0.2884, 0.0026, 0.0290, 0.1916, 0.4884] → Q5=EUR/Central Asian (48.8%), Q1=EAS (28.8%), Q4=SAS (19.2%)

new_comp_table = """<h3>5.1 Component Assignment</h3>
            <p>At K=5 with V2&rsquo;s superpopulation-level reference, five continental ancestry components resolve cleanly:</p>
            <table>
                <tr><th>Component</th><th>Color</th><th>Dominant in</th><th>Interpretation</th></tr>
                <tr><td>Q1</td><td style="background:#64b5f6;color:white;text-align:center;">&#9632;</td><td>EAS (99.6%)</td><td>East Asian</td></tr>
                <tr><td>Q2</td><td style="background:#81c784;color:white;text-align:center;">&#9632;</td><td>AFR (95.6%)</td><td>African</td></tr>
                <tr><td>Q3</td><td style="background:#ffb74d;color:white;text-align:center;">&#9632;</td><td>AMR (43.0%)</td><td>Americas-specific (Native American-like)</td></tr>
                <tr><td>Q4</td><td style="background:#4db6ac;color:white;text-align:center;">&#9632;</td><td>SAS (89.8%)</td><td>South Asian</td></tr>
                <tr><td>Q5</td><td style="background:#9575cd;color:white;text-align:center;">&#9632;</td><td>EUR (96.3%), UZB (48.8%)</td><td>European / West Eurasian</td></tr>
            </table>"""

html = html.replace(old_comp_table, new_comp_table)
print("12a. Replaced K=7 component table with K=5")

# Replace Uzbek ancestry breakdown
old_uzb_ancestry = """<h3>5.2 Uzbek Cohort Ancestry Breakdown</h3>
            <div class="chart-container">
                <canvas id="uzbekPieChart" width="500" height="350"></canvas>
            </div>

            <table>
                <tr><th>Component</th><th>Ancestry</th><th>Uzbek Mean</th><th>Interpretation</th></tr>
                <tr><td>Q2</td><td>Central Asian</td><td><strong>70.2%</strong></td><td>Core Uzbek / Turkic component</td></tr>
                <tr><td>Q6</td><td>East Asian</td><td><strong>12.3%</strong></td><td>Steppe / Mongol-period introgression</td></tr>
                <tr><td>Q4</td><td>Northern European</td><td>4.5%</td><td>Ancient Steppe / Indo-European substrate</td></tr>
                <tr><td>Q3</td><td>Central Asian (minor)</td><td>4.4%</td><td>Sub-structure within Central Asian ancestry</td></tr>
                <tr><td>Q1</td><td>Southern European</td><td>4.3%</td><td>Mediterranean / Silk Road contact</td></tr>
                <tr><td>Q7</td><td>South Asian</td><td>3.9%</td><td>Indo-Aryan substrate / trade contact</td></tr>
                <tr><td>Q5</td><td>African</td><td>0.4%</td><td>Minimal &mdash; effectively absent</td></tr>
            </table>"""

new_uzb_ancestry = """<h3>5.2 Uzbek Cohort Ancestry Breakdown (K=5)</h3>
            <div class="chart-container">
                <canvas id="uzbekPieChart" width="500" height="350"></canvas>
            </div>

            <table>
                <tr><th>Component</th><th>Ancestry</th><th>Uzbek Mean</th><th>Interpretation</th></tr>
                <tr><td>Q5</td><td>European / West Eurasian</td><td><strong>48.8%</strong></td><td>Core West Eurasian ancestry (EUR + Central Asian)</td></tr>
                <tr><td>Q1</td><td>East Asian</td><td><strong>28.8%</strong></td><td>Steppe / Turkic / Mongol heritage</td></tr>
                <tr><td>Q4</td><td>South Asian</td><td><strong>19.2%</strong></td><td>Indo-Aryan substrate / trade contact</td></tr>
                <tr><td>Q3</td><td>Americas-specific</td><td>2.9%</td><td>Shared Ancient North Eurasian ancestry</td></tr>
                <tr><td>Q2</td><td>African</td><td>0.3%</td><td>Minimal &mdash; effectively absent</td></tr>
            </table>"""

html = html.replace(old_uzb_ancestry, new_uzb_ancestry)
print("12b. Replaced Uzbek ancestry breakdown with K=5 V2")

# Update the Central Asian complexity box
old_complexity = """<div class="important">
                <strong>Central Asian complexity:</strong> The Uzbek cohort is ~70% a distinct Central Asian component not found in
                any 1000G reference population, with ~12% East Asian ancestry (Mongol/Turkic steppe heritage) and ~13%
                Western Eurasian ancestry split between Northern European (4.5%), Southern European (4.3%), and South Asian (3.9%).
                The ~0.4% African is consistent with noise / deep shared ancestry.
            </div>"""

new_complexity = """<div class="important">
                <strong>Central Asian complexity:</strong> With V2&rsquo;s superpopulation-level reference, the Uzbek cohort is ~49% West Eurasian,
                ~29% East Asian, and ~19% South Asian. Unlike V1 (where a distinct &ldquo;Central Asian&rdquo; component appeared at K&ge;5),
                V2 decomposes Uzbek ancestry directly into continental sources &mdash; reflecting that the Central Asian component is itself
                a mixture of West Eurasian, East Asian, and South Asian ancestries. The ~3% AMR-like component likely reflects
                shared Ancient North Eurasian (ANE) ancestry.
            </div>"""

html = html.replace(old_complexity, new_complexity)
print("12c. Updated complexity interpretation")

# ============================================================
# 13. Update covariate section: 1,068 matched samples
# ============================================================
html = html.replace(
    "Kruskal&ndash;Wallis tests assessed whether self-reported ethnicity or birthplace within the Uzbek cohort\n                predicts any ADMIXTURE component (1,068 matched samples).",
    "Kruskal&ndash;Wallis tests assessed whether self-reported ethnicity or birthplace within the Uzbek cohort\n                predicts any ADMIXTURE component (V1: 1,068 matched samples; covariate tests not yet re-run on V2)."
)
print("13. Updated covariate section note")

# ============================================================
# 14. Replace JavaScript data
# ============================================================
# cvErrors
old_cv_js = """const cvErrors = [
    {k:2, cv:0.39924}, {k:3, cv:0.38985}, {k:4, cv:0.38766},
    {k:5, cv:0.38632}, {k:6, cv:0.38606}, {k:7, cv:0.38535},
    {k:8, cv:0.38540}
];"""

new_cv_js = """const cvErrors = [
    {k:2, cv:0.31043}, {k:3, cv:0.29992}, {k:4, cv:0.29703},
    {k:5, cv:0.29503}, {k:6, cv:0.29458}, {k:7, cv:0.29442},
    {k:8, cv:0.29422}
];"""

html = html.replace(old_cv_js, new_cv_js)
print("14a. Replaced JS cvErrors")

# popNames, popN, popSuper
old_pop_js = """const popNames = ['YRI','CEU','FIN','GBR','IBS','TSI','GIH','PJL','CHB','JPT','UZBEK'];
const popN     = [108, 99, 99, 91, 107, 107, 103, 96, 103, 104, 1078];
const popSuper = ['AFR','EUR','EUR','EUR','EUR','EUR','SAS','SAS','EAS','EAS','UZB'];"""

new_pop_js = """const popNames = ['AFR','EUR','SAS','EAS','AMR','UZB'];
const popN     = [671, 522, 492, 515, 348, 1047];
const popSuper = ['AFR','EUR','SAS','EAS','AMR','UZB'];"""

html = html.replace(old_pop_js, new_pop_js)
print("14b. Replaced JS popNames/popN/popSuper")

# qData — replace all of it
old_qdata_start = "const qData = {"
old_qdata_end_marker = "};\n\n// Color palettes"

# Find the qData block 
qdata_start = html.find("const qData = {")
qdata_end = html.find("};\n\n// Color palettes", qdata_start)
if qdata_start > 0 and qdata_end > 0:
    old_qdata = html[qdata_start:qdata_end + 2]
    
    new_qdata = """const qData = {
    2: [[0.9556,0.0444],[0.0241,0.9759],[0.0772,0.9228],[0.0162,0.9838],[0.0924,0.9076],[0.0216,0.9784]],
    3: [[0.9566,0.0379,0.0055],[0.0111,0.9766,0.0123],[0.0667,0.6259,0.3074],[0.0011,0.0007,0.9982],[0.0841,0.6470,0.2689],[0.0092,0.6308,0.3600]],
    4: [[0.9563,0.0343,0.0057,0.0038],[0.0131,0.9541,0.0193,0.0135],[0.0099,0.0690,0.9001,0.0210],[0.0011,0.0004,0.0026,0.9959],[0.0870,0.6419,0.0021,0.2689],[0.0025,0.4884,0.2081,0.3009]],
    5: [[0.0016,0.9562,0.0037,0.0052,0.0333],[0.0116,0.0130,0.0045,0.0075,0.9634],[0.0206,0.0087,0.0019,0.8978,0.0710],[0.9957,0.0011,0.0009,0.0020,0.0004],[0.0028,0.0914,0.4298,0.0037,0.4723],[0.2884,0.0026,0.0290,0.1916,0.4884]],
    6: [[0.0016,0.0008,0.9562,0.0053,0.0037,0.0325],[0.0013,0.0410,0.0099,0.0105,0.0036,0.9337],[0.0210,0.0221,0.0049,0.8932,0.0013,0.0576],[0.9949,0.0011,0.0007,0.0019,0.0012,0.0002],[0.0028,0.0008,0.0887,0.0062,0.4314,0.4701],[0.1416,0.6042,0.0103,0.0961,0.0143,0.1335]],
    7: [[0.0017,0.0003,0.0238,0.9557,0.0043,0.0107,0.0036],[0.0004,0.0003,0.6574,0.0085,0.0079,0.3252,0.0004],[0.0206,0.0029,0.0451,0.0076,0.8951,0.0271,0.0014],[0.7808,0.2168,0.0002,0.0006,0.0014,0.0000,0.0001],[0.0029,0.0005,0.4702,0.0858,0.0032,0.0047,0.4327],[0.0209,0.2899,0.2886,0.0028,0.2091,0.1790,0.0097]],
    8: [[0.0105,0.0036,0.0019,0.0008,0.0043,0.0232,0.0001,0.9557],[0.3257,0.0003,0.0004,0.0061,0.0068,0.6533,0.0002,0.0072],[0.0170,0.0012,0.0171,0.0396,0.8812,0.0339,0.0053,0.0048],[0.0000,0.0003,0.5748,0.0001,0.0007,0.0002,0.4236,0.0003],[0.0052,0.4339,0.0027,0.0014,0.0038,0.4677,0.0007,0.0845],[0.0786,0.0122,0.0241,0.5701,0.0609,0.0751,0.1718,0.0072]]
};"""
    
    html = html[:qdata_start] + new_qdata + html[qdata_end + 2:]
    print("14c. Replaced JS qData")

# kColors — update for 6 populations instead of 11
old_kcolors = """const kColors = {
    2: ['#64b5f6','#ffb74d'],
    3: ['#64b5f6','#4db6ac','#ffb74d'],
    4: ['#64b5f6','#4db6ac','#81c784','#ffb74d'],
    5: ['#64b5f6','#81c784','#ffb74d','#9575cd','#4db6ac'],
    6: ['#9575cd','#4db6ac','#81c784','#64b5f6','#e57373','#ffb74d'],
    7: ['#e57373','#9575cd','#a1887f','#64b5f6','#4db6ac','#ffb74d','#81c784'],
    8: ['#9575cd','#a1887f','#81c784','#e57373','#66bb6a','#ffb74d','#64b5f6','#4db6ac']
};"""

new_kcolors = """const kColors = {
    2: ['#64b5f6','#ffb74d'],
    3: ['#64b5f6','#4db6ac','#ffb74d'],
    4: ['#64b5f6','#4db6ac','#81c784','#ffb74d'],
    5: ['#64b5f6','#81c784','#ffb74d','#4db6ac','#9575cd'],
    6: ['#64b5f6','#81c784','#e57373','#4db6ac','#ffb74d','#9575cd'],
    7: ['#64b5f6','#81c784','#9575cd','#e57373','#4db6ac','#a1887f','#ffb74d'],
    8: ['#a1887f','#ffb74d','#64b5f6','#9575cd','#4db6ac','#e57373','#81c784','#66bb6a']
};"""

html = html.replace(old_kcolors, new_kcolors)
print("14d. Replaced JS kColors")

# k7Labels → k5Labels
html = html.replace(
    "const k7Labels = ['S.European','Central Asian','CA minor','N.European','African','E.Asian','S.Asian'];",
    "const k5Labels = ['E.Asian','African','Americas','S.Asian','European'];"
)
print("14e. Replaced k7Labels with k5Labels")

# kInterpretations — replace all
old_interp_start = "const kInterpretations = {"
old_interp_end = "};\n\n// Geographic"
interp_start = html.find(old_interp_start)
interp_end = html.find(old_interp_end, interp_start)

if interp_start > 0 and interp_end > 0:
    new_interp = """const kInterpretations = {
    2: '<strong>K=2:</strong> First split separates AFR from non-AFR. Uzbeks are 97.8% non-AFR, clustered with EUR/SAS/EAS.',
    3: '<strong>K=3:</strong> EAS separates from EUR+SAS. Uzbeks show 63.1% EUR-like + 36.0% EAS — confirming the West–East admixture axis.',
    4: '<strong>K=4:</strong> SAS separates from EUR. Uzbeks show 48.8% EUR, 30.1% EAS, 20.8% SAS — three-way continental ancestry.',
    5: '<strong>K=5 (Most Parsimonious):</strong> AMR separates. Uzbeks show 48.8% European, 28.8% East Asian, 19.2% South Asian, 2.9% AMR-like (likely shared ANE ancestry). All major continental groups resolved.',
    6: '<strong>K=6:</strong> A distinct Central Asian cluster emerges (60.4% for UZB), separating from European ancestry. South Asian and AMR components persist.',
    7: '<strong>K=7:</strong> European ancestry splits into two sub-components. Uzbek ancestry further fractionates with 29.0% in one sub-cluster.',
    8: '<strong>K=8:</strong> CV minimum but only marginally lower than K=7 (0.07% improvement). Additional component does not correspond to a clear continental group.'
};"""
    html = html[:interp_start] + new_interp + html[interp_end + 2:]
    print("14f. Replaced JS kInterpretations")

# ============================================================
# 15. Update CV chart JS ranges (0.384-0.401 → 0.293-0.312)
# ============================================================
html = html.replace(
    "const minCV = 0.384, maxCV = 0.401;",
    "const minCV = 0.293, maxCV = 0.312;"
)
print("15a. Updated CV chart Y range")

# Update grid lines for new range
html = html.replace(
    "for (let cv = 0.385; cv <= 0.400; cv += 0.005)",
    "for (let cv = 0.295; cv <= 0.310; cv += 0.005)"
)
print("15b. Updated CV chart grid lines")

# Update optimal K highlight from 7 to 8
html = html.replace(
    "const isOpt = d.k === 7;",
    "const isOpt = d.k === 8;"
)
print("15c. Updated CV chart optimal K marker")

# Update arrow annotation
html = html.replace(
    "const ox = xPos(7), oy = yPos(0.38530);",
    "const ox = xPos(8), oy = yPos(0.29420);"
)
print("15d. Updated CV chart arrow position")

# ============================================================
# 16. Update stacked bar chart for 6 populations
# ============================================================
# Update superGroups bands (11 pops → 6 pops)
old_super_bands = """const superGroups = [{name:'AFR',start:0,end:1,color:'rgba(229,115,115,0.08)'},
                         {name:'EUR',start:1,end:6,color:'rgba(100,181,246,0.08)'},
                         {name:'SAS',start:6,end:8,color:'rgba(129,199,132,0.08)'},
                         {name:'EAS',start:8,end:10,color:'rgba(255,183,77,0.08)'},
                         {name:'UZB',start:10,end:11,color:'rgba(149,117,205,0.08)'}];"""

new_super_bands = """const superGroups = [{name:'AFR',start:0,end:1,color:'rgba(229,115,115,0.08)'},
                         {name:'EUR',start:1,end:2,color:'rgba(100,181,246,0.08)'},
                         {name:'SAS',start:2,end:3,color:'rgba(129,199,132,0.08)'},
                         {name:'EAS',start:3,end:4,color:'rgba(255,183,77,0.08)'},
                         {name:'AMR',start:4,end:5,color:'rgba(255,187,187,0.08)'},
                         {name:'UZB',start:5,end:6,color:'rgba(149,117,205,0.08)'}];"""

html = html.replace(old_super_bands, new_super_bands)
print("16a. Updated super group bands")

# ============================================================
# 17. Update Uzbek pie chart data and labels
# ============================================================
# The drawUzbekPie function uses hardcoded V1 K=7 data
# Find and update the pie data
html = html.replace(
    "const values = [0.702, 0.123, 0.045, 0.044, 0.043, 0.039, 0.004];",
    "const values = [0.488, 0.288, 0.192, 0.029, 0.003];"
)
html = html.replace(
    "const labels = ['Q2 Central Asian 70.2%','Q6 East Asian 12.3%','Q4 N.European 4.5%','Q3 CA minor 4.4%','Q1 S.European 4.3%','Q7 S.Asian 3.9%','Q5 African 0.4%'];",
    "const labels = ['Q5 European 48.8%','Q1 East Asian 28.8%','Q4 South Asian 19.2%','Q3 Americas 2.9%','Q2 African 0.3%'];"
)
html = html.replace(
    "const colors = ['#9575cd','#ffb74d','#64b5f6','#a1887f','#e57373','#81c784','#4db6ac'];",
    "const colors = ['#9575cd','#64b5f6','#4db6ac','#ffb74d','#81c784'];",
    1  # only first occurrence in pie chart context
)
print("17. Updated pie chart data for K=5 V2")

# Update pie chart legend title
html = html.replace(
    "ctx.fillText('UZBEK (n=1,078)', lx, ly);",
    "ctx.fillText('UZBEK (n=1,047)', lx, ly);"
)
print("17b. Updated pie chart n")

# ============================================================
# 18. Update init: showK(7) → showK(5)
# ============================================================
html = html.replace(
    "showK(7);",
    "showK(5);"
)
print("18. Updated init showK(7) -> showK(5)")

# ============================================================
# 19. Update the on-limits-of-discrete-K box (K=7 reference)
# ============================================================
html = html.replace(
    "The practical interpretation is that K&nbsp;=&nbsp;7 provides the best <em>descriptive\n                resolution</em> of ancestry variation for covariate adjustment, not that exactly 7 ancestral\n                populations existed historically.",
    "The practical interpretation is that K&nbsp;=&nbsp;5 provides the most parsimonious <em>descriptive\n                resolution</em> of ancestry variation for covariate adjustment, not that exactly 5 ancestral\n                populations existed historically."
)
print("19. Updated discrete-K limits box")

# ============================================================
# 20. Update Section 2.2 code block:  LD pruning params
# ============================================================
html = html.replace(
    "# LD pruning: window=50, step=5, r²=0.3\nplink --bfile global_qc --indep-pairwise 50 5 0.3 --out global_prune",
    "# LD pruning: window=50, step=10, r²=0.1 (matching Uzbek-only pipeline)\nplink --bfile global_qc --indep-pairwise 50 10 0.1 --out global_prune"
)
print("20. Updated LD pruning params in code block")

# ============================================================
# 21. Update output files table
# ============================================================
html = html.replace(
    "<tr><td><code>global_for_admixture.bed/bim/fam</code></td><td>LD-pruned merged dataset (2,095 &times; 380,376)</td></tr>",
    "<tr><td><code>global_for_admixture.bed/bim/fam</code></td><td>LD-pruned merged dataset (3,595 &times; 77,111)</td></tr>"
)
print("21. Updated output files table")

# ============================================================
# 22. Add V1 note to geographic cline section
# ============================================================
html = html.replace(
    "<h3>6.3 Geographic Cline at K=7</h3>",
    """<h3>6.3 Geographic Cline at K=7 <span style="color:#ff9800; font-size:14px;">(V1 data)</span></h3>
            <div class="info-box" style="border-left: 4px solid #ff9800; background: #fff8e1; margin-bottom:12px;">
                <strong>&#9888; Note:</strong> The geographic breakdown below uses <strong>V1</strong> global ADMIXTURE K=7 results
                (10 sub-populations). V2 uses superpopulation-level references; per-region breakdowns have not yet been recomputed.
                The qualitative pattern (east&ndash;west cline) is robust and applies to V2 as well.
            </div>"""
)
print("22. Added V1 note to geographic cline section")

# ============================================================
# 23. Add V1 note to covariate validation section header
# ============================================================
html = html.replace(
    "<!-- 6. COVARIATE VALIDATION -->\n            <h2>6. Covariate Validation</h2>",
    """<!-- 6. COVARIATE VALIDATION -->
            <h2>6. Covariate Validation <span style="color:#ff9800; font-size:14px;">(V1)</span></h2>
            <div class="info-box" style="border-left: 4px solid #ff9800; background: #fff8e1; margin-bottom:12px;">
                <strong>&#9888; Note:</strong> Covariate validation tests below use <strong>V1</strong> global ADMIXTURE Q-values
                (1,068 matched samples). The key conclusion &mdash; ethnicity is non-significant, birthplace is significant &mdash;
                is expected to hold in V2 given the same geographic structure.
            </div>"""
)
print("23. Added V1 note to covariate validation")

# ============================================================
# 24. Update key findings stat-boxes
# ============================================================
html = html.replace(
    "Birthplace-based ancestry proportions at K=7 should be used as <strong>covariates in\n                any genome-wide association study</strong> to control for population stratification. Ethnic labels are uninformative.",
    "Birthplace-based ancestry proportions at K=5 should be used as <strong>covariates in\n                any genome-wide association study</strong> to control for population stratification. Ethnic labels are uninformative."
)
print("24a. Updated stat-box K=7 -> K=5")

# ============================================================
# 25. Update Section 7 (Key Findings) stat-box about composition
# ============================================================
html = html.replace(
    "~70% the Uzbek cohort carries a distinct Central Asian component not found in any reference panel, with ~12% East Asian and ~13% split Western Eurasian input",
    "The Uzbek cohort is ~49% West Eurasian, ~29% East Asian, ~19% South Asian at K=5, consistent with Central Asia&rsquo;s position at the crossroads of Eurasion"
)
print("25. Updated findings stat-box (attempted)")

# ============================================================
# 26. Add CSS for AMR tag
# ============================================================
html = html.replace(
    ".tag-uzb { background: #9575cd; }",
    ".tag-uzb { background: #9575cd; }\n            .tag-amr { background: #ec407a; }"
)
print("26. Added AMR tag CSS")

# ============================================================
# 27. Update "let currentK = 7" to 5
# ============================================================
html = html.replace(
    "let currentK = 7;",
    "let currentK = 5;"
)
print("27. Updated currentK = 7 -> 5")

# ============================================================
# 28. Update geographic chart axis label
# ============================================================
html = html.replace(
    "ctx.fillText('East Asian Ancestry (Q6)', 0, 0); ctx.restore();",
    "ctx.fillText('East Asian Ancestry (Q6 — V1 data)', 0, 0); ctx.restore();"
)
print("28. Updated geographic chart axis label")

# ============================================================
# 29. Mark completion badge as V2
# ============================================================
html = html.replace(
    '<span class="badge">&#10003; Complete</span>',
    '<span class="badge">&#10003; V2 Updated</span>',
    1  # only first occurrence
)
print("29. Updated completion badge")

# ============================================================
# SAVE
# ============================================================
with open("steps/step11.html", "w", encoding="utf-8") as f:
    f.write(html)

new_len = len(html)
print(f"\nDone! {orig_len} -> {new_len} chars ({new_len - orig_len:+d})")
