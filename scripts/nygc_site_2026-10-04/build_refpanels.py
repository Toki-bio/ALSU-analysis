from sg import *

t = read("steps/reference_panels_expansion.html")
K = R["kaz"]; F = K["fst"]; cv = K["cv"]
order = ["UZB", "KAZ", "EUR", "SAS", "EAS", "AFR"]


def fw(a, b):
    for k in (f"{a}_{b}", f"{b}_{a}"):
        if k in F:
            return F[k]["w"], F[k]["m"]


fst_rows = []
for a, b in [("KAZ", "UZB"), ("KAZ", "EAS"), ("KAZ", "SAS"), ("KAZ", "EUR"), ("KAZ", "AFR"), ("UZB", "EAS"), ("UZB", "SAS"), ("UZB", "EUR"), ("UZB", "AFR"), ("EUR", "SAS"), ("EAS", "SAS"), ("EAS", "EUR"), ("AFR", "SAS"), ("AFR", "EUR"), ("AFR", "EAS")]:
    w, m = fw(a, b)
    fst_rows.append([f"{a} vs {b}", f"{w:.4f}", f"{m:.4f}"])
k4 = K["k4_means"]
k4_rows = [[g, f"{K['n'][g]:,}"] + [f"{k4[g][c]:.3f}" for c in ("EUR", "EAS", "SAS", "AFR")] for g in order]
k5 = K["k5_means"]
k5_rows = [[g] + [f"{x:.3f}" for x in (k5[g][1], k5[g][2], k5[g][3], k5[g][0], k5[g][4])] for g in order]
pm = K["pca_means"]
pc_rows = [[g] + [f"{pm[g][f'PC{i}'][0]:+.4f} &plusmn; {pm[g][f'PC{i}'][1]:.4f}" for i in range(1, 5)] for g in order]
cv_rows = [[k, f"{cv[str(k)]:.5f}"] for k in range(2, 8)]

blk = f"""
<h2 id="kaz-results" style="border-bottom:3px solid #16a34a;">Central Asian reference added: Kazakh GSA samples (results 2026-10-09)</h2>
{box("key", f"<strong>Result.</strong> Adding 224 Kazakh array genotypes (public data, EVA PRJEB89820) to the NYGC + Uzbek panel shows that the Uzbek cohort is far closer to Kazakhs than to any continental reference (F<sub>ST</sub> UZB&ndash;KAZ {fw('UZB','KAZ')[0]:.4f}, against {fw('UZB','SAS')[0]:.4f} to SAS and {fw('UZB','EUR')[0]:.4f} to EUR), and that when Kazakhs are present ADMIXTURE at K=5 finds a fifth component that both populations carry (Uzbek {k5['UZB'][4]*100:.0f}%, Kazakh {k5['KAZ'][4]*100:.0f}%) and no continental reference does (&le; 3%). This is the Central Asian signal that the &ldquo;50% EUR / 30% EAS / 19% SAS&rdquo; description had to squeeze into three continental boxes; the percentages are therefore a description relative to the references, not ancestral proportions.")}

<h3>1. Data and method</h3>
<ul>
<li><strong>Kazakh data:</strong> 224 per-sample VCFs (Illumina GSA-24v2, called with DRAGEN Array; genome build GRCh37; public, CC BY-NC-ND) downloaded from the EBI EVA FTP to DRAGEN (<code>.../reference_panels/kazakh_gsa_prjeb89820/vcf_all/</code>, 3.1 GB). Typed genotypes, not imputed. No metadata (region, clan, sex) were used.</li>
<li><strong>Processing:</strong> merge of the 224 VCFs, conversion to PLINK, UCSC liftOver GRCh37&rarr;GRCh38 (622,863 of 623,053 sites of the first sample lifted), harmonisation of alleles to the existing ADMIXTURE input panel (<code>adm_in</code>, 64,959 LD-pruned SNPs): 35,039 positions shared, 47 strand-flipped; SNP call rate &ge; 95% (420 SNPs dropped), sample call rate &ge; 95% (no sample dropped), no pair with PI_HAT &gt; 0.2 among the 224. Merge with the 3,968 existing samples gives <strong>4,192 samples and 34,619 SNPs</strong>. Script <code>kaz_pipeline.sh</code>, toy-tested on 6 and 40 files before the full run.</li>
<li><strong>Analyses (same SNPs for all):</strong> PLINK 1.9 PCA, Weir&ndash;Cockerham F<sub>ST</sub> for all 15 pairs, ADMIXTURE K=2&ndash;7 (5-fold cross-validation, one seed), unsupervised: the Kazakh samples carry no label. Script <code>kaz_analysis.sh</code>, output <code>.../kazakh_gsa_prjeb89820/merge_full/</code> on DRAGEN.</li>
<li><strong>Consistency check:</strong> on this smaller panel (34,619 SNPs) the Uzbek K=4 means are {k4['UZB']['EUR']:.2f} EUR-like / {k4['UZB']['EAS']:.2f} EAS-like / {k4['UZB']['SAS']:.2f} SAS-like, against 0.50 / 0.30 / 0.19 on the full 64,959-SNP run (<a href="step11.html">step 11</a>), and the reference groups remain nearly pure (EUR {k4['EUR']['EUR']:.2f}, EAS {k4['EAS']['EAS']:.2f}, SAS {k4['SAS']['SAS']:.2f}, AFR {k4['AFR']['AFR']:.2f}).</li>
</ul>

<h3>2. F<sub>ST</sub>, all pairs on the 34,619-SNP merged panel</h3>
{table(["Pair", "Weighted F<sub>ST</sub>", "Mean per-SNP F<sub>ST</sub>"], fst_rows)}
<p>The Kazakh samples are nearest to the Uzbek cohort ({fw('KAZ','UZB')[0]:.4f}) and then to East Asians ({fw('KAZ','EAS')[0]:.4f}), South Asians ({fw('KAZ','SAS')[0]:.4f}) and Europeans ({fw('KAZ','EUR')[0]:.4f}). Values for the other pairs are close to the 82,744-SNP matrix on <a href="step14.html">step 14</a> (for example UZB&ndash;EUR {fw('UZB','EUR')[0]:.4f}, UZB&ndash;SAS {fw('UZB','SAS')[0]:.4f}, UZB&ndash;EAS {fw('UZB','EAS')[0]:.4f}); a pair involving Kazakhs depends on 224 samples and on a typed rather than an imputed genotype source.</p>

<h3>3. PCA</h3>
{img("nygc_kaz_pca.png", "PCA of the merged panel. Uzbek (green) and Kazakh (purple) samples with the NYGC reference groups.")}
{table(["Group", "PC1", "PC2", "PC3", "PC4"], pc_rows)}
<p>On PC2 (the EUR&ndash;EAS axis) Kazakhs lie at {pm['KAZ']['PC2'][0]:+.4f}, between the Uzbek mean ({pm['UZB']['PC2'][0]:+.4f}) and the East Asian mean ({pm['EAS']['PC2'][0]:+.4f}), with much less spread (SD {pm['KAZ']['PC2'][1]:.4f} vs {pm['UZB']['PC2'][1]:.4f}). By nearest-reference label (25 neighbours, PC1&ndash;PC4) {K['knn_KAZ'].get('EAS', 0)} of 224 Kazakh samples fall nearest to the East Asian cluster and {K['knn_KAZ'].get('EUR', 0)} to the European one; this label describes position only, not ancestry. <strong>PC4 again separates the array-based cohorts from the sequenced references</strong>: Kazakh mean {pm['KAZ']['PC4'][0]:+.4f}, Uzbek {pm['UZB']['PC4'][0]:+.4f}, reference groups {min(pm[g]['PC4'][0] for g in ('EUR','SAS','EAS','AFR')):+.4f} to {max(pm[g]['PC4'][0] for g in ('EUR','SAS','EAS','AFR')):+.4f}. This was seen for the Uzbek cohort on <a href="step8.html">step 8</a> and had two possible explanations; the Kazakh samples show the same offset in the same direction, which is expected under either explanation (shared Central Asian ancestry, or a difference between array data and sequencing-based reference), so it does not decide between them.</p>

<h3>4. ADMIXTURE with Kazakhs included</h3>
{table(["K", "Cross-validation error"], cv_rows)}
<p>The error again has no clear minimum: it falls from {cv['2']:.4f} (K=2) to {cv['4']:.4f} (K=4) and then changes by less than 0.001 up to K=7.</p>
<h3>K=4: mean membership (components named by the reference groups)</h3>
{table(["Group", "n", "EUR-like", "EAS-like", "SAS-like", "AFR-like"], k4_rows)}
<p>Kazakhs are on average {k4['KAZ']['EAS']*100:.0f}% EAS-like, {k4['KAZ']['EUR']*100:.0f}% EUR-like and {k4['KAZ']['SAS']*100:.0f}% SAS-like, with a much narrower spread between individuals than Uzbeks (SD of the EAS-like share {K['k4_sd_KAZ']['EAS']:.2f} vs {K['k4_sd_UZB']['EAS']:.2f}).</p>
<h3>K=5: a component that Uzbeks and Kazakhs share</h3>
{table(["Group", "EUR-like", "EAS-like", "SAS-like", "AFR-like", "Central-Asian-like (new)"], k5_rows)}
{img("nygc_kaz_admixture.png", "ADMIXTURE K=4 (top) and K=5 (bottom) on the merged panel; Uzbek and Kazakh samples sorted by their largest component. Colours are arbitrary per K.")}
<p>At K=5 the fifth component is {k5['UZB'][4]*100:.0f}% of the average Uzbek and {k5['KAZ'][4]*100:.0f}% of the average Kazakh genome and at most 3% in any continental reference group. It is probably the same signal as the &ldquo;Uzbek-specific component&rdquo; that appeared at K=7 without Kazakhs (<a href="step11.html">step 11</a>). With it, the average Uzbek is {k5['UZB'][4]*100:.0f}% Central-Asian-like, {k5['UZB'][1]*100:.0f}% EUR-like, {k5['UZB'][3]*100:.0f}% SAS-like and {k5['UZB'][2]*100:.0f}% EAS-like. The method does not say <em>what</em> this component is. Kazakhs and Uzbeks are themselves mixtures of western and eastern Eurasian sources, so it is best read as &ldquo;the part of the genome that is common to these two Central Asian populations and not captured by the continental references&rdquo;, which is what a missing reference group looks like. A technical contribution cannot be excluded (both are GSA-array data, the references are sequencing-based; see PC4).</p>

<h3>5. What this changes</h3>
<ul>
<li>The statement &ldquo;Uzbeks are about 50% European-like, 30% East Asian-like and 19% South Asian-like&rdquo; remains a correct description <em>relative to the four continental references</em> and is reproduced on this panel, but it should not be quoted as ancestral proportions: a Central Asian reference group absorbs about half of the genome at K=5.</li>
<li>&ldquo;Closest to SAS or to EUR&rdquo; (steps 9 and 14) is superseded as a headline: the closest group in this panel is Kazakh, by a factor of about 3 in F<sub>ST</sub>.</li>
<li>The open item &ldquo;add Central Asian reference populations&rdquo; is partly closed: one Central Asian group (Kazakh, 224) is now in; HGDP Central Asian populations from the gnomAD HGDP+1KG release (second section below) have not been added.</li>
</ul>
<h3>6. Limitations</h3>
<ul>
<li>One Kazakh data set, 224 samples, regional and clan composition unknown; one ADMIXTURE seed per K (the other panel used two seeds that agreed to 1e-5).</li>
<li>Typed (Kazakh) versus imputed (Uzbek) genotypes, GSA-24v2 versus GSA-24v3, GRCh37 liftover: all can shift allele frequencies slightly; the analysis uses only SNPs shared by all three sources and 95% call rate.</li>
<li>The panel has 34,619 SNPs (54% of the 64,959 pruned panel), which is enough for continental structure and the Kazakh/Uzbek separation shown but gives less resolution than the full panel.</li>
<li>Unsupervised ADMIXTURE and PCA describe shared structure; they do not date admixture or identify source populations.</li>
</ul>
<h3>7. Files</h3>
<p>DRAGEN: <code>/staging/ALSU-analysis/admixture_analysis/reference_panels/kazakh_gsa_prjeb89820/merge_full/</code> (<code>panel_kaz.*</code>, <code>pca_kaz.*</code>, <code>fst_pairs.tsv</code>, <code>panel_kaz.[2-7].Q</code>, <code>adm_K*.log</code>); scripts <code>kaz_pipeline.sh</code>, <code>kaz_analysis.sh</code> in <code>scripts/nygc_site_2026-10-04/</code>. Data source: EVA PRJEB89820 (Kazakh GSA genotypes; licence CC BY-NC-ND, used here for analysis only, not redistributed).</p>
"""
t = old_wrap(t, "refpanels", '\n    <h2 id="overview">1. Why this branch</h2>', None,
             "ARCHIVED status notes of 2026-09-07/20 (gnomAD HGDP+1KG chr22 pilot that extracted 0 variants, Kazakh sample_1 pilot, merge targets); kept for the record",
             end_regex=r"\s*</div>\s*</div>\s*</body>")
t = put(t, "refpanels", blk, before="<!--NYGC:refpanels:OLDBEGIN-->")
t = t.replace("Started 2026-09-07", "Kazakh reference added 2026-10-09")
t = banner(t, "<strong>Updated 2026-10-09:</strong> the Kazakh GSA samples (224) were merged with the NYGC + Uzbek panel and analysed; the gnomAD HGDP+1KG branch is not done. <a href=\"deprecated_old_runs.html\">Current and deprecated runs</a>.")
write("steps/reference_panels_expansion.html", t)
print("refpanels ok")
