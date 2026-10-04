from sg import *

t = read("steps/step13.html")
L = R["ld"]; N = R["ld_n"]
pops = ["UZB", "EUR", "SAS", "EAS", "AFR"]
bins = [(r[0], r[1]) for r in L["UZB"]]


def cell(p, i, adj=False):
    v = L[p][i][2]
    return f"{(v - 1 / N[p]) if adj else v:.4f}"


rows_raw = [[f"{b[0]}&ndash;{b[1]}"] + [cell(p, i) for p in pops] for i, b in enumerate(bins)]
rows_adj = [[f"{b[0]}&ndash;{b[1]}"] + [cell(p, i, True) for p in pops] for i, b in enumerate(bins)]
rows_n = [[f"{b[0]}&ndash;{b[1]}"] + [f"{L[p][i][3]:,}" for p in pops] for i, b in enumerate(bins)]
hd = ["Distance (kb)"] + [f"{p} (n={N[p]:,})" for p in pops]


def cross(p, thr):
    v = [r[2] for r in L[p]]
    for i in range(len(v)):
        if v[i] < thr:
            return f"{bins[i][0]}&ndash;{bins[i][1]} kb"
    return "beyond 2 Mb"


half = [[p, cross(p, 0.2), cross(p, 0.1), cross(p, 0.05), cross(p, 0.01)] for p in pops]
u0 = L["UZB"][0][2]; e0 = L["EUR"][0][2]; s0 = L["SAS"][0][2]; a0 = L["EAS"][0][2]; f0 = L["AFR"][0][2]
ul = L["UZB"][-1][2] - 1 / N["UZB"]; el = L["EUR"][-1][2] - 1 / N["EUR"]; sl = L["SAS"][-1][2] - 1 / N["SAS"]

blk = f"""
<h2>1. Overview</h2>
<p>Linkage disequilibrium (LD) is the non-random association between alleles at nearby variants, measured here as r&sup2; between pairs of SNPs. How fast it decays with distance reflects population history: effective population size, bottlenecks and admixture all leave a mark. This step characterises LD in the Uzbek cohort and compares it with the four NYGC reference groups on exactly the same SNPs.</p>
{box("warn", "<strong>This page was rebuilt on 2026-10-04.</strong> The earlier version of this step asked whether the &ldquo;PBS candidate&rdquo; SNPs of step 10 were independent loci (PLINK clumping, pairwise r&sup2;). Those candidates are withdrawn as artefacts of the old 1000G reference file (see <a href='step10_pbs_candidates.html'>step10_pbs_candidates</a>), so that analysis has no object any more; its text is kept collapsed at the bottom. What remains and is new here is the LD-decay comparison, recomputed on the NYGC reference.")}
{cards([("30,000", "random SNVs (seed 42)"), ("1,256", "Uzbek samples"), ("2,712", "NYGC reference samples"), ("&le; 2 Mb", "pair distance"), ("0.05", "minimum MAF in each group")])}

<h2>2. Data and method</h2>
<ul>
<li><strong>Uzbek data:</strong> imputed, quality-filtered genotypes of the expanded cohort (<code>imputation_results/hq_filtered/UZB_imputed_HQ_qc</code>, 1,256 samples, 5,358,770 variants) on DRAGEN. From its 3,298,098 single-nucleotide variants with Uzbek MAF &ge; 0.05, 30,000 were drawn at random (fixed seed, all chromosomes).</li>
<li><strong>Reference data:</strong> the same 30,000 positions were extracted from the NYGC 30x GRCh38 VCFs (biallelic SNVs only) for EUR (633), EAS (585), SAS (601) and AFR (893). Each group was analysed separately.</li>
<li><strong>r&sup2;:</strong> PLINK 1.9 <code>--r2 --ld-window-kb 2000 --ld-window 99999 --ld-window-r2 0</code>, all SNP pairs on the same chromosome within 2 Mb, with <code>--maf 0.05</code> applied within each group. Pairs were binned by distance and r&sup2; averaged per bin. Because the SNP set is the same, differences between groups are not due to different sites; because MAF is filtered within each group, the pair sets overlap but are not identical.</li>
<li><strong>Sample-size baseline:</strong> r&sup2; computed from n individuals has an expected value of about 1/n even for unlinked SNPs (0.0008 for n = 1,256; 0.0016 for n = 633). The second table subtracts 1/n from every mean so that groups of different size can be compared at long distance.</li>
<li>Script <code>/staging/tmp/ld13.sh</code>; output <code>/staging/tmp/scratch/ld13_full/</code> on DRAGEN (<code>ld_bins.json</code> and per-chromosome r&sup2; files); log <code>/staging/tmp/ld13_full.out</code>. It was toy-tested on chromosome 22 before the full run.</li>
</ul>

<h2>3. Results</h2>
{img("nygc_ld_decay.png", "Mean r&sup2; by distance, five groups, same 30,000 SNPs. Left: raw. Right: after subtracting the 1/n sample-size baseline.")}
<h3>Mean r&sup2; by distance (raw)</h3>
{table(hd, rows_raw)}
<h3>Mean r&sup2; minus 1/n</h3>
{table(hd, rows_adj)}
<h3>Number of SNP pairs per bin</h3>
{table(hd, rows_n)}
<h3>Distance at which mean r&sup2; first falls below a threshold</h3>
{table(["Group", "&lt; 0.2", "&lt; 0.1", "&lt; 0.05", "&lt; 0.01"], half)}

<h2>4. Interpretation</h2>
{box("key", f"<strong>The Uzbek LD curve lies with the Eurasian groups, not with the African group.</strong> At 0&ndash;10 kb the mean r&sup2; is {u0:.3f} (UZB), {e0:.3f} (EUR), {s0:.3f} (SAS), {a0:.3f} (EAS) and {f0:.3f} (AFR). The African group decays fastest at every distance (shorter haplotypes, as expected for a population with a larger long-term effective size). The Uzbek curve is closest to South Asian and European, slightly below East Asian at short range.")}
<ul>
<li><strong>No excess long-range LD.</strong> At 1&ndash;2 Mb the baseline-corrected r&sup2; is {ul:.4f} for UZB, {el:.4f} for EUR and {sl:.4f} for SAS: indistinguishable. Recent or ongoing admixture between differentiated sources would raise long-range LD above that of the parental groups; this is not seen in a 30,000-SNP sample. This does not rule out older admixture (generations ago LD between ancestry blocks has long since decayed) and it is not a date for it.</li>
<li><strong>Short-range LD in the Uzbek group is slightly lower than in EUR, SAS and EAS.</strong> The difference is about 5&ndash;15% (for example {u0:.3f} vs {e0:.3f} at 0&ndash;10 kb). Possible causes that this analysis cannot separate: a larger effective population size or mixed ancestry (more haplotype diversity), and the fact that the Uzbek genotypes are imputed from an array while the reference is sequencing-based (imputation errors lower r&sup2; at short distance). It should not be read as a population-history result by itself.</li>
<li>The earlier version of this page reported r&sup2; = 0.207 in the first 50 kb bin and 0.069 in the next, from 3,000 SNPs with no MAF filter. The numbers are not comparable: rare SNPs have lower r&sup2;, and this analysis keeps only SNPs with MAF &ge; 0.05 in each group (the corresponding 25&ndash;50 kb bin here is {L['UZB'][2][2]:.3f}, the 50&ndash;100 kb bin {L['UZB'][3][2]:.3f}).</li>
</ul>

<h2>5. What happened to the PBS-candidate LD analysis</h2>
{table(["Earlier statement", "Status now"], [["&ldquo;5 PBS candidates, 5 clumps, 0 pairs with r&sup2; &ge; 0.1: all independent loci&rdquo; (expanded cohort); &ldquo;8 candidates&rdquo; (original cohort)", "Withdrawn. The candidate SNPs come from the old 1000G reference file and are not differentiated on NYGC (see <a href='step10_pbs_candidates.html'>step10_pbs_candidates</a>). The finding &ldquo;independent loci&rdquo; was true but uninformative: the sites are on different chromosomes or tens of Mb apart."], ["&ldquo;LD decay: 0.207 &rarr; 0.069, background by 1,500 kb&rdquo;", "Replaced by the table above (different SNP set and MAF filter; see section 4)."], ["&ldquo;Pattern consistent with admixed Central Asian population history&rdquo;", "Not supported by this analysis: no excess long-range LD relative to EUR/SAS. See section 4."]])}

<h2>6. Limitations</h2>
<ul>
<li>Imputed (Uzbek) versus sequenced (reference) genotypes are not equivalent; see section 4.</li>
<li>The baseline correction uses the approximation E[r&sup2;] &asymp; 1/n and ignores relatedness: the NYGC reference includes related individuals (not removed) while the Uzbek set went through the project&rsquo;s own relatedness QC; unequal relatedness can shift LD estimates slightly.</li>
<li>One random sample of 30,000 SNPs was analysed; the sampling error of the bin means was not estimated, but each bin contains thousands of pairs.</li>
</ul>
<h2>7. Next steps</h2>
<ul>
<li><a href="step14.html">Step 14: F<sub>ST</sub> &amp; MDS</a> &mdash; the pairwise F<sub>ST</sub> matrix and its MDS.</li>
<li><a href="step15.html">Step 15: ROH &amp; IBD</a> &mdash; runs of homozygosity and identity by descent.</li>
</ul>
"""
t = old_wrap(t, "step13", "\n            <h2>1. Overview</h2>", None,
             "ARCHIVED, not current: previous version of this page (clumping of PBS candidates, 3,000-SNP LD decay, August 2026). Candidate analysis withdrawn; do not cite.",
             end_regex=r"\s*</div>\s*</div>\s*</div>\s*<script>")
t = put(t, "step13", blk, before="<!--NYGC:step13:OLDBEGIN-->")
t = t.replace("<p>Linkage disequilibrium clumping of PBS candidates and genome-wide LD decay estimation</p>", "<p>Genome-wide LD decay in the Uzbek cohort compared with the NYGC reference populations</p>")
t = t.replace("✓ Expanded cohort — August 2026 (1,256 samples)", "✓ NYGC 30x reference, expanded cohort (1,256 samples) — rebuilt 2026-10-04")
t = t.replace('<span class="step-badge step-badge-old">Original cohort (1,047 samples) — March 2026</span>', '')
if '<div class="info-box" style="margin:20px 30px 0;">' in t:
    i = t.index('<div class="info-box" style="margin:20px 30px 0;">'); j = t.index("</div>", i) + 6
    t = t[:i] + t[j:]
write("steps/step13.html", t)
print("step13 ok")
