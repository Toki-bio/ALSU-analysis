from sg import *
import pandas as pd
from scipy.stats import norm

t = read("steps/step10.html")
S_ = R["pbs"]; U = R["ue_dist"]; P = R["pairs"]
pb = pd.read_csv(D + "pbs_all.tsv", sep="\t")
n = len(pb)
zmax = (S_["max"] - S_["mean"]) / S_["sd"]
zexp = float(norm.isf(1.0 / n))
exp_max = S_["mean"] + zexp * S_["sd"]
ge = {x: int((pb.PBS_UZB >= x).sum()) for x in (0.02, 0.03, 0.05, 0.07)}


def gene(g):
    named = [x for x in (g or []) if not x.startswith("ENSG")]
    return ", ".join(named) if named else ("gene without symbol" if g else "no gene overlap")


top15 = [[i + 1, r["snp"], f(r["pbs"], 3), f(r["uzb"], 2), f(r["eur"], 2), f(r["eas"], 2), f(r["sas"], 2), f(r["afr"], 2), f(r["daf"], 3), esc(gene(r["genes"]))] for i, r in enumerate(R["pbs_top15"])]
rng = [(min(r["eur"], r["eas"], r["sas"], r["afr"]), max(r["eur"], r["eas"], r["sas"], r["afr"]), r["uzb"]) for r in R["pbs_top15"]]
within15 = sum(1 for lo, hi, u in rng if lo <= u <= hi)
outm = [max(lo - u, u - hi) for lo, hi, u in rng if not (lo <= u <= hi)]
out15_txt = ("at most " + f"{max(outm):.3f}") if outm else "none"
ctl = R["pbs_control"]
ctl_rows = []
for col, lab in (("PBS_EAS", "East Asian branch"), ("PBS_EUR", "European branch")):
    for r in ctl[col][:5]:
        ctl_rows.append([lab, r["snp"], f(r["pbs"], 3), f(r["uzb"], 2), f(r["eur"], 2), f(r["eas"], 2), esc(gene(r["genes"]))])
pbs_eur_max = max(r["pbs"] for r in ctl["PBS_EUR"]); pbs_eas_max = max(r["pbs"] for r in ctl["PBS_EAS"])

blk = f"""
<h2>1. Overview</h2>
<p>The <strong>Population Branch Statistic (PBS)</strong> measures how much allele frequencies have changed along one population&rsquo;s own branch of a three-population tree. For the Uzbek branch the tree is UZB&ndash;EUR&ndash;EAS. Each pairwise F<sub>ST</sub> is turned into a branch length T = &minus;ln(1 &minus; F<sub>ST</sub>) and
PBS<sub>UZB</sub> = (T<sub>UZB-EUR</sub> + T<sub>UZB-EAS</sub> &minus; T<sub>EUR-EAS</sub>) / 2. A large value means that the Uzbek frequency has moved away from <em>both</em> outgroups, which is expected under local selection, strong drift, or an allele that arose in the lineage.</p>
{box("warn", "<strong>This page was rebuilt on 2026-10-04.</strong> The &ldquo;PBS Tier 1&rdquo; candidates reported earlier (8 SNPs, then 5 in the expanded cohort) came from the old 1000G GRCh38 reference file and do not exist on the NYGC 30x reference. On NYGC no SNP reaches any tier, and the highest PBS is " + f(S_['max'], 3) + ". Conclusion: <strong>no population-specific selection signal is supported by this analysis.</strong> Details of the earlier candidates are on <a href='step10_pbs_candidates.html'>step10_pbs_candidates</a>.")}
{cards([(f"{S_['n']:,}", "SNPs analysed"), (f(S_['max'], 3), "maximum PBS<sub>UZB</sub>"), (f(S_['mean'], 4), "mean PBS<sub>UZB</sub>"), ("0 / 0 / 0", "Tier 1 / 2 / 3 SNPs"), ("0", "near-private alleles")])}
<h3>Tier definitions (unchanged from the earlier version)</h3>
{table(["Tier", "Rule", "SNPs on NYGC"], [["1", "PBS<sub>UZB</sub> &ge; 0.3", S_['t1']], ["2", "smallest allele-frequency difference between UZB and any reference group &ge; 0.3 (&Delta;AF)", f"{S_['t2']} (largest &Delta;AF in the panel: {f(S_['maxdaf'], 3)})"], ["3", "near-private: Uzbek MAF &ge; 5% and all four references &le; 1%", S_['np']]])}

<h2>2. Data and method</h2>
<ul>
<li><strong>Samples and markers:</strong> the same as <a href="step9.html">step 9</a>: UZB 1,256, EUR 633, EAS 585, SAS 601, AFR 893; 82,744 LD-pruned SNPs shared with the NYGC 30x GRCh38 reference, of which {S_['n']:,} remain after dropping the 8,573 sites with an overlapping indel, multiallelic or structural-variant record in the NYGC call set (list on DRAGEN: <code>nygc30x_1256_filt/flagged_overlap.txt</code>).</li>
<li><strong>Frequencies:</strong> PLINK 1.9 <code>--freq --keep-allele-order</code> per group, so the same allele is compared in every population. In the tables the allele shown is the Uzbek minor allele.</li>
<li><strong>F<sub>ST</sub>:</strong> per-SNP Weir &amp; Cockerham for UZB&ndash;EUR, UZB&ndash;EAS and EUR&ndash;EAS (<code>--fst --within</code>); negative values are set to 0 before the transform and F<sub>ST</sub> is capped at 0.999.</li>
<li><strong>PBS and tiers:</strong> script <code>pbs_refqc.py</code>. SAS and AFR are not part of the tree; they enter only the tier 2 and 3 rules.</li>
<li><strong>Run:</strong> <code>/staging/ALSU-analysis/spring2026/full_expanded_cohort/nygc30x_1256_filt/run/</code> on DRAGEN (<code>pbs_all.tsv</code>, <code>pbs_stats.json</code>, per-SNP F<sub>ST</sub> and frequency files); scripts <code>pbs_refqc.py</code>, <code>pbs_refqc_run.sh</code> in <code>pbs_refqc_2026-10/</code> and <code>/staging/tmp/scratch/filt.sh</code>.</li>
</ul>

<h2>3. Results</h2>
{img("nygc_pbs_manhattan.png", "PBS of the Uzbek branch for each SNP. The red line is the Tier 1 threshold (0.3); every SNP lies far below it.")}
{table(["Statistic", "Value"], [["SNPs", f"{S_['n']:,}"], ["Mean", f(S_['mean'], 5)], ["Median", f(S_['median'], 5)], ["Standard deviation", f(S_['sd'], 4)], ["95th / 99th / 99.9th percentile", f"{f(S_['p95'], 4)} / {f(S_['p99'], 4)} / {f(S_['p999'], 4)}"], ["Maximum", f(S_['max'], 5)], ["SNPs with PBS &ge; 0.02 / 0.03 / 0.05 / 0.07", f"{ge[0.02]:,} / {ge[0.03]:,} / {ge[0.05]} / {ge[0.07]}"], ["SNPs with PBS &ge; 0.1, &ge; 0.15, &ge; 0.3", f"{S_['gt01']} / {S_['gt015']} / {S_['gt03']}"]])}
{box("key", f"<strong>How far is the maximum from what chance alone produces?</strong> The highest PBS is {f(S_['max'], 3)}, which is {zmax:.1f} standard deviations above the mean. For {n:,} values drawn from a normal distribution with the same mean and standard deviation the expected maximum is about {zexp:.1f} standard deviations, i.e. {f(exp_max, 3)}. The observed maximum is where the bulk of the distribution would put it; there is no outlier tail. (PBS is not exactly normal, so this is a rough yardstick, not a p-value.)")}
<h3>The 15 highest PBS<sub>UZB</sub> values</h3>
{table(["#", "SNP (GRCh38)", "PBS", "UZB", "EUR", "EAS", "SAS", "AFR", "&Delta;AF", "Gene overlap (Ensembl)"], top15)}
<p>None of these is a candidate. &Delta;AF (the smallest difference between the Uzbek frequency and any single reference group) is at most {f(max(r['daf'] for r in R['pbs_top15']), 3)} among these 15, against a Tier 2 cut-off of 0.3, and the Uzbek frequency lies inside the range spanned by the four references at {within15} of the 15 SNPs (at the remaining ones it is outside by {out15_txt}). Functional annotation of this list is on <a href="step12.html">step 12</a>.</p>

<h2>4. Does the pipeline produce large values when a population is really differentiated?</h2>
<p>As a check on dynamic range, PBS was computed for the other two branches of the same tree from the same files. The East Asian branch reaches {f(pbs_eas_max, 3)} and the European branch {f(pbs_eur_max, 3)} in this panel, against {f(S_['max'], 3)} for the Uzbek branch.</p>
{table(["Branch", "SNP", "PBS", "UZB", "EUR", "EAS", "Gene overlap"], ctl_rows)}
<p>These large values are not selection findings either: they are SNPs on which EUR and EAS have opposite frequencies and Uzbeks sit in between, so both outgroup branches look long (the same SNPs head the F<sub>ST</sub> ranking on <a href="step9.html">step 9</a>). The check shows that the pipeline returns large PBS wherever a lineage-specific frequency difference exists in this panel. It is <em>not</em> a power calculation for selection on the Uzbek lineage, and the statistic by construction cannot flag anything that is simply an intermediate frequency between the two outgroups, which is what admixture produces.</p>

<h2>5. How the earlier &ldquo;candidates&rdquo; arose, in two stages</h2>
{table(["Stage", "What happened", "Result"], [
 ["1. Old reference file (until 2026-10-02)", "In <code>ALL.chr*.shapeit2_integrated_v1a.GRCh38</code> the genotype at a few dozen sites depends on the sequencing batch of the sample, not on population (for example 12:22967890: EUR 97% in that file, about 2.6% in NYGC and in gnomAD).", "Fake PBS up to 3.0; Tier 1 lists of 8 (original) and 5 (expanded) SNPs; &ldquo;chr12 cluster&rdquo;."],
 ["2. First NYGC run (2026-10-02)", "Five new PBS peaks appeared. They were artefacts of my own extraction: single-base records were compared at repeat or indel sites where the NYGC call set has overlapping records (for example CA-repeat deletions at 2:30673774, a 12-bp insertion at 22:48293814); the Uzbek array and 1000G phase 3 agree with each other at these sites.", "Removed by excluding all 8,573 sites with an overlapping indel, multiallelic or SV record within 1 bp."],
 ["3. Current result", "All-site check on the filtered panel: the Uzbek frequency lies outside the range spanned by the four reference groups by more than 0.05 at 66 of 74,171 SNPs, by more than 0.08 at 6, by more than 0.10 at 2 (largest margin 0.11), and by more than 0.15 at none.", "No Tier 1, 2 or 3 SNP; maximum PBS " + f(S_['max'], 3) + "."]])}

<h2>6. Limitations</h2>
<ul>
<li>The panel is an LD-pruned set of about 74,000 array SNPs. It samples the genome sparsely and can miss a selected locus that no panel SNP tags. &ldquo;No signal in this panel&rdquo; is not &ldquo;no selection&rdquo;.</li>
<li>PBS with an admixed focal population is hard to interpret: if the Uzbek ancestry is a mixture of the two outgroups, PBS<sub>UZB</sub> is close to zero or negative almost everywhere (mean {f(S_['mean'], 3)}), and the statistic is least sensitive exactly where selection would matter for the admixed gene pool. Methods that model admixture (for example local-ancestry-based tests) were not run.</li>
<li>Reference individuals were not filtered for relatedness; the effect on per-SNP F<sub>ST</sub> was not tested.</li>
</ul>
<h2>7. Next steps (status)</h2>
<ul>
<li>Functional annotation of the new top list: <a href="step12.html">done (step 12)</a>.</li>
<li>Selection analysis that models admixture: not done; see the list of open items in <a href="next_steps.html">next steps</a>.</li>
</ul>
"""
t = old_wrap(t, "step10", "\n            <h2>1. Overview</h2>", None,
             "ARCHIVED, not current: previous version of this page (August 2026, old 1000G reference; &ldquo;8 / 5 Tier-1 candidates&rdquo;). Withdrawn; do not cite.",
             end_regex=r"\s*</div>\s*</div>\s*</div>\s*<script>")
t = put(t, "step10", blk, before="<!--NYGC:step10:OLDBEGIN-->")
t = t.replace("✓ Expanded cohort — August 2026 (1,256 samples)", "✓ NYGC 30x reference, expanded cohort (1,256 samples) — rebuilt 2026-10-04")
t = t.replace('<span class="step-badge step-badge-old">Spring 2026 (1,047 samples) — April 11, 2026</span>', '')
if '<div class="info-box" style="margin:20px 30px 0;">' in t:
    i = t.index('<div class="info-box" style="margin:20px 30px 0;">'); j = t.index("</div>", i) + 6
    t = t[:i] + t[j:]
t = banner(t, "<strong>Updated 2026-10-04:</strong> rebuilt on the NYGC 30x reference; the earlier &ldquo;Tier 1 candidates&rdquo; are withdrawn. See the <a href=\"deprecated_old_runs.html\">list of current and deprecated runs</a>.")
write("steps/step10.html", t)
print("step10 ok")
