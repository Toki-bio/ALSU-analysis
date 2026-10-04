from sg import *
import re

t = read("steps/step15.html")
Z = R["roh"]; E = R["eth"]
sp = Z["spearman"]; tt = Z["tertiles"]; uo = Z["uzbek_only"]; tail = Z["tail"]; q = Z["quad"]


def fp(p):
    return "&lt; 1e-300" if p == 0 else f"{p:.1e}"


tert_rows = []
for c, lab in (("EUR", "EUR-like"), ("EAS", "EAS-like"), ("SAS", "SAS-like")):
    r = tt[c]
    tert_rows.append([lab, f"{r['bounds'][0]:.2f} / {r['bounds'][1]:.2f}"] + [f"{r[k]['median']:.4f} (mean {r[k]['mean']:.4f}, n={r[k]['n']})" for k in ("low", "mid", "high")] + [f"{sp[c]['rho']:+.2f} ({fp(sp[c]['p'])})", f"{uo[c]['rho']:+.2f} ({fp(uo[c]['p'])})"])

blk = f"""<h2 id="admix-roh">7. ADMIXTURE &times; ROH Cross-Reference (rebuilt 2026-10-04 on NYGC K=4)</h2>
<p>The earlier version of this section described an expected pattern but showed no result, and it used an Uzbek-only K=2 run. It is replaced by a computed comparison of each person&rsquo;s F<sub>ROH</sub> with their ancestry components from the NYGC-reference ADMIXTURE K=4 run (expanded cohort, seed 11, <a href="step11.html">step 11</a>), joined by sample ID ({Z['n']:,} of {Z['n']:,} cohort samples have both). Check on the join: the F<sub>ROH</sub> computed here has mean {Z['froh_mean']:.4f} and median {Z['froh_median']:.4f}, the values given in section 2, with {Z['n_gt_0625']} samples above 0.0625 as stated there.</p>
{table(["Ancestry component", "Tertile bounds (low/mid cut-points)", "Median F<sub>ROH</sub>, lowest third", "middle third", "highest third", "Spearman &rho; (all 1,256)", "Spearman &rho; (self-reported Uzbeks, n = " + str(uo['n']) + ")"], tert_rows)}
{box("key", f"<strong>Result.</strong> F<sub>ROH</sub> is weakly but consistently associated with ancestry: it is lower in people with more EUR-like ancestry (&rho; = {sp['EUR']['rho']:+.2f}), higher in people with more EAS-like ancestry (&rho; = {sp['EAS']['rho']:+.2f}) and slightly lower with more SAS-like ancestry (&rho; = {sp['SAS']['rho']:+.2f}). The same direction holds among self-reported Uzbeks alone ({uo['EUR']['rho']:+.2f}, {uo['EAS']['rho']:+.2f}, {uo['SAS']['rho']:+.2f}). The effect is small: median F<sub>ROH</sub> moves from {tt['EAS']['low']['median']:.4f} to {tt['EAS']['high']['median']:.4f} between the lowest and highest third of EAS-like ancestry, against a spread (interquartile range 0.0126&ndash;0.0180) many times larger. The consanguineous tail (F<sub>ROH</sub> &gt; 0.0625, n = {Z['n_gt_0625']}) is not clearly enriched for any component; the nominal p = 0.03 for EUR-like (lower in the tail) does not survive correction for three tests (mean EUR-like {tail['EUR']['hi']:.2f} vs {tail['EUR']['lo']:.2f} in the rest, p = {tail['EUR']['p']:.2f}; EAS-like {tail['EAS']['hi']:.2f} vs {tail['EAS']['lo']:.2f}, p = {tail['EAS']['p']:.2f}; SAS-like {tail['SAS']['hi']:.2f} vs {tail['SAS']['lo']:.2f}, p = {tail['SAS']['p']:.2f}).")}
<p><strong>The U-shape that was predicted is not seen.</strong> The earlier text expected higher F<sub>ROH</sub> at both extremes of the ancestry distribution. A quadratic term in a regression of F<sub>ROH</sub> on each component is not significant (t = {q['EUR']['t_quad']:.1f} for EUR-like, {q['EAS']['t_quad']:.1f} for EAS-like, {q['SAS']['t_quad']:.1f} for SAS-like); the relation is monotone and weak.</p>
<p><strong>Caution on interpretation.</strong> This is an association, not an explanation. ROH detection on imputed array data depends on how well each ancestry&rsquo;s haplotypes are imputed, so the EAS-like trend could partly reflect imputation quality or reference-panel coverage rather than mating patterns; this was not tested. Because ancestry and F<sub>ROH</sub> are correlated (|&rho;| up to {max(abs(sp[c]['rho']) for c in sp):.2f}), models that use F<sub>ROH</sub> as a covariate should also include ancestry or PCs.</p>
"""
t = put(t, "step15", blk, replace_range=('<h2 id="admix-roh">', '<h2 id="clinical">'))
# anchor repair: put() removed the old h2 along with the old section body; keep the clinical heading as-is.

# remove unsourced external comparison table in section 8
i = t.index("<h3>F<sub>ROH</sub> vs. published populations</h3>")
j = t.index("</table>", i) + len("</table>")
note = ('<h3>F<sub>ROH</sub> vs. published populations</h3><div class="callout callout-info" style="margin:12px 0;padding:12px 16px">'
        'The comparison table that stood here (UK Biobank, 1000G SAS, Qatar/Saudi, Amish) is removed: its values have no source in the project records, '
        'and one row was internally inconsistent (&ldquo;3&ndash;4&times; lower&rdquo; for 0.008&ndash;0.012 against 0.0147). F<sub>ROH</sub> depends strongly on the '
        'SNP set, the minimum ROH length and the calling parameters, so values from other studies are not comparable without recomputing them with the same pipeline.</div>')
t = t[:i] + note + t[j:]
# soften unsourced external comparisons in the text
t = t.replace("is comparable to outbred European populations\n  (~0.01–0.02) and nearly identical to the original cohort (0.0148).", "is nearly identical to the original cohort (0.0148); no external comparison is made (see section 8).")
t = re.sub(r"is comparable to outbred European populations\s*\(~0\.01.0\.02\) and nearly identical to the original cohort \(0\.0148\)\.", "is nearly identical to the original cohort (0.0148); no external comparison is made (see section 8).", t)
t = re.sub(r"Median F<sub>ROH</sub> = 0\.0147 — comparable to outbred European and South Asian populations, nearly identical to the original cohort \(0\.0148\)", "Median F<sub>ROH</sub> = 0.0147, nearly identical to the original cohort (0.0148)", t)
t = t.replace("<li>Consistent with historical preference for endogamous marriages in Uzbek communities</li>", "<li>Whether the tail reflects endogamous marriage customs was not tested here (no pedigree or marriage data were analysed)</li>")
# individual identifiers on a public page -> ranks
for k, v in {"<strong>14-79m</strong>": "<strong>#1</strong>", "<td>312-AM</td>": "<td>#2</td>", "<td>14-20m</td>": "<td>#3</td>", "<td>03-160</td>": "<td>#4</td>", "<td>02-35</td>": "<td>#5</td>", "<strong>08-227</strong>": "<strong>#1 (lowest)</strong>", "<td>08-416</td>": "<td>#2 (lowest)</td>"}.items():
    t = t.replace(k, v)
t = t.replace("<th>Individual</th>", "<th>Rank (individual codes not shown)</th>")
t = banner(t, "<strong>Updated 2026-10-04:</strong> section 7 (ancestry &times; ROH) was recomputed on the NYGC-reference ADMIXTURE; ROH and IBD themselves use Uzbek data only and are unaffected by the reference change. External comparison values were removed because they were unsourced. See the <a href=\"deprecated_old_runs.html\">list of current and deprecated runs</a>.")
write("steps/step15.html", t)
print("step15 ok")
