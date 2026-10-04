from sg import *
import re

t = read("steps/step12.html")
A = json.load(open(D + "annot_cache.json"))
S_ = R["pbs"]; O = json.load(open(D + "old8.json"))


def gene(g):
    named = [x for x in (g or []) if not x.startswith("ENSG")]
    return ", ".join(named) if named else ("gene without symbol" if g else "no gene overlap")


def arow(i, r, pbs=False):
    a = A[r["snp"]]
    rs = ", ".join(a["rsids"]) if a["rsids"] else "no rsID in Ensembl"
    cons = ", ".join(sorted(set(c for c in a["cons"] if c))) or "&mdash;"
    if a.get("gwas_status") != "ok":
        gw = "lookup not completed"
    elif a["gwas_n"]:
        gw = f"{a['gwas_n']} association(s): " + esc("; ".join(a["gwas_traits"][:3]))
    else:
        gw = "none listed"
    clin = ", ".join(a["clin"]) if a["clin"] else "none"
    val = f(r["pbs"], 3) if pbs else f(r["fst"], 3)
    return [i + 1, r["snp"], val, esc(rs), esc(gene(r["genes"])), esc(cons), esc(gw), esc(clin)]


top_pbs = [arow(i, r, True) for i, r in enumerate(R["pbs_top15"])]
top_fst = [arow(i, r) for i, r in enumerate(R["ue_top30"])]
n_ok = sum(1 for v in A.values() if v.get("gwas_status") == "ok")
n_hit = sum(1 for v in A.values() if v.get("gwas_n"))
n_total = len(A)
hits = [(s, v) for s, v in A.items() if v.get("gwas_n")]
hit_txt = "; ".join(f"{s} ({', '.join(v['rsids'])}): {', '.join(v['gwas_traits'][:2])}" for s, v in hits) if hits else "none"
n_clin = sum(1 for v in A.values() if v['clin'])
clin_txt = ', '.join(sorted(set(c for v in A.values() for c in v['clin'])))
cons_cnt = {}
for v in A.values():
    for c in set(x for x in v["cons"] if x):
        cons_cnt[c] = cons_cnt.get(c, 0) + 1
cons_rows = [[k, v] for k, v in sorted(cons_cnt.items(), key=lambda kv: -kv[1])]

blk = f"""
<h2 style="border-bottom:3px solid #16a34a;">Annotation of the highest-ranking SNPs on the NYGC reference (rebuilt 2026-10-04)</h2>
<div class="box warn" style="background:#fff3e0;border-left:4px solid #ef6c00;padding:13px 17px;margin:14px 0;border-radius:8px"><strong>What changed.</strong> This step used to annotate the PBS &ldquo;candidates&rdquo; (8 in the original cohort, 5 in the expanded cohort, with the highlighted rs56186913 / uterine leiomyoma lookup). Those SNPs were artefacts of the old 1000G reference file; on the NYGC 30x reference there are no PBS candidates (maximum PBS<sub>UZB</sub> {S_['max']:.3f}; <a href="step10.html">step 10</a>). Annotating nothing would be the honest minimum; instead the same four-source approach is applied to the SNPs that rank highest on the new data, so that the reader can see what they are. None of them is a selection candidate; see the discussion below.</div>

<h3>1. Method</h3>
<ul>
<li><strong>Positions and genes:</strong> Ensembl REST (<code>overlap/region</code>, GRCh38, queried 2026-10-04): overlapping genes, rsIDs and the dbSNP consequence type of the variant at exactly that position.</li>
<li><strong>GWAS Catalog:</strong> EBI GWAS Catalog REST (<code>singleNucleotidePolymorphisms/&lt;rsID&gt;/associations</code>), up to two rsIDs per SNP. The query helper was checked with a positive control (rs1800414, a known pigmentation SNP, must return associations) before use, and every lookup records whether the call succeeded, so that &ldquo;none listed&rdquo; never stands for a failed request. {n_ok} of {n_total} SNP lookups completed.</li>
<li><strong>ClinVar:</strong> the clinical-significance field returned by Ensembl for the same variant.</li>
<li><strong>Not done here:</strong> VEP consequence prediction with SIFT/PolyPhen, and GTEx eQTL lookups (done earlier for the withdrawn candidates). For a list of ancestry-informative SNPs they would add nothing a reader could interpret.</li>
</ul>

<h3>2. The 15 highest PBS<sub>UZB</sub> values</h3>
{table(["#", "SNP (GRCh38)", "PBS", "rsID", "Gene overlap", "Consequence (dbSNP)", "GWAS Catalog", "ClinVar"], top_pbs)}
<h3>3. The 30 highest UZB&ndash;EUR F<sub>ST</sub> values</h3>
{table(["#", "SNP (GRCh38)", "F<sub>ST</sub>", "rsID", "Gene overlap", "Consequence (dbSNP)", "GWAS Catalog", "ClinVar"], top_fst)}
<h3>4. Summary</h3>
<ul>
<li>Of the {n_total} SNPs, {n_hit} have a GWAS Catalog association listed: {esc(hit_txt)}.</li>
<li>Consequence types among the {n_total} SNPs (a SNP can have more than one): {", ".join(f"{k} {v}" for k, v in sorted(cons_cnt.items(), key=lambda kv: -kv[1]))}. No SNP is annotated as missense or loss of function. {n_clin} have a ClinVar-significance value in Ensembl ({clin_txt}); none is pathogenic.</li>
<li><strong>Interpretation.</strong> The top F<sub>ST</sub> SNPs are ancestry-informative markers (Uzbek frequency between EUR and EAS; <a href="step9.html">step 9</a>), and the top PBS SNPs have maximum PBS {S_['max']:.3f}; a gene annotation does not turn either into a selection finding. The mostly intronic and intergenic consequences are what one expects of common neutral variants.</li>
</ul>

<h3>5. What became of the earlier annotation claims</h3>
{table(["Earlier claim", "Status"], [
    ["rs56186913 (11:207698, BET1L/RIC8A 5&prime; UTR) is genome-wide significant for uterine leiomyoma and age at menarche and is &ldquo;the most biologically interesting of the 5 candidates&rdquo;", "Withdrawn as a PBS candidate. The GWAS association of the variant itself is a property of the variant and is not in question, but the premise (UZB MAF 47.8% vs &le; 9.1% in every reference) was wrong: on NYGC the allele is at 59% in EUR, 60% in SAS, 22% in EAS and 27% in AFR (<a href='step10_pbs_candidates.html'>PBS candidates</a>); the Uzbek frequency was not recomputed here because the site is only in the imputed data."],
    ["Annotation of the 8 / 5 candidates: no GWAS Catalog, ClinVar or GTEx hits for the other SNPs", "Not carried over; the SNPs are artefacts."],
    ["&ldquo;Annotate the 490 relaxed candidates&rdquo; (initial run)", "Stale V1 result; not supported by any current data."]])}
<h3>6. Limitations</h3>
<ul>
<li>The lists are a small slice of the ranking, chosen by rank; annotation of a ranked top list is descriptive, not a test of enrichment.</li>
<li>Ensembl and the GWAS Catalog were queried on 2026-10-04; their content changes over time. Results are stored with the query date in the project folder (<code>annot_cache.json</code>, scripts <code>annotate.py</code> and <code>annotate_gwas.py</code>).</li>
<li>The GWAS Catalog lookup used rsIDs found at the exact position; SNPs without an rsID or with several overlapping records may be under-annotated.</li>
</ul>
"""
i = t.index('<h2 style="border-bottom:3px solid #4caf50;">0. Expanded Cohort Results (August 2026)</h2>')
t = old_wrap(t, "step12", t[i:i + 80], None,
             "ARCHIVED, not current: previous version of this page (annotation of the withdrawn PBS candidates, March and August 2026). Do not cite.",
             end_regex=r"\s*</div>\s*</div>\s*</div>\s*<script>")
t = put(t, "step12", blk, before="<!--NYGC:step12:OLDBEGIN-->")
t = t.replace("&#9989; Expanded cohort — August 2026 (5 candidates)", "&#9989; NYGC 30x reference — rebuilt 2026-10-04 (top-ranking SNPs)")
t = t.replace('<span class="step-badge step-badge-old">Original cohort (8 candidates) — March 2026</span>', '')
t = banner(t, "<strong>Updated 2026-10-04:</strong> the PBS candidates that this step used to annotate are withdrawn; the step now annotates the highest-ranking SNPs on the NYGC reference. <a href=\"deprecated_old_runs.html\">Current and deprecated runs</a>.")
write("steps/step12.html", t)
print("step12 ok", n_ok, n_total)
