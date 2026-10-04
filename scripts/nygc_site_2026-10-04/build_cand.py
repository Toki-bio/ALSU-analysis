from sg import *
import re

t = read("steps/step10_pbs_candidates.html")
O = json.load(open(D + "old8.json"))
S_ = R["pbs"]
old_pbs = {"11:207698": 0.38, "5:53879140": 0.53, "3:133749168": 0.49, "10:8512594": 0.32, "11:20665570": 1.14, "12:22967890": 2.99, "12:125520190": 2.69, "12:5664803": 2.67}


def pc(x):
    return "n/a" if x is None else f"{100 * x:.1f}%"


rows = ""
for snp, ch, nyg, uzb_arr, note in O["old"]:
    op = O["old_page"][snp]
    rows += (f'<tr><td class="snp">{snp} {ch}</td><td class="pbs dim">{old_pbs[snp]:.2f}</td>'
             f'<td>UZB {op[0]}% &middot; EUR {op[1]}% &middot; EAS {op[2]}% &middot; SAS {op[3]}% &middot; AFR {op[4]}%</td>'
             f'<td>EUR {pc(nyg["EUR"])} &middot; EAS {pc(nyg["EAS"])} &middot; SAS {pc(nyg["SAS"])} &middot; AFR {pc(nyg["AFR"])}</td>'
             f'<td>{pc(uzb_arr)}</td></tr>\n')

rng = [(min(r["eur"], r["eas"], r["sas"], r["afr"]), max(r["eur"], r["eas"], r["sas"], r["afr"]), r["uzb"]) for r in R["pbs_top15"]]
within15 = sum(1 for lo, hi, u in rng if lo <= u <= hi)
top = ""
for i, r in enumerate(R["pbs_top15"]):
    g = [x for x in (r["genes"] or []) if not x.startswith("ENSG")]
    top += (f'<tr><td>{i + 1}</td><td class="snp">{r["snp"]}</td><td class="pbs">{r["pbs"]:.3f}</td>'
            f'<td>UZB {r["uzb"]:.2f} &middot; EUR {r["eur"]:.2f} &middot; EAS {r["eas"]:.2f} &middot; SAS {r["sas"]:.2f} &middot; AFR {r["afr"]:.2f}</td>'
            f'<td>{", ".join(g) if g else "no gene overlap"}</td></tr>\n')

body = f"""<body>

<p><a href="step10.html" style="color:#555;text-decoration:none;">&larr; Back to Step 10: Multi-Population PBS Analysis</a></p>
<h1>PBS candidates: the earlier Tier 1 list (withdrawn) and the current result</h1>
<p class="sub">NYGC 30x reference &nbsp;&middot;&nbsp; {S_['n']:,} SNPs tested &nbsp;&middot;&nbsp; rebuilt 2026-10-04</p>

<div class="correction-banner" style="background:#fff4e5;border:2px solid #d97706;padding:12px 16px;margin:12px 0;font-size:14px;line-height:1.45;color:#1f2937"><strong>Result.</strong> On the NYGC 30x reference no SNP reaches Tier 1, 2 or 3 and the highest PBS<sub>UZB</sub> is {S_['max']:.3f} (Tier 1 needs 0.3). The eight &ldquo;Tier 1 candidates&rdquo; listed on this page until 2026-10-04 were artefacts of the old 1000G GRCh38 reference file; their claimed biology (founder effect, ancestral retention in Africa, Uzbek-specific selection, the rs56186913 / uterine leiomyoma story) is withdrawn. The full analysis is on <a href="step10.html">step 10</a>.</div>

<div class="stats">
  <div class="stat ok"><div class="n">0</div><div class="l">Tier 1</div></div>
  <div class="stat ok"><div class="n">0</div><div class="l">Tier 2</div></div>
  <div class="stat ok"><div class="n">0</div><div class="l">Tier 3</div></div>
  <div class="stat hi"><div class="n">{S_['max']:.3f}</div><div class="l">max PBS</div></div>
  <div class="stat dim"><div class="n">{S_['n']:,}</div><div class="l">SNPs tested</div></div>
</div>

<h2 style="font-size:1.05em;margin:22px 0 8px;">1. What the eight earlier candidates look like on the NYGC reference</h2>
<p style="font-size:.88em;color:#555;line-height:1.5;margin-bottom:10px">&ldquo;Old file&rdquo; is what this page showed before (frequencies from the 1000G GRCh38 shapeit2 file; as printed then, for the allele after &gt;). &ldquo;NYGC&rdquo; is the frequency of the same allele computed directly from the NYGC 30x VCFs (bcftools, SNV record at that position). &ldquo;Uzbek array&rdquo; is the frequency of the same allele in the raw array data (1,247 samples, PLINK <code>--freq</code>, strand-converted where needed); it is n/a for the four sites that are present only in the imputed data.</p>
<table>
<thead><tr><th>SNP, alleles</th><th class="r">old PBS</th><th>Old file (%)</th><th>NYGC (all four references)</th><th>Uzbek array</th></tr></thead>
<tbody>
{rows}</tbody>
</table>
<p class="note" style="margin-top:10px"><b>Reading.</b> Where an Uzbek array frequency exists it agrees with the old value (5:53879140: 49% then, 48.6% now; 10:8512594: 48% then, 47.1% now), so the Uzbek side was not the problem. The reference side was: in the old file EUR carried 1% at 5:53879140 and 3:133749168, whereas NYGC has 53% and 44%; at the three chr12 sites the old file gave AFR 49% and SAS 22&ndash;24%, whereas NYGC has AFR 0.1&ndash;0.4% and SAS 0.4&ndash;3.5%. For 11:20665570 the old file gave EUR 49%, NYGC 2.5%. On NYGC none of the eight sites is differentiated between the Uzbek cohort and the references: three of the four sites with Uzbek array data lie inside the range spanned by the references, and 12:22967890 is 2.2 percentage points above the highest reference value (4.8% vs 2.6% in EUR). For the other four sites (present only in the imputed data) no Uzbek frequency was recomputed here, so only the reference side is shown. Three sites also have an indel or duplication record at the same position in the NYGC call set (10:8512594 CTG&gt;C, 11:20665570 GAT&gt;G, 11:207698 a DUP record), the kind of site that was excluded from the PBS run.</p>

<h2 style="font-size:1.05em;margin:22px 0 8px;">2. The 15 highest PBS<sub>UZB</sub> values on NYGC (none is a candidate)</h2>
<table>
<thead><tr><th>#</th><th>SNP (GRCh38)</th><th class="r">PBS</th><th>Allele frequencies (Uzbek minor allele)</th><th>Gene overlap</th></tr></thead>
<tbody>
{top}</tbody>
</table>
<p class="note" style="margin-top:10px">&Delta;AF (the smallest difference between the Uzbek frequency and any single reference group) is at most {max(r['daf'] for r in R['pbs_top15']):.3f} among these 15, against a Tier 2 cut-off of 0.3; the Uzbek frequency lies inside the range spanned by the four references at {within15} of the 15 SNPs. Annotation (rsID, consequence, GWAS Catalog traits) is on <a href="step12.html">step 12</a>.</p>

</body>
</html>
"""
i = t.index("<body>")
t = t[:i] + body
t = t.replace("<title>PBS Tier 1 Candidates</title>", "<title>PBS candidates (earlier Tier 1 list withdrawn) - ALSU Pipeline</title>")
t = t.replace("tr.pri{border-radius:6px}", "tr.pri{border-radius:6px}\nh2{font-weight:600;color:#333}\nth.r{text-align:right}\ncode{background:#f3f3f3;padding:1px 4px;border-radius:3px;font-size:.9em}")
write("steps/step10_pbs_candidates.html", t)
print("cand ok")
