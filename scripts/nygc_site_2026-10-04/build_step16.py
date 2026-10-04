from sg import *

t = read("steps/step16.html")
S_ = R["pbs"]
h = "<h2>3. Tier 1: PBS Candidate Association</h2>"
assert h in t
t = t.replace(h, "<h2>3. Tier 1: PBS Candidate Association <span style=\"color:#991b1b;font-size:.7em\">[WITHDRAWN 2026-10-04]</span></h2>"
              "<div style=\"background:#fff4e5;border-left:4px solid #d97706;padding:12px 16px;margin:12px 0;border-radius:6px\"><strong>Withdrawn.</strong> The eight SNPs tested in this section were artefacts of the old 1000 Genomes GRCh38 reference file; on the NYGC 30x reference no SNP reaches any PBS tier (maximum PBS<sub>UZB</sub> " + f"{S_['max']:.3f}" + ", see <a href=\"step10.html\">step 10</a> and <a href=\"step10_pbs_candidates.html\">PBS candidates</a>). The association results below are kept for the record only. The Tier 2 and Tier 3 analyses and the sensitivity runs do not depend on the PBS list.</div>", 1)
a = "(Uzbeks' nearest 1000G neighbour, F<sub>ST</sub>=0.018) would shift the candidate set."
t = wsub(t, a, "(on NYGC, UZB&ndash;SAS F<sub>ST</sub> = 0.0143 and UZB&ndash;EUR 0.0152, see step 14) would change the tree and the statistic; moot now that no candidates remain.")
write("steps/step16.html", t)
print("step16 ok")
