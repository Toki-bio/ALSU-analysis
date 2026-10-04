from sg import *

t = read("steps/robustness_analyses.html")
S_ = R["pbs"]
note = ("<div style=\"background:#fff4e5;border-left:4px solid #d97706;padding:12px 16px;margin:12px 0;border-radius:6px\"><strong>Withdrawn 2026-10-04.</strong> "
        "The PBS candidates analysed here are artefacts of the old 1000 Genomes GRCh38 reference file; on the NYGC 30x reference no SNP reaches any PBS tier (maximum PBS<sub>UZB</sub> " + f"{S_['max']:.3f}" + "; "
        "<a href=\"step10.html\">step 10</a>). A confidence interval that excludes zero describes sampling variance around a value computed from wrong reference genotypes, so it supports nothing about selection. Kept for the record only.</div>")
t = wsub(t, '<h2 id="pbs-boot">2. Bootstrap confidence intervals on PBS candidates</h2>', '<h2 id="pbs-boot">2. Bootstrap confidence intervals on PBS candidates [WITHDRAWN]</h2>' + note)
t = wsub(t, "<h3>PBS candidates &mdash; FROH-adjusted RPL association</h3>", "<h3>PBS candidates &mdash; FROH-adjusted RPL association [WITHDRAWN: candidates are artefacts; see section 2]</h3>") if "<h3>PBS candidates &mdash; FROH-adjusted RPL association</h3>" in t else t
t = banner(t, "<strong>Updated 2026-10-04:</strong> the PBS-candidate parts (section 2 and the FROH-adjusted candidate test in section 1) are withdrawn because the candidates were reference-file artefacts. The FROH-covariate GWAS and the conditional analyses do not depend on PBS and use Uzbek-only PCs, so they are unaffected by the reference change. <a href=\"deprecated_old_runs.html\">Current and deprecated runs</a>.")
write("steps/robustness_analyses.html", t)
print("robust ok")
