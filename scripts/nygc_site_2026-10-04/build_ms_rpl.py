from sg import *
import re

t = read("manuscript_rpl.html")
S_ = R["pbs"]
nbsp = "&nbsp;"
# abstract
t = wsub(t, "(Tier&nbsp;1) targeted test of 8 Uzbek-specific selection candidates from a paired population genomics study (this cohort, paper 1); (Tier&nbsp;2)", "(Tier&nbsp;1, withdrawn 2026-10-04) a targeted test of 8 SNPs from an earlier PBS list, which is no longer supported because those SNPs are artefacts of an old reference file; (Tier&nbsp;2)")
# aims
t = wsub(t, "Test selection candidates from a paired Uzbek population-genomics study (paper 1) for association with RPL.", "<del>Test selection candidates from a paired Uzbek population-genomics study (paper 1) for association with RPL.</del> Withdrawn 2026-10-04: the paired study found no supported selection candidates on the NYGC reference (maximum PBS " + f"{S_['max']:.3f}" + "), so there is nothing to test.")
# covariates: sensitivity pointer
t = wsub(t, "confirms adequate control of population structure.", "confirms adequate control of population structure. A re-check in 2026-10 adding self-reported ethnicity, using PC1&ndash;PC10, and restricting to self-reported Uzbeks did not change the picture and nothing reached genome-wide significance (see <a href=\"steps/step16.html\">step 16</a>). Ancestry covariates come from the Uzbek-only PCA, which does not depend on the reference panel.")
# Tier 1 method paragraph
i = t.index("<strong>Tier 1.</strong>"); j = t.index("</p>", i) + 4
t = t[:i] + ("<strong>Tier 1 (withdrawn 2026-10-04).</strong> The earlier plan was a targeted test of 8 SNPs with PBS &ge; 0.3 from the paired population-genomics study. Those SNPs came from batch-dependent genotypes in an old 1000 Genomes GRCh38 file and do not replicate on the NYGC 30x reference (see <a href=\"steps/step10_pbs_candidates.html\">PBS candidates</a>); on the NYGC reference no SNP reaches any PBS tier. The test was run before this was known; its result is kept below for the record but has no interpretation.") + t[j:]
# 3.2 result paragraph
i = t.index("<h3>3.2. Tier 1"); i = t.index("<p>", i); j = t.index("</p>", i) + 4
t = t[:i] + ("<p>\n    Withdrawn. The eight SNPs tested here (earlier &ldquo;Uzbek-specific selection candidates&rdquo;) are not differentiated between Uzbeks and the references on the NYGC 30x reference, so a test of their association with RPL addresses no hypothesis. For the record, none reached Bonferroni significance (smallest p&nbsp;&gt;&nbsp;0.25); this is not a negative result for a selection&nbsp;&harr;&nbsp;reproductive-fitness hypothesis, because the premise no longer holds.\n</p>") + t[j:]
t = banner(t, "<strong>Updated 2026-10-04:</strong> Tier 1 (the eight PBS &ldquo;candidates&rdquo;) is withdrawn because the candidates were artefacts of the old 1000 Genomes reference file; see <a href=\"steps/step10_pbs_candidates.html\">PBS candidates</a>. Ancestry covariates come from the Uzbek-only PCA and are unaffected by the reference change. The genome-wide results (Tier 3) and the known-variant screen (Tier 2) are not changed by this correction. <a href=\"steps/deprecated_old_runs.html\">Current and deprecated runs</a>.")
write("manuscript_rpl.html", t)
print("ms_rpl ok")
