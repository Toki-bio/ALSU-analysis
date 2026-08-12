# Step 11 ADMIXTURE Correction — Technical Notes (August 2026)

## Summary

Step 11's global ADMIXTURE analysis had a pending "B2+B4 fixes" re-analysis noted
as in-progress on the site. Investigation found the completed run on disk had
used the **wrong input panel**, not just an incomplete run. This document records
the investigation, the fix, and the verification performed, as the source of
truth for the technical-provenance comment embedded in `steps/step11.html`.

## What was wrong

The completed "fixed" run in `/staging/ALSU-analysis/admixture_analysis/global_admixture/`
used a stale 2,095-sample merge dated Feb 19 2026 — only 10 of the 26 available
1000 Genomes populations (~1,017 reference samples) — instead of the documented,
correct 3,595-sample spring2026 panel (1,047 Uzbek + 2,548 1000G across all 26
populations / 5 superpopulations) used everywhere else on the site.

Found by directly comparing sample/variant counts and file modification dates on
the DRAGEN server against the site's own documented pipeline recipe (`step11.html`'s
own pipeline description), not assumed from log text.

## The fix

Relaunched the B2 (strand-ambiguous SNP removal: A/T, C/G pairs) + B4 (long-range
LD exclusion: MHC 6p21, 8p23 inversion, 17q21 inversion, LCT 2q21, plus 17 other
standard regions) correction starting from the **correct** spring2026 panel.

- Script: `rerun_global_fixed_v2.sh`
- Output: `/staging/ALSU-analysis/admixture_analysis/global_admixture_v2/`
- Pipeline: `global_qc` (60,485 SNPs) → B2 (−4,954 → 55,531) → B4 (−800 → 54,731)
  → LD pruning (−13 → 54,718 final markers), same 3,595 samples throughout.

Every stage was verified against actual file contents/counts on the DRAGEN server
(host `Biotech2024`), not just log output.

Two unplanned server reboots interrupted the run mid-K (once during K=6, twice
total). Root cause was not conclusively identified — the BMC/IPMI hardware event
log (`ipmitool sel`) was silent at both reboot times, ruling out the power-supply
flapping event found separately later (~8 hours after) as the cause. No scheduled
reboot job exists (crontab/systemd timers checked). The run script was made
resumable per-K (skips any K whose log already shows a completed CV error) to
survive further interruptions.

## Verification layers applied

1. **Direct file verification** at every pipeline stage (sample/marker counts
   recomputed from the actual `.bim`/`.fam` files, not read from log text).
2. **Replicate-seed stability check.** Full K=2–8 sweep run 5 times total
   (1 main run + 4 independent replicate seeds). Result:

   | K | Range across 5 seeds |
   |---|---|
   | 2 | 0.00200 (real spread) |
   | 3–7 | 0.00001–0.00003 (essentially identical every run) |
   | 8 | 0.00019 (noticeably less stable than 5–7) |

   K=8's nominally "best" CV score is the least reproducible of the useful range.

3. **sNMF independent cross-check.** A completely different statistical method
   from ADMIXTURE, rerun on the same corrected marker set (K=2–10, 10 replicates
   each). Component correlation vs. ADMIXTURE:

   | K | Mean correlation | Range |
   |---|---|---|
   | 3 | 0.9989 | 0.9987–0.9990 |
   | 5 | 0.9983 | 0.9966–0.9993 |
   | 7 | 0.7075 | 0.2703–0.9973 |

   The two independent methods agree almost perfectly through K=5, then diverge
   sharply at K=7 — evidence the two methods stop describing the same structure
   past 5 groups.

4. **Population-mean sanity check.** K=5 means recomputed directly from the
   corrected Q-matrix and compared to the original (unfixed) run: every value
   changed by under 0.6 percentage points. The corrections were real and worth
   doing, but did not meaningfully change the substantive ancestry conclusions.

## Why K=5 is the working number (not the technically-lowest-scoring K=8)

Three independent lines of evidence converge on 5, not one:

- **Matches the known reference groups.** 5 is the number of external reference
  groups (AFR/EUR/SAS/EAS/AMR) this study is anchored to, and all 5 components
  map cleanly 1:1 onto those groups with no ambiguity (see Section 5 of
  `step11.html`).
- **Most stable across replicate runs** (see table above).
- **Where ADMIXTURE and sNMF still agree** (see table above).

The raw CV-error / cross-entropy scores never produce a clean "elbow" on their
own — this is expected when real population history is a continuous gradient
rather than a small number of discrete founding groups, a limitation explicitly
discussed in Section 3.3 of `step11.html` ("On the limits of discrete-K models").

## A genuine, unresolved finding — flagged, not explained away

At K=6, a component appears that is enriched specifically within the Uzbek
cohort (~55% mean) and barely present (under 5%) in any of the 5 reference
groups — i.e., not a split of an existing continental group, but a distinct
signal. Verified this is real and was **already present in the original
(pre-correction) April 2026 data too** — it was simply never investigated in
detail before, because the site's own analysis stopped at K=5 early on and
never examined component identities beyond that. Documented in Section 5 of
`step11.html` but not further explained; most plausibly a genuine
Central-Asian-specific signal not captured by any of the five external
reference populations. Worth a dedicated follow-up investigation.

## Cross-reference bug found and fixed separately

`step11.html`'s Overview previously misattributed the "Uzbek-only ADMIXTURE /
West-East cline" finding to "Step 10" with a link to `step10.html` — but Step 10
is "Multi-Population PBS Analysis" and contains no such content anywhere on that
page. The actual finding lives in `step11.html`'s own Section 3.3. Fixed to link
to `#uzbek-only-validation` on the same page instead.

## Still outstanding

- The interactive chart's per-K `qData` in Section 4 is real, recalculated data,
  but ADMIXTURE's raw column order per K is **not guaranteed consistent between
  runs** ("label switching") — the Qn numbering should be read as descriptive
  per-K, not as a stable identity across different K values.
- The ~44-sample discrepancy between this site's 2,548-sample 1000G reference
  panel and the standard "2,504 unrelated" Phase 3 release was investigated but
  not conclusively resolved. Ruled out relatedness as the cause — cross-checked
  all 45 extra / 1 missing samples against the official 1000G pedigree file
  (`20130606_g1k.ped`); none are flagged as related to anyone.
