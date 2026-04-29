# Rescue Tactics for the 69 Identity-Resolution Samples

Date: 2026-04-29

## Bottom Line

Do not full-regenotype all 69 as the first move. The 69-sample list is primarily an identity-resolution panel, not a genotype-quality rescue panel.

Optimal tactic:

1. Run a 96-SNP identity fingerprint panel on all 69 original DNA tubes, keyed by Sample_ID plus Sentrix barcode and position.
2. Compare each fingerprint against all existing GSA genotypes, not only the nominated pair.
3. Full-regenotype only conditional cases: failed-chip samples and original tubes whose fingerprint proves the current GSA genotype belongs to someone else.

This is the highest-value path because most of the 69 have clean genotypes already. Their problem is label-to-phenotype uncertainty, not missing genotype data.

## Evidence Checked

Local evidence:

- `data/investigation_data_v2.json`: full sample verdicts, PI_HAT pairs, identity audit, fingerprint plan.
- `data/fingerprint_samples.txt` and `data/fingerprint_samples_lab.txt`: 69-sample lab list.
- `data/king_data.json`: KING validation of PLINK duplicate pairs.
- `data/scanner_qc_data.json`, `data/idat_controls_v2.txt`: chip and IDAT control evidence.

Server evidence:

- Raw PLINK files exist at `/staging/ALSU-analysis/winter2025/PLINK_301125_0312/`:
  - `ConvSK_raw.bed`, `ConvSK_raw.bim`, `ConvSK_raw.fam`
  - `ConvSK_raw_miss.imiss`
  - `ConvSK_mind20_ibd.genome`
- All 69 samples are present in server `ConvSK_raw.fam` and `ConvSK_raw_miss.imiss`.
- Server raw `F_MISS` values match the local investigation JSON exactly for all 69.
- Server `ConvSK_mind20_ibd.genome` confirms the expected PI_HAT edges among the 69.
- Original merged Illumina sample sheet is present at `/staging/ALSU-analysis/Conversion/1248_merged_sample_sheet_12022025.csv.csv`.
- The indexed GenomeStudio `sample_sheet_indexed.csv` under `duplicate_genotypes` is present but zero bytes, so use the original merged sheet plus barcode/position from the PLINK/investigation data.
- The merged sample sheet contains 54 of the 69 IDs directly. The 15 missing IDs are d/t suffix entries added later to disambiguate duplicate Sample_ID rows; request those by base ID plus Sentrix barcode/position, not by suffix alone.

## What the 69 Samples Are

The list contains:

- 30 `KEEP` samples with unverified identity.
- 33 removed partners from PLINK PI_HAT duplicate clusters.
- 3 KING-only clean checks.
- 3 failed-chip / KING-linked samples on chip `208993030109`.

Pair structure:

- 36 total pairs to resolve.
- 33 PLINK PI_HAT identity pairs.
- 3 KING-only duplicate-like pairs.
- Among the 33 PLINK identity pairs: 12 same-chip and 21 cross-chip.
- KING confirms all 55 unexpected PLINK duplicate pairs in the broader dataset; no PLINK duplicate pair was refuted.

## Tactical Classes

### Tier 1 - Fingerprint First, No Full Array Initially

These are the core identity-resolution samples. They already have usable GSA genotypes, so full regenotyping is wasteful unless fingerprinting proves the current genotype is attached to the wrong label.

KEEP but identity-unverified samples:

```text
01-18, 01-29, 02-104, 02-36, 02-52, 02-59, 02-90, 03-37,
04-14, 04-22, 04-25, 04-45, 06-04, 06-30, 07-04, 07-10,
08-124, 08-160, 08-181, 08-265, 08-436, 08-498, 08-509,
08-541, 08-744, 08-774, 08-799, 08-817, 09-101, 12-05
```

Removed partner samples to fingerprint as the comparator tube:

```text
01-53, 02-29, 02-39, 02-45d, 02-45t, 02-49d, 02-52d,
04-20d, 04-22d, 04-36, 04-40, 04-54d, 04-55d, 06-06d,
06-23d, 06-41d, 07-02d, 07-04d, 07-15, 07-16, 07-17,
07-19, 08-128, 08-179, 08-194, 08-267, 08-45, 08-493,
08-795d, 08-81d, 09-37, 09-76, 12-04
```

Expected outcome: verify up to 30 currently excluded GWAS samples. If confirmed, these move from `KEEP/unverified` to `KEEP/verified` and can be used for phenotype association.

### Tier 2 - KING-Only Clean Checks

These are not part of the PLINK PI_HAT >= 0.98 duplicate set, but KING flags duplicate-like kinship. Fingerprint both tubes first.

```text
04-50, 08-123, 20-08
```

Specific caution:

- `20-08` and `08-123` are both currently kept but have moderate missingness (`F_MISS` about 0.089 and 0.124). If their fingerprints match each other, treat them as the same individual and keep or regenotype only the better-supported identity.
- `04-50` is paired by KING with `07-04d`, which is already a removed d/t partner. Fingerprinting should decide whether this is identity ambiguity or a quality-driven KING-only duplicate signal.

### Tier 3 - Failed-Chip Samples: Full Regenotype Recommended

These three are on chip `208993030109`, part of the failed `208993030xxx` batch. Their current GSA calls should not be treated as a reliable truth set without caution.

```text
08-742, 08-749, 08-754
```

Recommended tactic:

- Include them in the 96-SNP panel if the 69-sample order is already being run.
- Also plan full GSA regenotyping from original tubes, ideally as part of a broader failed-chip rescue for chip `208993030109` or the full failed-chip batch.
- `08-742` has a PLINK PI_HAT edge to `09-101` despite being on the failed chip. Fingerprint both first; if `08-742` tube does not match the existing failed-chip GSA profile, full regenotype the original `08-742` tube to recover the actual individual.
- `08-749` and `08-754` are a same-chip KING-only pair on the failed chip. Because both are low quality, full regenotyping both original tubes is more defensible than trying to adjudicate identity from the degraded GSA alone.

## Pair-Level Decision Logic

For every pair, run the fingerprint on both original tubes and compare against all current GSA genotypes.

If tube A matches the current genotype and tube B does not:

- Verify A.
- B is not recovered by fingerprinting; full-regenotype B only if that phenotype/sample is important.

If tube B matches the current genotype and tube A does not:

- The current genotype belongs to B, not A.
- Relabel or quarantine the GSA genotype accordingly.
- Full-regenotype A if A is needed as an independent individual.

If both tubes match the same current genotype:

- This is a true duplicate/tube duplication or repeated DNA aliquot.
- Keep one genotype only; do not count both as independent individuals.

If neither tube matches:

- Do not use either label for GWAS.
- Search the fingerprint against all 1,247 GSA genotypes for a third-sample match.
- Request the original plate-loading manifest and consider full regenotyping from the physical tubes.

If the fingerprints reveal two distinct real individuals:

- Keep the matched GSA genotype under the correct identity.
- Full-regenotype the other original tube if its individual is absent from the dataset.

## Same-Chip Pairs

Same-chip pairs are enriched for well/sample-sheet mapping problems or physical duplicate loading. Fingerprint both wells/tubes, and request the plate-loading manifest if available.

```text
01-29 / 02-29
02-52 / 01-53
07-10 / 04-54d
08-160 / 07-15
08-436 / 07-17
02-90 / 07-16
08-498 / 08-795d
08-124 / 08-128
08-744 / 08-493
01-18 / 09-76
08-265 / 08-267
12-05 / 12-04
08-754 / 08-749  [KING-only, failed chip]
```

## Cross-Chip Pairs

Cross-chip pairs are more consistent with tube swaps, duplicate aliquots, or label propagation errors across runs. Fingerprinting is the correct first assay.

```text
02-104 / 09-37
02-36 / 07-19
02-59 / 02-45d
02-59 / 02-45t
04-14 / 02-49d
04-22 / 06-23d
04-25 / 07-02d
04-45 / 06-06d
06-04 / 02-52d
06-30 / 02-39
07-04 / 04-20d
08-799 / 06-41d
08-541 / 04-22d
08-509 / 04-36
08-509 / 04-40
08-774 / 04-55d
03-37 / 08-179
03-37 / 08-194
08-817 / 08-81d
08-181 / 08-45
09-101 / 08-742  [partner is failed-chip]
07-04d / 04-50  [KING-only]
20-08 / 08-123  [KING-only]
```

## Lab Order Notes

For every sample, provide:

- FID/IID.
- Original Sample_ID.
- Sentrix barcode.
- Sentrix position.
- Current action/status.
- Pair partner.

Important d/t warning:

The d/t suffixes are not original sample-sheet IDs. They were added during/after GenomeStudio handling to disambiguate duplicate Sample_ID rows. For d/t entries, the lab request must include base Sample_ID plus Sentrix barcode and position. Otherwise the lab may pull the wrong duplicate tube.

## Recommended Assay Strategy

1. Run 96-SNP identity fingerprinting for all 69.
2. For d/t entries, identify tubes by base ID plus barcode/position.
3. Compare each fingerprint to all 1,247 GSA genotypes.
4. Promote verified `KEEP/unverified` samples into GWAS eligibility.
5. Quarantine any sample whose original tube does not match its current GSA genotype.
6. Full-regenotype only:
   - the 3 failed-chip samples (`08-742`, `08-749`, `08-754`), and
   - any original tube whose fingerprint proves that its true individual is missing from the current GSA dataset.

This avoids wasting full-array slots on samples whose genotype is already clean, while preserving the option to recover genuinely missing individuals after the identity evidence is known.
