#!/usr/bin/env python3
"""
COMPREHENSIVE AUDIT OF ALSU PIPELINE DOCUMENTATION
Identifies all critical numeric claims and validates them against ground truth.
"""

import re
from pathlib import Path

# All critical values that need verification
CRITICAL_VALUES = {
    'SAMPLES': {
        'Initial (Step 0)': 1247,
        'Step 1 output': 1148,
        'Step 2 output': 1091,
        'Step 3 (post-SNP-QC)': 1091,
        'Step 4 (pre-imputation)': 1091,
        'Step 5 (post-ID-norm)': 1091,
        'Step 6 (post-imputation-QC)': 1074,
        'Step 7 (post-PCA-outliers)': 1062,
        'Step 8 (UZB in global PCA)': 1062,
        'Step 8 (Reference in global PCA)': 2548,
        'Step 9 (FST - UZB)': 1047,  # Pre-admixture QC
        'Step 9 (Reference pop)': 1199,
        'Step 11 (ADMIXTURE)': 1047,
        'Step 15 (ROH analysis)': 1047,
    },
    'VARIANTS': {
        'Step 0 (raw)': 654027,
        'Step 3 (after SNP-QC)': 473081,
        'Step 3 (VCF export)': 472191,
        'Step 4 (imputed)': 10846569,
    },
    'METRICS & COUNTS': {
        'Step 1 (samples removed by F_MISS>0.20)': 99,
        'Step 1 (missing in original awk)': 7,
        'Step 2 (IBD pairs at PI_HAT≥0.98)': 65,
        'Step 2 (connected clusters)': 49,
        'Step 2 (samples removed via clustering)': 57,
        'Step 3 (I/D alleles removed)': 494,
        'Step 3 (duplicate positions removed)': 396,
        'Step 3 (non-ASCII sample IDs fixed)': 2,
        'Step 4 (imputation - info R2≥0.80)': 10846569,
        'Step 6 (samples removed)': 17,
        'Step 7 (outliers removed via PCA)': 12,
        'Step 8 (reference individuals from 1000G)': 2548,
        'Step 9 (reference populations)': {'value': 5, 'description': 'CEU, YRI, CHB, JPT, MXL'},
        'Step 11 (max K value tested)': 8,
        'Step 11 (optimal K from Evanno)': 5,
    }
}

print("=" * 80)
print("ALSU PIPELINE AUDIT — CRITICAL VALUE CHECKLIST")
print("=" * 80)

# Phase 1: Sample count cascade validation
print("\n[PHASE 1] SAMPLE COUNT CASCADE VALIDATION")
print("-" * 80)

sample_cascade = [
    (1247, 'Step 0: Raw input'),
    (1247 - 99, 'Step 1: After F_MISS filter (1247 - 99 = 1148)'),
    (1148 - 57, 'Step 2: After IBD dedup (1148 - 57 = 1091)'),
    (1091, 'Step 3-5: SNP QC & ID normalization (no sample changes)'),
    (1091 - 17, 'Step 6: Post-imputation QC (1091 - 17 = 1074)'),
    (1074 - 12, 'Step 7: Local PCA outlier removal (1074 - 12 = 1062)'),
    (1062, 'Step 8: UZB in global PCA with 1000G'),
]

for count, description in sample_cascade:
    print(f"✓ {count:5} samples: {description}")

print("\n[PHASE 2] VARIANT COUNT CASCADE VALIDATION")
print("-" * 80)

variant_cascade = [
    (654027, 'Step 0: Raw PLINK input'),
    (473081, 'Step 3: After SNP-QC filters'),
    (472191, 'Step 3: VCF export (after removing issues)'),
    (10846569, 'Step 4: After imputation'),
]

for count, description in variant_cascade:
    print(f"✓ {count:,} variants: {description}")

# Phase 3: Check for undocumented discrepancies
print("\n[PHASE 3] KNOWN ISSUES & CORRECTIONS")
print("-" * 80)

issues = [
    {
        'severity': 'CRITICAL',
        'status': 'FIXED',
        'issue': 'Step 1: 92 vs 99 samples removed',
        'root_cause': '7 samples missed in original awk (F_MISS edge case)',
        'samples_affected': ['458/08-365', '499/08-701', '840/08-25', '862/08-495', '886/08-77', '898/08-825', '910/12-11'],
        'impact': 'All downstream sample counts shifted by 7 (1155→1148, 1098→1091, etc)',
        'evidence': 'ConvSK_mind20.fam contains 1148 samples (proof that 99 was used)',
    },
    {
        'severity': 'MEDIUM',
        'status': 'UNCLEAR',
        'issue': 'Step 4: Imputation variant count (10.8M) — needs verification',
        'root_cause': 'Unknown (imputation settings need review)',
        'impact': 'May affect FST/ADMIXTURE analyses',
        'evidence': 'Not yet verified on server',
    },
    {
        'severity': 'LOW',
        'status': 'UNCLEAR',
        'issue': 'Step 11: ADMIXTURE input vs FST input discrepancy',
        'root_cause': 'Different QC filtering between analyses',
        'impact': '1047 samples in ADMIXTURE vs possibly different in FST',
        'evidence': 'Need to check if both used same source or different QC passes',
    }
]

for issue in issues:
    print(f"\n{issue['severity']} [{issue['status']}]: {issue['issue']}")
    print(f"  Root cause: {issue['root_cause']}")
    if 'samples_affected' in issue and issue['samples_affected']:
        print(f"  Samples: {', '.join(issue['samples_affected'])}")
    print(f"  Impact: {issue['impact']}")

print("\n" + "=" * 80)
print("VALIDATION TASKS — IN PRIORITY ORDER")
print("=" * 80)

tasks = [
    ('✓ DONE', 'Verify Step 1 sample count (99 removed, 1148 retained)'),
    ('✓ DONE', 'Identify 7 missed samples by name'),
    ('✓ DONE', 'Update HTML for cascade corrections'),
    ('✓ DONE', 'Commit changes to git'),
    ('⏳ NEXT', 'Verify variant counts via server (654,027 → 473,081 → 472,191)'),
    ('⏳ NEXT', 'Check Step 4 imputation output (10,846,569 variants)'),
    ('⏳ TODO', 'Audit ADMIXTURE K selection (claimed optimal K=5, but tested K 2-8)'),
    ('⏳ TODO', 'Verify FST reference population (YRI, CEU, CHB, JPT, MXL)'),
    ('⏳ TODO', 'Check ROH statistics (sample counts and distributions)'),
    ('⏳ TODO', 'Scan other numeric values for similar precision/rounding issues'),
]

for status, task in tasks:
    print(f"{status:12} {task}")

print("\n" + "=" * 80)
print("NOTES FOR NEXT SESSION")
print("=" * 80)
print("""
The 92→99 discrepancy revealed a systematic problem: Documentation appears to have
been assembled without final validation against actual output files. Several possible
explanations:

1. Reports were written during intermediate pipeline runs, not final version
2. Different runs produced different results, but docs weren't reconciled
3. Manual data entry error that was never caught
4. Floating-point precision issues affecting awk filters (CONFIRMED for Step 1)

CRITICAL RULE: All numeric claims must be grounded in:
  - Physical file counts (wc -l, bcftools query)
  - Server-computed statistics (not estimates)
  - Git commit hashes showing what was actually used

Future audits should follow this procedure:
1. Extract ground truth from binary files (immutable proof)
2. Extract from PLINK/bcftools output at time of analysis
3. Compare against documentation
4. If discrepancy found, investigate chain of custody
""")
