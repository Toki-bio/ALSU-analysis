"""
Prepare GWAS phenotype and covariate files from clinical + genotype data.

Outputs (PLINK2-compatible, tab-separated):
  data/gwas_pheno.txt      — FID  IID  RPL   (1=case, 0=control, NA=excluded)
  data/gwas_covar.txt       — FID  IID  AGE  BMI  BIRTHPLACE  PC1..PC10  Q1
  data/gwas_sample_map.txt  — FID  IID  PHENO_ID  CLASS  DETAILS

Uses spring2026 data (FID=0, IID=phenotype_id format).

Usage:
  python scripts/prepare_gwas_covariates.py
"""

import csv, json, sys, os
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ── Paths ──────────────────────────────────────────────────────────
PHENO_CSV  = os.path.join(ROOT, 'GWAS от 27.08 - Sheet1.csv')
PCA_FILE   = os.path.join(ROOT, 'data', 'spring2026_pca.eigenvec')
Q_FILE     = os.path.join(ROOT, 'data', 'spring2026_admix.2.Q')
Q_FAM      = os.path.join(ROOT, 'data', 'spring2026_admix.fam')

OUT_PHENO  = os.path.join(ROOT, 'data', 'gwas_pheno.txt')
OUT_COVAR  = os.path.join(ROOT, 'data', 'gwas_covar.txt')
OUT_MAP    = os.path.join(ROOT, 'data', 'gwas_sample_map.txt')
OUT_STATS  = os.path.join(ROOT, 'data', 'gwas_prep_stats.json')

# ── Classification parameters (same as build_pregnancy_viz.py) ─────
SLOT_START = 28; SLOT_SIZE = 7; MAX_SLOTS = 20
LOSS_CODES = {'2', '3', '6'}

# ── Birthplace mapping (from CSV coding) ──────────────────────────
# 1.x = Uzbekistan regions, 2 = Russia, 3+ = other countries
BIRTHPLACE_MAP = {
    '1.1':  'Tashkent_city',
    '1.2':  'Karakalpakstan',
    '1.3':  'Andijan',
    '1.4':  'Bukhara',
    '1.5':  'Jizzakh',
    '1.6':  'Kashkadarya',
    '1.7':  'Navoi',
    '1.8':  'Namangan',
    '1.9':  'Samarkand',
    '1.10': 'Surkhandarya',
    '1.11': 'Syrdarya',
    '1.12': 'Tashkent_region',
    '1.13': 'Fergana',
    '1.14': 'Khorezm',
    '2':    'Russia',
    '3':    'Other',
}

# Genome length for F_ROH calculation (hg38 autosomal ~2.88 Gb)
GENOME_LENGTH_KB = 2_881_033.286


def load_phenotype_csv():
    """Load clinical CSV → dict keyed by phenotype ID."""
    samples = {}
    with open(PHENO_CSV, encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            sid = row[0].strip()
            if not sid:
                continue

            age_raw = row[7].strip()
            try:
                age = int(age_raw)
            except (ValueError, IndexError):
                age = None

            bmi_raw = row[12].strip()
            try:
                bmi = round(float(bmi_raw), 2)
            except (ValueError, IndexError):
                bmi = None

            bplace_raw = row[8].strip()
            # Normalize: "1.1" → Tashkent_city, etc.
            bplace = BIRTHPLACE_MAP.get(bplace_raw)
            if bplace is None:
                # Try truncating to integer part for broad categories
                try:
                    broad = str(int(float(bplace_raw)))
                    bplace = BIRTHPLACE_MAP.get(broad, 'Unknown')
                except (ValueError, TypeError):
                    bplace = 'Unknown'

            nation_raw = row[13].strip()

            # ── Pregnancy event parsing (identical to build_pregnancy_viz.py) ──
            events = []
            for slot in range(MAX_SLOTS):
                base = SLOT_START + slot * SLOT_SIZE
                if base + 4 >= len(row):
                    break
                outcome = row[base + 3].strip()
                ga_raw = row[base + 4].strip()
                if not outcome or outcome == 'NA':
                    continue
                try:
                    ga = float(ga_raw) if ga_raw else None
                except ValueError:
                    ga = None
                events.append({'outcome': outcome, 'ga': ga})

            loss_strict = sum(1 for e in events if e['outcome'] in LOSS_CODES
                              and e['ga'] is not None and e['ga'] < 20)
            loss_no_ga  = sum(1 for e in events if e['outcome'] in LOSS_CODES
                              and e['ga'] is None)
            loss_any    = sum(1 for e in events if e['outcome'] in LOSS_CODES)
            loss_ge20   = sum(1 for e in events if e['outcome'] in LOSS_CODES
                              and e['ga'] is not None and e['ga'] >= 20)
            live        = sum(1 for e in events if e['outcome'] == '1')
            ectopic     = sum(1 for e in events if e['outcome'] == '5')
            medical     = sum(1 for e in events if e['outcome'] == '7')
            preg_now    = sum(1 for e in events if e['outcome'] == '8')
            n_preg      = len(events)

            # ── Classification ──
            age_ok = age is not None and 18 <= age <= 45
            only_ectopic = ectopic > 0 and loss_any == 0
            only_medical = medical > 0 and loss_any == 0

            if only_ectopic:
                cls = 'EXCLUDE_ECTOPIC'
            elif only_medical and loss_any == 0:
                cls = 'EXCLUDE_MEDICAL'
            elif not age_ok:
                cls = 'EXCLUDE_AGE'
            elif loss_ge20 > 0 and loss_strict < 2:
                cls = 'REVIEW'
            elif loss_strict >= 2:
                cls = 'CASE'
            elif loss_strict + loss_no_ga >= 2 and loss_strict < 2:
                cls = 'CASE_PROBABLE'
            elif loss_any == 0 and live >= 1 and not preg_now:
                cls = 'CONTROL'
            elif loss_any == 0 and live >= 1 and preg_now:
                cls = 'CONTROL_PROBABLE'
            elif loss_any == 1:
                cls = 'GRAY'
            else:
                cls = 'UNCLASSIFIED'

            samples[sid] = {
                'age': age, 'bmi': bmi, 'bplace': bplace,
                'nationality': nation_raw, 'cls': cls,
                'loss_strict': loss_strict, 'live': live, 'n_preg': n_preg,
            }
    return samples


def load_pca():
    """Load PCA eigenvectors → dict keyed by IID."""
    pca = {}
    with open(PCA_FILE) as f:
        for line in f:
            parts = line.strip().split()
            if parts[0].startswith('#') or parts[0] == 'FID':
                continue
            iid = parts[1]
            try:
                pcs = [float(x) for x in parts[2:12]]
            except ValueError:
                continue
            pca[iid] = pcs
    return pca


def load_admixture_q():
    """Load K=2 Q values. Sample order from ADMIXTURE .fam file."""
    sample_order = []
    with open(Q_FAM) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sample_order.append(parts[1])

    q_vals = {}
    with open(Q_FILE) as f:
        for i, line in enumerate(f):
            if i >= len(sample_order):
                break
            parts = line.strip().split()
            q1 = float(parts[0])
            q_vals[sample_order[i]] = q1
    return q_vals





def main():
    print('Loading clinical phenotype data...')
    pheno_data = load_phenotype_csv()
    print(f'  {len(pheno_data)} samples in CSV')

    print('Loading PCA eigenvectors...')
    pca_data = load_pca()
    print(f'  {len(pca_data)} genotyped samples')

    print('Loading ADMIXTURE K=2 Q-values...')
    q_data = load_admixture_q()
    print(f'  {len(q_data)} Q-values')

    # ── Merge genotype + phenotype ──
    matched = 0
    unmatched_geno = []
    stats = {
        'pheno_total': len(pheno_data),
        'geno_total': len(pca_data),
        'matched': 0,
        'case': 0, 'control': 0, 'gray': 0, 'review': 0, 'excluded': 0,
        'missing_bmi': 0, 'missing_age': 0,
    }

    pheno_rows = []
    covar_rows = []
    map_rows = []

    for iid in sorted(pca_data.keys()):
        fid = '0'  # spring2026 format: FID=0, IID=phenotype_id
        pheno_id = iid

        if pheno_id not in pheno_data:
            unmatched_geno.append(iid)
            continue

        matched += 1
        p = pheno_data[pheno_id]
        cls = p['cls']

        # Binary phenotype: CASE=1, CONTROL=0, everything else=NA
        if cls == 'CASE':
            rpl = 1
            stats['case'] += 1
        elif cls == 'CONTROL':
            rpl = 0
            stats['control'] += 1
        elif cls == 'GRAY':
            rpl = 'NA'
            stats['gray'] += 1
        elif cls == 'REVIEW':
            rpl = 'NA'
            stats['review'] += 1
        else:
            rpl = 'NA'
            stats['excluded'] += 1

        pheno_rows.append(f'{fid}\t{iid}\t{rpl}')

        # Covariates
        age = p['age'] if p['age'] is not None else 'NA'
        bmi = p['bmi'] if p['bmi'] is not None else 'NA'
        bplace = p['bplace'] if p['bplace'] != 'Unknown' else 'NA'

        if age == 'NA':
            stats['missing_age'] += 1
        if bmi == 'NA':
            stats['missing_bmi'] += 1

        pcs = pca_data.get(iid, [0]*10)
        q1 = q_data.get(iid, 'NA')

        pc_str = '\t'.join(f'{x:.8f}' for x in pcs)
        covar_rows.append(f'{fid}\t{iid}\t{age}\t{bmi}\t{bplace}\t{pc_str}\t{q1}')

        map_rows.append(f'{fid}\t{iid}\t{pheno_id}\t{cls}\t'
                        f'loss={p["loss_strict"]} live={p["live"]} preg={p["n_preg"]}')

    stats['matched'] = matched
    stats['unmatched'] = len(unmatched_geno)

    # ── Write outputs ──
    with open(OUT_PHENO, 'w') as f:
        f.write('#FID\tIID\tRPL\n')
        f.write('\n'.join(pheno_rows) + '\n')
    print(f'\nWrote {OUT_PHENO} ({len(pheno_rows)} samples)')

    with open(OUT_COVAR, 'w') as f:
        f.write('#FID\tIID\tAGE\tBMI\tBIRTHPLACE\t'
                'PC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\tQ1\n')
        f.write('\n'.join(covar_rows) + '\n')
    print(f'Wrote {OUT_COVAR} ({len(covar_rows)} samples)')

    with open(OUT_MAP, 'w') as f:
        f.write('#FID\tIID\tPHENO_ID\tCLASS\tDETAILS\n')
        f.write('\n'.join(map_rows) + '\n')
    print(f'Wrote {OUT_MAP} ({len(map_rows)} samples)')

    with open(OUT_STATS, 'w') as f:
        json.dump(stats, f, indent=2)

    print(f'\n=== GWAS Preparation Summary ===')
    print(f'Phenotype samples:  {stats["pheno_total"]}')
    print(f'Genotyped samples:  {stats["geno_total"]}')
    print(f'Matched:            {stats["matched"]}')
    print(f'Unmatched genotype: {stats["unmatched"]}')
    print(f'  CASE:     {stats["case"]}')
    print(f'  CONTROL:  {stats["control"]}')
    print(f'  GRAY:     {stats["gray"]}  (excluded from GWAS)')
    print(f'  REVIEW:   {stats["review"]}  (excluded from GWAS)')
    print(f'  Other:    {stats["excluded"]}  (excluded from GWAS)')
    print(f'  GWAS-ready: {stats["case"]} cases + {stats["control"]} controls = {stats["case"]+stats["control"]}')
    print(f'  Missing BMI: {stats["missing_bmi"]}')
    print(f'  Missing age: {stats["missing_age"]}')

    if unmatched_geno:
        print(f'\nUnmatched genotype IDs (first 10): {unmatched_geno[:10]}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
