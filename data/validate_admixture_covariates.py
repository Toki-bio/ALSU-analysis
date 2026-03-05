#!/usr/bin/env python3
"""
Validate ADMIXTURE ancestry components against self-reported covariates.
- Ethnicity (nationality) vs Q-values
- Birthplace region vs Q-values
- Statistical tests (Kruskal-Wallis, ANOVA)
- Outputs JSON summary for web visualization + TSV tables

Usage:
    python3 validate_admixture_covariates.py [--global]

Without --global: validates Uzbek-only ADMIXTURE (K=2)
With --global:    validates global ADMIXTURE (Uzbek + 1000G refs)
"""

import os, sys, json, csv
import numpy as np
from collections import defaultdict, OrderedDict

# ─── Configuration ───────────────────────────────────────────────────
BASE = "/staging/ALSU-analysis"
ADMIX_DIR = os.path.join(BASE, "admixture_analysis")
GLOBAL_DIR = os.path.join(ADMIX_DIR, "global_admixture")

# Covariate files
BIRTH_LIST   = os.path.join(BASE, "table/birth_list.tsv")
NAT_FILE     = os.path.join(BASE, "table/nationalities/nationalities.txt")
BIRTH_FILE   = os.path.join(BASE, "table/birthblace/birthplaces_all.txt")

# Ethnicity code labels (from nat.last survey)
ETHNICITY_LABELS = {
    "1": "Uzbek",
    "2": "Kazakh",
    "3": "Tajik",
    "4": "Russian",
    "5": "Tatar",
    "6": "Korean",
    "7": "Karakalpak",
    "8": "Other",
    "9": "Mongol",
}

# Birthplace code labels (from survey)
BIRTHPLACE_LABELS = {
    "1.1":  "Tashkent city",
    "1.2":  "Andijan",
    "1.3":  "Bukhara",
    "1.4":  "Fergana",
    "1.5":  "Jizzakh",
    "1.6":  "Kashkadarya",
    "1.7":  "Khorezm",
    "1.8":  "Namangan",
    "1.9":  "Samarkand",
    "1.10": "Surkhandarya",
    "1.11": "Syrdarya",
    "1.12": "Tashkent region",
    "1.13": "Navoi",
    "1.14": "Karakalpakstan",
}


def load_covariate_ids():
    """Load sample IDs from birth_list.tsv (phenotype ID list)."""
    ids = []
    with open(BIRTH_LIST, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                ids.append(parts[0])
    return ids


def load_nationalities():
    """Load nationality code for each sample from nationalities.txt.
    Returns dict: phenotype_id -> nationality_code (string)."""
    nat_map = {}
    pheno_ids = load_covariate_ids()
    with open(NAT_FILE, 'r') as f:
        for i, line in enumerate(f):
            parts = line.strip().split('\t')
            if i < len(pheno_ids) and parts:
                code = parts[-1].strip()
                if code and code != '':
                    nat_map[pheno_ids[i]] = code
    return nat_map


def load_birthplaces():
    """Load birthplace code for each sample from birthplaces_all.txt.
    Returns dict: phenotype_id -> birthplace_code (string, e.g. '1.1')."""
    bp_map = {}
    with open(BIRTH_FILE, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sample_id = parts[0].strip()
                bp_code = parts[1].strip()
                if sample_id and bp_code:
                    bp_map[sample_id] = bp_code
    return bp_map


def load_admixture_q(fam_file, q_file):
    """Load ADMIXTURE Q matrix and corresponding sample IDs from .fam file.
    Returns: list of (genotype_id, q_values) tuples."""
    # Read FAM for IDs
    sample_ids = []
    with open(fam_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sample_ids.append(parts[1])  # IID

    # Read Q matrix
    q_matrix = []
    with open(q_file, 'r') as f:
        for line in f:
            vals = [float(x) for x in line.strip().split()]
            q_matrix.append(vals)

    if len(sample_ids) != len(q_matrix):
        print(f"WARNING: FAM has {len(sample_ids)} samples, Q has {len(q_matrix)} rows")
        n = min(len(sample_ids), len(q_matrix))
        sample_ids = sample_ids[:n]
        q_matrix = q_matrix[:n]

    return list(zip(sample_ids, q_matrix))


def genotype_to_phenotype_id(geno_id):
    """Convert genotype ID (e.g., '2_01-02') to phenotype ID (e.g., '01-02').
    Strip the numeric prefix before the first underscore."""
    if '_' in geno_id:
        return geno_id.split('_', 1)[1]
    return geno_id


def kruskal_wallis_manual(groups):
    """Manual Kruskal-Wallis H test (no scipy dependency).
    groups: list of arrays, each array = values for one group.
    Returns H statistic and approximate p-value (chi-sq with k-1 df).
    """
    # Pool all values and rank
    all_vals = []
    group_labels = []
    for gi, g in enumerate(groups):
        for v in g:
            all_vals.append(v)
            group_labels.append(gi)

    N = len(all_vals)
    if N < 3:
        return 0.0, 1.0

    # Rank
    sorted_indices = np.argsort(all_vals)
    ranks = np.empty(N)
    ranks[sorted_indices] = np.arange(1, N + 1)

    # Handle ties: average ranks
    vals_arr = np.array(all_vals)
    sorted_vals = vals_arr[sorted_indices]
    i = 0
    while i < N:
        j = i
        while j < N and sorted_vals[j] == sorted_vals[i]:
            j += 1
        if j > i + 1:
            avg_rank = (i + 1 + j) / 2.0
            for k in range(i, j):
                ranks[sorted_indices[k]] = avg_rank
        i = j

    # Compute H
    k = len(groups)
    H = 0.0
    for gi, g in enumerate(groups):
        ni = len(g)
        if ni == 0:
            continue
        group_rank_sum = sum(ranks[idx] for idx, gl in enumerate(group_labels) if gl == gi)
        H += (group_rank_sum ** 2) / ni
    H = (12.0 / (N * (N + 1))) * H - 3 * (N + 1)

    # P-value approximation using chi-squared with k-1 df
    df = k - 1
    if df <= 0:
        return H, 1.0
    # Simple chi-squared survival function approximation
    p = chi2_survival(H, df)

    return H, p


def chi2_survival(x, k):
    """Approximate chi-squared survival function P(X > x) for k degrees of freedom.
    Uses the regularized incomplete gamma function approximation."""
    if x <= 0:
        return 1.0
    # Use series expansion of regularized lower incomplete gamma
    a = k / 2.0
    z = x / 2.0
    # Continued fraction / series
    if z < a + 1:
        # Series expansion
        s = 1.0 / a
        term = 1.0 / a
        for n in range(1, 200):
            term *= z / (a + n)
            s += term
            if abs(term) < 1e-12:
                break
        ln_gamma_a = log_gamma(a)
        p_lower = np.exp(a * np.log(z) - z - ln_gamma_a + np.log(s))
        return max(0.0, 1.0 - p_lower)
    else:
        # Continued fraction (Lentz's method)
        ln_gamma_a = log_gamma(a)
        f = 1.0
        b0 = z + 1 - a
        c = 1e30
        d = 1.0 / b0 if b0 != 0 else 1e30
        f = d
        for n in range(1, 200):
            an = n * (a - n)
            bn = z + 2 * n + 1 - a
            d = bn + an * d
            if abs(d) < 1e-30:
                d = 1e-30
            c = bn + an / c
            if abs(c) < 1e-30:
                c = 1e-30
            d = 1.0 / d
            delta = d * c
            f *= delta
            if abs(delta - 1.0) < 1e-12:
                break
        p_upper = np.exp(a * np.log(z) - z - ln_gamma_a + np.log(f))
        return max(0.0, min(1.0, p_upper))


def log_gamma(x):
    """Stirling's approximation for log(Gamma(x))."""
    if x <= 0:
        return 0.0
    # Lanczos approximation
    g = 7
    coefs = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7
    ]
    if x < 0.5:
        return np.log(np.pi / np.sin(np.pi * x)) - log_gamma(1 - x)
    x -= 1
    a = coefs[0]
    t = x + g + 0.5
    for i in range(1, len(coefs)):
        a += coefs[i] / (x + i)
    return 0.5 * np.log(2 * np.pi) + (x + 0.5) * np.log(t) - t + np.log(a)


def compute_effect_size(groups):
    """Compute eta-squared (proportion of variance explained by group membership)."""
    all_vals = np.concatenate(groups)
    grand_mean = np.mean(all_vals)
    ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups)
    ss_total = np.sum((all_vals - grand_mean) ** 2)
    if ss_total == 0:
        return 0.0
    return ss_between / ss_total


def run_validation(mode='uzbek'):
    """Run the covariate validation analysis."""
    print(f"\n{'='*70}")
    print(f"  ADMIXTURE Covariate Validation — {mode.upper()} mode")
    print(f"{'='*70}\n")

    # ── 1. Determine files ──
    if mode == 'global':
        fam_file = os.path.join(GLOBAL_DIR, "global_for_admixture.fam")
        q_base = os.path.join(GLOBAL_DIR, "admix_results", "global_for_admixture")
        pop_labels_file = os.path.join(GLOBAL_DIR, "global_pop_labels.txt")
        output_dir = os.path.join(GLOBAL_DIR, "validation")
    else:
        fam_file = os.path.join(ADMIX_DIR, "UZB_for_admixture.fam")
        q_base = os.path.join(ADMIX_DIR, "UZB_for_admixture")
        pop_labels_file = None
        output_dir = os.path.join(ADMIX_DIR, "validation")

    os.makedirs(output_dir, exist_ok=True)

    # ── 2. Check K=2 Q file exists ──
    q2_file = f"{q_base}.2.Q"
    if not os.path.exists(q2_file):
        print(f"ERROR: Q file not found: {q2_file}")
        print("ADMIXTURE may still be running. Exiting.")
        sys.exit(1)

    if not os.path.exists(fam_file):
        print(f"ERROR: FAM file not found: {fam_file}")
        sys.exit(1)

    # ── 3. Load covariates ──
    print("Loading covariates...")
    nat_map = load_nationalities()
    bp_map = load_birthplaces()
    print(f"  Nationalities loaded: {len(nat_map)} samples")
    print(f"  Birthplaces loaded:   {len(bp_map)} samples")

    # ── 4. Load ADMIXTURE results (K=2 through K=8) ──
    results = {}

    for K in range(2, 9):
        q_file = f"{q_base}.{K}.Q"
        if not os.path.exists(q_file):
            print(f"  K={K}: Q file not found, skipping")
            continue

        print(f"\n{'─'*50}")
        print(f"  Processing K={K}...")

        samples = load_admixture_q(fam_file, q_file)
        print(f"  Loaded {len(samples)} samples with {K} components")

        # ── 5. Map genotype IDs → covariates ──
        matched_nat = []
        matched_bp = []

        for geno_id, q_vals in samples:
            pheno_id = genotype_to_phenotype_id(geno_id)

            nat_code = nat_map.get(pheno_id, None)
            bp_code = bp_map.get(pheno_id, None)

            if nat_code:
                matched_nat.append((geno_id, pheno_id, nat_code, q_vals))
            if bp_code:
                matched_bp.append((geno_id, pheno_id, bp_code, q_vals))

        print(f"  Matched to nationality: {len(matched_nat)} samples")
        print(f"  Matched to birthplace:  {len(matched_bp)} samples")

        if len(matched_nat) == 0:
            print("  WARNING: No nationality matches! Check ID mapping.")
            continue

        # ── 6. Ethnicity analysis ──
        print(f"\n  --- Ethnicity vs Q-values (K={K}) ---")
        eth_groups = defaultdict(list)
        for geno_id, pheno_id, nat_code, q_vals in matched_nat:
            label = ETHNICITY_LABELS.get(nat_code, f"Code_{nat_code}")
            eth_groups[label].append(q_vals)

        # Print mean Q per ethnicity
        eth_summary = {}
        print(f"  {'Ethnicity':<15} {'N':>5}  ", end='')
        for c in range(K):
            print(f"{'Q'+str(c+1):>8}", end='')
        print()

        for eth in sorted(eth_groups.keys(), key=lambda x: -len(eth_groups[x])):
            n = len(eth_groups[eth])
            q_arr = np.array(eth_groups[eth])
            means = q_arr.mean(axis=0)
            stds = q_arr.std(axis=0)
            print(f"  {eth:<15} {n:>5}  ", end='')
            for c in range(K):
                print(f"  {means[c]:.4f}", end='')
            print()

            eth_summary[eth] = {
                'n': n,
                'mean_q': means.tolist(),
                'std_q': stds.tolist(),
                'median_q': np.median(q_arr, axis=0).tolist(),
            }

        # Kruskal-Wallis test for each component
        eth_tests = {}
        for c in range(K):
            component_groups = []
            group_labels_list = []
            for eth in sorted(eth_groups.keys()):
                vals = [q[c] for q in eth_groups[eth]]
                if len(vals) >= 2:
                    component_groups.append(np.array(vals))
                    group_labels_list.append(eth)

            if len(component_groups) >= 2:
                H, p = kruskal_wallis_manual(component_groups)
                eta2 = compute_effect_size(component_groups)
                eth_tests[f'Q{c+1}'] = {
                    'H_statistic': round(H, 4),
                    'p_value': p,
                    'p_value_str': f"{p:.2e}" if p < 0.001 else f"{p:.4f}",
                    'eta_squared': round(eta2, 4),
                    'n_groups': len(component_groups),
                    'significant': p < 0.05,
                }
                sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
                print(f"  Kruskal-Wallis Q{c+1}: H={H:.2f}, p={p:.2e}, eta²={eta2:.4f} {sig}")

        # ── 7. Birthplace analysis ──
        bp_summary = {}
        bp_tests = {}
        if len(matched_bp) > 0:
            print(f"\n  --- Birthplace vs Q-values (K={K}) ---")
            bp_groups = defaultdict(list)
            for geno_id, pheno_id, bp_code, q_vals in matched_bp:
                label = BIRTHPLACE_LABELS.get(bp_code, f"Region_{bp_code}")
                bp_groups[label].append(q_vals)

            print(f"  {'Region':<20} {'N':>5}  ", end='')
            for c in range(K):
                print(f"{'Q'+str(c+1):>8}", end='')
            print()

            for reg in sorted(bp_groups.keys(), key=lambda x: -len(bp_groups[x])):
                n = len(bp_groups[reg])
                if n < 3:
                    continue
                q_arr = np.array(bp_groups[reg])
                means = q_arr.mean(axis=0)
                stds = q_arr.std(axis=0)
                print(f"  {reg:<20} {n:>5}  ", end='')
                for c in range(K):
                    print(f"  {means[c]:.4f}", end='')
                print()

                bp_summary[reg] = {
                    'n': n,
                    'mean_q': means.tolist(),
                    'std_q': stds.tolist(),
                }

            # KW test for each component
            for c in range(K):
                component_groups = []
                for reg in sorted(bp_groups.keys()):
                    vals = [q[c] for q in bp_groups[reg]]
                    if len(vals) >= 3:
                        component_groups.append(np.array(vals))
                if len(component_groups) >= 2:
                    H, p = kruskal_wallis_manual(component_groups)
                    eta2 = compute_effect_size(component_groups)
                    bp_tests[f'Q{c+1}'] = {
                        'H_statistic': round(H, 4),
                        'p_value': p,
                        'p_value_str': f"{p:.2e}" if p < 0.001 else f"{p:.4f}",
                        'eta_squared': round(eta2, 4),
                        'n_groups': len(component_groups),
                        'significant': p < 0.05,
                    }
                    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
                    print(f"  Kruskal-Wallis Q{c+1} by region: H={H:.2f}, p={p:.2e}, eta²={eta2:.4f} {sig}")

        # ── 8. Concordance: identify "outliers" ──
        concordance = {}
        if K == 2:
            q_arr = np.array([q for _, _, _, q in matched_nat])
            print(f"\n  --- K=2 Concordance analysis ---")

            # For K=2, find which component is "Western" vs "Eastern"
            # By convention: check which component is higher for Russians (Western)
            russian_idx = [i for i, (_, _, nat, _) in enumerate(matched_nat) if nat == '4']
            if russian_idx:
                russian_q = q_arr[russian_idx].mean(axis=0)
                west_comp = np.argmax(russian_q)
            else:
                west_comp = 0

            east_comp = 1 - west_comp
            print(f"  Western component: Q{west_comp+1}")
            print(f"  Eastern component: Q{east_comp+1}")

            # Check extreme individuals
            for eth_name, eth_data in eth_groups.items():
                q_west = [q[west_comp] for q in eth_data]
                mean_w = np.mean(q_west)
                if eth_name in ['Russian', 'Tatar']:
                    expected = 'high Western'
                    actual = f"{mean_w:.3f}"
                elif eth_name in ['Korean', 'Kazakh', 'Mongol', 'Karakalpak']:
                    expected = 'high Eastern'
                    actual = f"{1-mean_w:.3f}"
                else:
                    expected = 'mixed'
                    actual = f"W={mean_w:.3f}, E={1-mean_w:.3f}"

                concordance[eth_name] = {
                    'expected': expected,
                    'mean_western': round(mean_w, 4),
                    'mean_eastern': round(1 - mean_w, 4),
                    'n': len(eth_data),
                }
                print(f"  {eth_name:<15}: expected {expected}, actual {actual}")

        # ── 9. Store results for this K ──
        results[f'K{K}'] = {
            'K': K,
            'n_samples_total': len(samples),
            'n_matched_nationality': len(matched_nat),
            'n_matched_birthplace': len(matched_bp),
            'ethnicity': {
                'summary': eth_summary,
                'tests': eth_tests,
            },
            'birthplace': {
                'summary': bp_summary,
                'tests': bp_tests,
            },
        }
        if concordance:
            results[f'K{K}']['concordance'] = concordance

    # ── 10. Write TSV output tables ──
    # Ethnicity × Q table for K=2
    if 'K2' in results:
        k2 = results['K2']
        tsv_file = os.path.join(output_dir, "ethnicity_q_values.tsv")
        with open(tsv_file, 'w') as f:
            f.write("Ethnicity\tN\tQ1_mean\tQ1_sd\tQ2_mean\tQ2_sd\n")
            for eth, data in sorted(k2['ethnicity']['summary'].items(), key=lambda x: -x[1]['n']):
                f.write(f"{eth}\t{data['n']}\t{data['mean_q'][0]:.4f}\t{data['std_q'][0]:.4f}\t{data['mean_q'][1]:.4f}\t{data['std_q'][1]:.4f}\n")
        print(f"\n  Wrote: {tsv_file}")

        # Birthplace × Q table
        if k2['birthplace']['summary']:
            tsv_file = os.path.join(output_dir, "birthplace_q_values.tsv")
            with open(tsv_file, 'w') as f:
                f.write("Region\tN\tQ1_mean\tQ1_sd\tQ2_mean\tQ2_sd\n")
                for reg, data in sorted(k2['birthplace']['summary'].items(), key=lambda x: -x[1]['n']):
                    f.write(f"{reg}\t{data['n']}\t{data['mean_q'][0]:.4f}\t{data['std_q'][0]:.4f}\t{data['mean_q'][1]:.4f}\t{data['std_q'][1]:.4f}\n")
            print(f"  Wrote: {tsv_file}")

    # ── 11. Per-sample table: ID, ethnicity, birthplace, Q1, Q2 ──
    if 'K2' in results:
        q2_file_path = f"{q_base}.2.Q"
        samples_k2 = load_admixture_q(fam_file, q2_file_path)
        sample_file = os.path.join(output_dir, "per_sample_covariates.tsv")
        with open(sample_file, 'w') as f:
            f.write("Geno_ID\tPheno_ID\tEthnicity\tBirthplace\tQ1\tQ2\n")
            for geno_id, q_vals in samples_k2:
                pheno_id = genotype_to_phenotype_id(geno_id)
                nat_code = nat_map.get(pheno_id, '')
                nat_label = ETHNICITY_LABELS.get(nat_code, nat_code)
                bp_code = bp_map.get(pheno_id, '')
                bp_label = BIRTHPLACE_LABELS.get(bp_code, bp_code)
                f.write(f"{geno_id}\t{pheno_id}\t{nat_label}\t{bp_label}\t{q_vals[0]:.6f}\t{q_vals[1]:.6f}\n")
        print(f"  Wrote: {sample_file}")

    # ── 12. JSON output ──
    json_file = os.path.join(output_dir, "validation_results.json")
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"  Wrote: {json_file}")

    print(f"\n{'='*70}")
    print(f"  Validation complete! Results in: {output_dir}")
    print(f"{'='*70}\n")

    return results


if __name__ == '__main__':
    mode = 'global' if '--global' in sys.argv else 'uzbek'
    run_validation(mode)
