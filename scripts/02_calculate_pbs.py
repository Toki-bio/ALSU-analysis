#!/usr/bin/env python3
"""
Calculate Population Branch Statistic (PBS), delta-AF, and identify
Uzbek-specific SNPs from multi-population Fst and frequency data.

Usage:
    python3 calculate_pbs.py --outdir /staging/ALSU-analysis/pop_specific_analysis

Inputs (expected in outdir):
    - fst_UZB_vs_EUR.fst  (or symlink to genomewide_uzbek_vs_eur_fst.fst)
    - fst_UZB_vs_EAS.fst
    - fst_UZB_vs_SAS.fst
    - fst_UZB_vs_AFR.fst
    - fst_EUR_vs_EAS.fst
    - freq_UZB.frq
    - freq_EUR.frq
    - freq_EAS.frq
    - freq_SAS.frq
    - freq_AFR.frq

Outputs:
    - pbs_results.tsv          (PBS scores for all SNPs)
    - uzbek_specific_snps.tsv  (filtered Uzbek-specific candidates)
    - delta_af_all.tsv         (allele frequency differences across populations)
    - summary_stats.txt        (human-readable summary)
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import sys


def load_fst(filepath):
    """Load PLINK .fst file → DataFrame with CHR, SNP, POS, FST."""
    df = pd.read_csv(filepath, sep=r'\s+', comment='#')
    # PLINK .fst columns: CHR SNP POS NMISS FST
    df = df.rename(columns=lambda c: c.strip())
    # Clamp negative Fst to 0
    df['FST'] = df['FST'].clip(lower=0)
    return df[['CHR', 'SNP', 'POS', 'FST']].copy()


def load_freq(filepath):
    """Load PLINK .frq file → DataFrame with SNP, MAF, A1, A2."""
    df = pd.read_csv(filepath, sep=r'\s+')
    df = df.rename(columns=lambda c: c.strip())
    # .frq columns: CHR SNP A1 A2 MAF NCHROBS
    return df[['CHR', 'SNP', 'A1', 'A2', 'MAF', 'NCHROBS']].copy()


def fst_to_time(fst):
    """Convert Fst to divergence time: T = -log(1 - Fst).
    Clamp Fst to [0, 0.999] to avoid log(0)."""
    fst = np.clip(fst, 0, 0.999)
    return -np.log(1 - fst)


def calculate_pbs(t_ab, t_ac, t_bc):
    """
    Population Branch Statistic for population A (target).
    
    PBS_A = (T_AB + T_AC - T_BC) / 2
    
    Where:
        A = target population (UZB)
        B = outgroup 1 (EUR)
        C = outgroup 2 (EAS)
        T_XY = -log(1 - Fst_XY)
    
    Positive PBS → branch lengthening on UZB lineage → possible selection
    """
    return (t_ab + t_ac - t_bc) / 2


def main():
    parser = argparse.ArgumentParser(description='PBS and Uzbek-specific SNP analysis')
    parser.add_argument('--outdir', required=True, help='Directory with Fst/freq files')
    parser.add_argument('--fst-eur', default=None, 
                        help='Path to UZB-vs-EUR Fst file (if not in outdir)')
    parser.add_argument('--pbs-threshold', type=float, default=0.3,
                        help='PBS threshold for significance (default: 0.3)')
    parser.add_argument('--delta-af-threshold', type=float, default=0.3,
                        help='Min delta-AF vs ALL populations for Uzbek-specific (default: 0.3)')
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    
    print("=" * 70)
    print("PBS & Uzbek-Specific SNP Analysis")
    print("=" * 70)
    
    # ---- Load Fst data ----
    print("\n[1] Loading Fst data...")
    
    fst_eur_path = args.fst_eur or outdir / 'fst_UZB_vs_EUR.fst'
    fst_uzb_eur = load_fst(fst_eur_path)
    fst_uzb_eas = load_fst(outdir / 'fst_UZB_vs_EAS.fst')
    fst_uzb_sas = load_fst(outdir / 'fst_UZB_vs_SAS.fst')
    fst_uzb_afr = load_fst(outdir / 'fst_UZB_vs_AFR.fst')
    fst_eur_eas = load_fst(outdir / 'fst_EUR_vs_EAS.fst')
    
    print(f"  UZB-EUR: {len(fst_uzb_eur):,} SNPs")
    print(f"  UZB-EAS: {len(fst_uzb_eas):,} SNPs")
    print(f"  UZB-SAS: {len(fst_uzb_sas):,} SNPs")
    print(f"  UZB-AFR: {len(fst_uzb_afr):,} SNPs")
    print(f"  EUR-EAS: {len(fst_eur_eas):,} SNPs")
    
    # ---- Load frequency data ----
    print("\n[2] Loading allele frequencies...")
    
    freq_uzb = load_freq(outdir / 'freq_UZB.frq')
    freq_eur = load_freq(outdir / 'freq_EUR.frq')
    freq_eas = load_freq(outdir / 'freq_EAS.frq')
    freq_sas = load_freq(outdir / 'freq_SAS.frq')
    freq_afr = load_freq(outdir / 'freq_AFR.frq')
    
    print(f"  UZB: {len(freq_uzb):,} SNPs")
    print(f"  EUR: {len(freq_eur):,} SNPs")
    print(f"  EAS: {len(freq_eas):,} SNPs")
    print(f"  SAS: {len(freq_sas):,} SNPs")
    print(f"  AFR: {len(freq_afr):,} SNPs")
    
    # ---- Merge all frequency data by SNP ----
    print("\n[3] Merging frequency data across populations...")
    
    # Start with UZB frequencies
    merged = freq_uzb[['SNP', 'CHR', 'A1', 'A2', 'MAF']].rename(
        columns={'MAF': 'MAF_UZB'}
    )
    
    # Merge with other populations (handle allele flips)
    for pop_name, pop_freq in [('EUR', freq_eur), ('EAS', freq_eas), 
                                 ('SAS', freq_sas), ('AFR', freq_afr)]:
        pop_sub = pop_freq[['SNP', 'A1', 'A2', 'MAF']].rename(
            columns={'MAF': f'MAF_{pop_name}', 'A1': f'A1_{pop_name}', 'A2': f'A2_{pop_name}'}
        )
        merged = merged.merge(pop_sub, on='SNP', how='inner')
        
        # Handle allele orientation: if ref alleles are flipped, invert MAF
        flipped = merged['A1'] != merged[f'A1_{pop_name}']
        same_as_a2 = merged['A1'] == merged[f'A2_{pop_name}']
        flip_mask = flipped & same_as_a2
        merged.loc[flip_mask, f'MAF_{pop_name}'] = 1 - merged.loc[flip_mask, f'MAF_{pop_name}']
        
        n_flipped = flip_mask.sum()
        if n_flipped > 0:
            print(f"  {pop_name}: flipped {n_flipped:,} alleles to match UZB orientation")
        
        # Drop temporary allele columns
        merged = merged.drop(columns=[f'A1_{pop_name}', f'A2_{pop_name}'])
    
    print(f"  Merged dataset: {len(merged):,} SNPs with frequencies in all 5 populations")
    
    # ---- Compute delta-AF (allele frequency differences) ----
    print("\n[4] Computing delta-AF...")
    
    merged['dAF_EUR'] = (merged['MAF_UZB'] - merged['MAF_EUR']).abs()
    merged['dAF_EAS'] = (merged['MAF_UZB'] - merged['MAF_EAS']).abs()
    merged['dAF_SAS'] = (merged['MAF_UZB'] - merged['MAF_SAS']).abs()
    merged['dAF_AFR'] = (merged['MAF_UZB'] - merged['MAF_AFR']).abs()
    
    # Minimum delta-AF across all comparisons (high = different from everyone)
    merged['dAF_min'] = merged[['dAF_EUR', 'dAF_EAS', 'dAF_SAS', 'dAF_AFR']].min(axis=1)
    # Maximum delta-AF (most differentiated comparison)
    merged['dAF_max'] = merged[['dAF_EUR', 'dAF_EAS', 'dAF_SAS', 'dAF_AFR']].max(axis=1)
    # Mean delta-AF
    merged['dAF_mean'] = merged[['dAF_EUR', 'dAF_EAS', 'dAF_SAS', 'dAF_AFR']].mean(axis=1)
    
    print(f"  Mean dAF(UZB-EUR): {merged['dAF_EUR'].mean():.4f}")
    print(f"  Mean dAF(UZB-EAS): {merged['dAF_EAS'].mean():.4f}")
    print(f"  Mean dAF(UZB-SAS): {merged['dAF_SAS'].mean():.4f}")
    print(f"  Mean dAF(UZB-AFR): {merged['dAF_AFR'].mean():.4f}")
    
    # ---- Calculate PBS ----
    print("\n[5] Calculating Population Branch Statistic (PBS)...")
    
    # Merge Fst values by SNP
    # UZB-EUR-EAS triangle
    fst_merge = fst_uzb_eur[['SNP', 'FST']].rename(columns={'FST': 'FST_UZB_EUR'})
    fst_merge = fst_merge.merge(
        fst_uzb_eas[['SNP', 'FST']].rename(columns={'FST': 'FST_UZB_EAS'}),
        on='SNP', how='inner'
    )
    fst_merge = fst_merge.merge(
        fst_eur_eas[['SNP', 'FST']].rename(columns={'FST': 'FST_EUR_EAS'}),
        on='SNP', how='inner'
    )
    
    # Also add UZB-SAS and UZB-AFR Fst
    fst_merge = fst_merge.merge(
        fst_uzb_sas[['SNP', 'FST']].rename(columns={'FST': 'FST_UZB_SAS'}),
        on='SNP', how='left'
    )
    fst_merge = fst_merge.merge(
        fst_uzb_afr[['SNP', 'FST']].rename(columns={'FST': 'FST_UZB_AFR'}),
        on='SNP', how='left'
    )
    
    print(f"  SNPs with all three pairwise Fst (UZB-EUR-EAS): {len(fst_merge):,}")
    
    # Compute divergence times
    fst_merge['T_UZB_EUR'] = fst_to_time(fst_merge['FST_UZB_EUR'])
    fst_merge['T_UZB_EAS'] = fst_to_time(fst_merge['FST_UZB_EAS'])
    fst_merge['T_EUR_EAS'] = fst_to_time(fst_merge['FST_EUR_EAS'])
    
    # PBS for UZB branch
    fst_merge['PBS_UZB'] = calculate_pbs(
        fst_merge['T_UZB_EUR'], fst_merge['T_UZB_EAS'], fst_merge['T_EUR_EAS']
    )
    
    # Also compute PBS for EUR and EAS branches (for comparison)
    fst_merge['PBS_EUR'] = calculate_pbs(
        fst_merge['T_UZB_EUR'], fst_merge['T_EUR_EAS'], fst_merge['T_UZB_EAS']
    )
    fst_merge['PBS_EAS'] = calculate_pbs(
        fst_merge['T_UZB_EAS'], fst_merge['T_EUR_EAS'], fst_merge['T_UZB_EUR']
    )
    
    print(f"  Mean PBS_UZB: {fst_merge['PBS_UZB'].mean():.6f}")
    print(f"  Median PBS_UZB: {fst_merge['PBS_UZB'].median():.6f}")
    print(f"  Max PBS_UZB: {fst_merge['PBS_UZB'].max():.4f}")
    
    # PBS percentiles
    for pct in [95, 99, 99.5, 99.9]:
        val = fst_merge['PBS_UZB'].quantile(pct / 100)
        print(f"  PBS_UZB {pct}th percentile: {val:.4f}")
    
    # ---- Combine all data ----
    print("\n[6] Combining all metrics...")
    
    # Merge PBS with frequency data
    combined = merged.merge(
        fst_merge[['SNP', 'FST_UZB_EUR', 'FST_UZB_EAS', 'FST_UZB_SAS', 
                    'FST_UZB_AFR', 'FST_EUR_EAS', 'PBS_UZB', 'PBS_EUR', 'PBS_EAS']],
        on='SNP', how='inner'
    )
    
    print(f"  Combined dataset: {len(combined):,} SNPs")
    
    # ---- Identify Uzbek-specific SNPs ----
    print("\n[7] Identifying Uzbek-specific SNPs...")
    
    # Criterion 1: High PBS (Uzbek branch lengthening)
    pbs_sig = combined['PBS_UZB'] >= args.pbs_threshold
    print(f"  PBS >= {args.pbs_threshold}: {pbs_sig.sum():,} SNPs")
    
    # Criterion 2: High delta-AF vs ALL populations (min delta-AF)
    daf_sig = combined['dAF_min'] >= args.delta_af_threshold
    print(f"  min(dAF) >= {args.delta_af_threshold}: {daf_sig.sum():,} SNPs")
    
    # Criterion 3: UZB-specific by Fst (high vs all, not just one pair)
    fst_all_high = (
        (combined['FST_UZB_EUR'] >= 0.1) & 
        (combined['FST_UZB_EAS'] >= 0.1) & 
        (combined['FST_UZB_SAS'].fillna(0) >= 0.05)
    )
    print(f"  Fst >= 0.1 vs EUR & EAS, >= 0.05 vs SAS: {fst_all_high.sum():,} SNPs")
    
    # Combined: PBS high OR (dAF high AND Fst-all-high)
    uzbek_specific = combined[pbs_sig | (daf_sig & fst_all_high)].copy()
    
    # Sort by PBS
    uzbek_specific = uzbek_specific.sort_values('PBS_UZB', ascending=False)
    
    print(f"\n  ★ Uzbek-specific candidates: {len(uzbek_specific):,} SNPs")
    print(f"    (PBS >= {args.pbs_threshold} OR [min_dAF >= {args.delta_af_threshold} AND multi-pop Fst high])")
    
    # ---- Near-private variants (UZB common, rare elsewhere) ----
    print("\n[8] Identifying near-private variants...")
    
    uzb_common = combined['MAF_UZB'] >= 0.05
    eur_rare = combined['MAF_EUR'] <= 0.01
    eas_rare = combined['MAF_EAS'] <= 0.01
    sas_rare = combined['MAF_SAS'] <= 0.01
    afr_rare = combined['MAF_AFR'] <= 0.01
    
    # Private vs ALL (very strict)
    private_all = combined[uzb_common & eur_rare & eas_rare & sas_rare & afr_rare].copy()
    print(f"  UZB MAF>=5%, ALL others MAF<=1%: {len(private_all):,} SNPs")
    
    # Private vs non-African (more relevant for Central Asia)
    private_nonafr = combined[uzb_common & eur_rare & eas_rare & sas_rare].copy()
    print(f"  UZB MAF>=5%, EUR+EAS+SAS MAF<=1%: {len(private_nonafr):,} SNPs")
    
    # Relaxed: UZB elevated, others low
    uzb_elevated = combined['MAF_UZB'] >= 0.10
    others_low = (
        (combined['MAF_EUR'] <= 0.03) & 
        (combined['MAF_EAS'] <= 0.03) & 
        (combined['MAF_SAS'] <= 0.03)
    )
    enriched_uzb = combined[uzb_elevated & others_low].copy()
    print(f"  UZB MAF>=10%, EUR+EAS+SAS MAF<=3%: {len(enriched_uzb):,} SNPs")
    
    # ---- Save outputs ----
    print("\n[9] Saving results...")
    
    # Full PBS results
    pbs_out = combined[['CHR', 'SNP', 'MAF_UZB', 'MAF_EUR', 'MAF_EAS', 'MAF_SAS', 'MAF_AFR',
                        'dAF_EUR', 'dAF_EAS', 'dAF_SAS', 'dAF_AFR', 'dAF_min', 'dAF_mean',
                        'FST_UZB_EUR', 'FST_UZB_EAS', 'FST_UZB_SAS', 'FST_UZB_AFR',
                        'PBS_UZB', 'PBS_EUR', 'PBS_EAS']].copy()
    pbs_out = pbs_out.sort_values('PBS_UZB', ascending=False)
    pbs_out.to_csv(outdir / 'pbs_results.tsv', sep='\t', index=False, float_format='%.6f')
    print(f"  ✓ pbs_results.tsv ({len(pbs_out):,} rows)")
    
    # Uzbek-specific candidates
    uzbek_specific.to_csv(outdir / 'uzbek_specific_snps.tsv', sep='\t', index=False, float_format='%.6f')
    print(f"  ✓ uzbek_specific_snps.tsv ({len(uzbek_specific):,} rows)")
    
    # Delta-AF summary
    daf_out = combined[['CHR', 'SNP', 'MAF_UZB', 'MAF_EUR', 'MAF_EAS', 'MAF_SAS', 'MAF_AFR',
                        'dAF_EUR', 'dAF_EAS', 'dAF_SAS', 'dAF_AFR', 'dAF_min', 'dAF_mean']].copy()
    daf_out = daf_out.sort_values('dAF_min', ascending=False)
    daf_out.to_csv(outdir / 'delta_af_all.tsv', sep='\t', index=False, float_format='%.6f')
    print(f"  ✓ delta_af_all.tsv ({len(daf_out):,} rows)")
    
    # Near-private variants
    if len(private_nonafr) > 0:
        private_nonafr_out = private_nonafr[['CHR', 'SNP', 'MAF_UZB', 'MAF_EUR', 'MAF_EAS', 
                                              'MAF_SAS', 'MAF_AFR']].sort_values('MAF_UZB', ascending=False)
        private_nonafr_out.to_csv(outdir / 'near_private_variants.tsv', sep='\t', index=False, float_format='%.6f')
        print(f"  ✓ near_private_variants.tsv ({len(private_nonafr_out):,} rows)")
    
    # ---- Summary statistics ----
    print("\n[10] Writing summary...")
    
    summary_lines = [
        "=" * 70,
        "UZBEK-SPECIFIC SNP ANALYSIS SUMMARY",
        "=" * 70,
        "",
        "DATASET:",
        f"  Total SNPs analyzed: {len(combined):,}",
        f"  Populations: UZB (n=1199), EUR (n=503), EAS (n=504), SAS (n=489), AFR (n=661)",
        "",
        "PAIRWISE FST (weighted estimates from PLINK):",
        f"  UZB vs EUR: see existing analysis",
        f"  UZB vs EAS: see fst_UZB_vs_EAS.log",
        f"  UZB vs SAS: see fst_UZB_vs_SAS.log",
        f"  UZB vs AFR: see fst_UZB_vs_AFR.log",
        f"  EUR vs EAS: see fst_EUR_vs_EAS.log",
        "",
        "POPULATION BRANCH STATISTIC (PBS):",
        f"  Triangle: UZB – EUR – EAS",
        f"  PBS formula: PBS_UZB = (T_UZB_EUR + T_UZB_EAS - T_EUR_EAS) / 2",
        f"  where T = -log(1 - Fst)",
        f"  Mean PBS_UZB:   {fst_merge['PBS_UZB'].mean():.6f}",
        f"  Median PBS_UZB: {fst_merge['PBS_UZB'].median():.6f}",
        f"  Max PBS_UZB:    {fst_merge['PBS_UZB'].max():.4f}",
        f"  95th percentile:  {fst_merge['PBS_UZB'].quantile(0.95):.4f}",
        f"  99th percentile:  {fst_merge['PBS_UZB'].quantile(0.99):.4f}",
        f"  99.9th percentile: {fst_merge['PBS_UZB'].quantile(0.999):.4f}",
        "",
        "ALLELE FREQUENCY DIFFERENCES (delta-AF):",
        f"  Mean |dAF| UZB-EUR: {merged['dAF_EUR'].mean():.4f}",
        f"  Mean |dAF| UZB-EAS: {merged['dAF_EAS'].mean():.4f}",
        f"  Mean |dAF| UZB-SAS: {merged['dAF_SAS'].mean():.4f}",
        f"  Mean |dAF| UZB-AFR: {merged['dAF_AFR'].mean():.4f}",
        "",
        "UZBEK-SPECIFIC CANDIDATES:",
        f"  PBS >= {args.pbs_threshold}: {pbs_sig.sum():,} SNPs",
        f"  min(dAF) >= {args.delta_af_threshold} & multi-pop Fst: {(daf_sig & fst_all_high).sum():,} SNPs",
        f"  Combined (union): {len(uzbek_specific):,} SNPs",
        "",
        "NEAR-PRIVATE VARIANTS (UZB MAF≥5%):",
        f"  Rare in ALL populations (MAF≤1%): {len(private_all):,} SNPs",
        f"  Rare in EUR+EAS+SAS (MAF≤1%): {len(private_nonafr):,} SNPs",
        f"  Enriched in UZB (MAF≥10%, others MAF≤3%): {len(enriched_uzb):,} SNPs",
        "",
        "TOP 20 PBS HITS:",
        f"  {'SNP':<20} {'CHR':>4} {'PBS_UZB':>10} {'MAF_UZB':>9} {'MAF_EUR':>9} {'MAF_EAS':>9} {'MAF_SAS':>9}",
        "-" * 80,
    ]
    
    top20 = pbs_out.head(20)
    for _, row in top20.iterrows():
        summary_lines.append(
            f"  {row['SNP']:<20} {int(row['CHR']):>4} {row['PBS_UZB']:>10.4f} "
            f"{row['MAF_UZB']:>9.4f} {row['MAF_EUR']:>9.4f} {row['MAF_EAS']:>9.4f} "
            f"{row.get('MAF_SAS', float('nan')):>9.4f}"
        )
    
    summary_lines.extend([
        "",
        "=" * 70,
        "Output files:",
        f"  pbs_results.tsv          - Full PBS scores ({len(pbs_out):,} SNPs)",
        f"  uzbek_specific_snps.tsv  - Uzbek-specific candidates ({len(uzbek_specific):,} SNPs)",
        f"  delta_af_all.tsv         - Allele frequency differences ({len(daf_out):,} SNPs)",
        f"  near_private_variants.tsv - Near-private variants",
        f"  summary_stats.txt        - This file",
        "=" * 70,
    ])
    
    summary = '\n'.join(summary_lines)
    (outdir / 'summary_stats.txt').write_text(summary)
    print(f"  ✓ summary_stats.txt")
    
    print("\n" + summary)
    print("\nDone!")


if __name__ == '__main__':
    main()
