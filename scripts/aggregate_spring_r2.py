#!/usr/bin/env python3
"""Aggregate per-chromosome R² stats from spring_stats.txt into step4.html values.

Usage:
  1. Copy /tmp/spring_stats.txt from DRAGEN to data/spring_stats.txt
  2. python scripts/aggregate_spring_r2.py
  3. Copy output values into step4.html
"""
import sys, re, os

STATS_FILE = os.path.join(os.path.dirname(__file__), '..', 'data', 'spring_stats.txt')

# Parse CHR lines
totals = {'N': 0, 'SUM': 0.0}
bin_totals = [0] * 10
ge_totals = {'GE30': 0, 'GE50': 0, 'GE80': 0, 'GE90': 0}
chr_data = []

with open(STATS_FILE, encoding='utf-8') as f:
    for line in f:
        if not line.startswith('CHR '):
            continue
        parts = line.strip().split()
        d = {}
        i = 0
        while i < len(parts) - 1:
            d[parts[i]] = parts[i+1]
            i += 2
        n = int(d['N'])
        s = float(d['SUM'])
        bins = [int(d[f'B{j}']) for j in range(10)]
        ge30 = int(d['GE30'])
        ge50 = int(d['GE50'])
        ge80 = int(d['GE80'])
        ge90 = int(d['GE90'])
        
        chr_num = int(d['CHR'])
        chr_data.append((chr_num, n, s, bins, ge30, ge50, ge80, ge90))
        totals['N'] += n
        totals['SUM'] += s
        for j in range(10):
            bin_totals[j] += bins[j]
        ge_totals['GE30'] += ge30
        ge_totals['GE50'] += ge50
        ge_totals['GE80'] += ge80
        ge_totals['GE90'] += ge90

if not chr_data:
    print("ERROR: No CHR lines found in", STATS_FILE)
    sys.exit(1)

N = totals['N']
mean_r2 = totals['SUM'] / N
lt30 = N - ge_totals['GE30']
b30_80 = ge_totals['GE30'] - ge_totals['GE80']
b80_90 = ge_totals['GE80'] - ge_totals['GE90']

print(f"=== Spring 2026 Aggregated R² Statistics ({len(chr_data)} chromosomes) ===\n")
print(f"Total variants:  {N:,}")
print(f"Sum R²:          {totals['SUM']:,.2f}")
print(f"Mean R²:         {mean_r2:.6f}")
print()

print("--- Histogram bins (for step4.html JS) ---")
labels = ['0.0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5',
          '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1.0']
for i, (lbl, cnt) in enumerate(zip(labels, bin_totals)):
    print(f"  B{i} ({lbl}): {cnt:>12,}  ({cnt/N*100:5.1f}%)")
print()

print("--- Threshold counts ---")
print(f"  < 0.30:  {lt30:>12,}  ({lt30/N*100:.1f}%)")
print(f"  0.30-0.80: {b30_80:>10,}  ({b30_80/N*100:.1f}%)")
print(f"  0.80-0.90: {b80_90:>10,}  ({b80_90/N*100:.1f}%)")
print(f"  >= 0.90: {ge_totals['GE90']:>12,}  ({ge_totals['GE90']/N*100:.1f}%)")
print(f"  >= 0.80: {ge_totals['GE80']:>12,}  ({ge_totals['GE80']/N*100:.1f}%)")
print(f"  >= 0.30: {ge_totals['GE30']:>12,}  ({ge_totals['GE30']/N*100:.1f}%)")
print()

print("--- step4.html JS bin values (copy-paste) ---")
for i, cnt in enumerate(bin_totals):
    lo = i * 0.1
    hi = (i + 1) * 0.1 if i < 9 else 1.01
    print(f"            {{lo: {lo:.2f}, hi: {hi:.2f}, count: {cnt}, label: '{lo:.1f}'}},")
print()

density = N / 456684
print(f"--- Stat card values ---")
print(f"  Total output: {N:,}  ({N/1e6:.2f}M)")
print(f"  Mean INFO:    {mean_r2:.4f}")
print(f"  Density increase: {density:.0f}×")
print()

print(f"--- Chart subtitle ---")
print(f"  N = {N:,} imputed variants | Mean R² = {mean_r2:.4f}")
print()

print("--- Per-chromosome breakdown ---")
for chr_num, n, s, bins, ge30, ge50, ge80, ge90 in sorted(chr_data):
    mr2 = s / n if n > 0 else 0
    print(f"  chr{chr_num:2d}: {n:>10,} variants, mean R² = {mr2:.4f}, "
          f"≥0.80: {ge80:>8,} ({ge80/n*100:.1f}%)")
print(f"  TOTAL: {N:>10,} variants, mean R² = {mean_r2:.4f}, "
      f"≥0.80: {ge_totals['GE80']:>8,} ({ge_totals['GE80']/N*100:.1f}%)")
