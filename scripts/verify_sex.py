import json, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
d = json.load(open('data/investigation_data_v2.json'))
es = d['enriched_samples']
vals = [s['sex_f'] for s in es if s.get('sex_f') is not None]
print(f'Total with sex_f: {len(vals)}')
print(f'Min: {min(vals):.4f}')
print(f'Max: {max(vals):.4f}')
print(f'Mean: {sum(vals)/len(vals):.4f}')
print(f'F < 0.2 (FEMALE): {sum(1 for v in vals if v < 0.2)}')
print(f'F > 0.8 (MALE):   {sum(1 for v in vals if v > 0.8)}')
print(f'0.2-0.8 (AMBIG):  {sum(1 for v in vals if 0.2 <= v <= 0.8)}')
print()
print("=== AMBIGUOUS SAMPLES (F >= 0.2) ===")
print(f"{'IID':14s} {'F(X)':>8s} {'F_MISS':>8s} {'HET':>8s} {'CHIP':>10s}")
ambig = [s for s in es if s.get('sex_f') is not None and s['sex_f'] >= 0.2]
for s in sorted(ambig, key=lambda x: x['sex_f'], reverse=True):
    print(f"{s['iid']:14s} {s['sex_f']:8.4f} {s.get('fmiss',0):8.4f} {s.get('het_rate',0):8.4f} {str(s.get('barcode',''))[-8:]:>10s}")

print()
print("=== F(X) DISTRIBUTION (histogram bins) ===")
import collections
bins = collections.Counter()
for v in vals:
    b = round(v * 20) / 20  # 0.05 bins
    bins[b] += 1
for b in sorted(bins):
    print(f"  {b:.2f}: {'#' * bins[b]} ({bins[b]})")
