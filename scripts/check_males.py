import json, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
d = json.load(open('data/investigation_data_v2.json'))
es = d['enriched_samples']
verdicts = {s['iid']: s for s in d['sample_verdicts']}

males = [s for s in es if s.get('sex_f') is not None and s['sex_f'] < 0.2]
females = [s for s in es if s.get('sex_f') is not None and s['sex_f'] > 0.8]
ambig = [s for s in es if s.get('sex_f') is not None and 0.2 <= s['sex_f'] <= 0.8]

print(f"Female (F>0.8): {len(females)}")
print(f"Male (F<0.2):   {len(males)}")
print(f"Ambiguous:       {len(ambig)}")

keep_males = [s for s in males if verdicts.get(s['iid'], {}).get('action') == 'KEEP']
print(f"\nMales KEEP: {len(keep_males)}")
print(f"Males REMOVE: {len(males) - len(keep_males)}")

print(f"\nALL MALES:")
print("IID            F(X)     fmiss    het      chip   action   id_status")
for s in sorted(males, key=lambda x: x['sex_f']):
    v = verdicts.get(s['iid'], {})
    print(f"{s['iid']:14s} {s['sex_f']:.4f}  {s.get('fmiss',0):.4f}  {s.get('het_rate',0):.4f}  {str(s.get('barcode',''))[-4:]:6s} {v.get('action','?'):8s} {v.get('identity_status','?')}")
