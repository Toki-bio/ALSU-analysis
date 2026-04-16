"""Generate JSON data for pre-GWAS pregnancy outcome visualization."""
import csv, json, sys
from collections import Counter, defaultdict

with open('GWAS от 27.08 - Sheet1.csv', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

SLOT_START = 28; SLOT_SIZE = 7; MAX_SLOTS = 20
LOSS_CODES = {'2','3','6'}
CODE_NAMES = {
    '1':'Live birth','2':'Stillbirth','3':'Miscarriage',
    '4':'Elective abortion','5':'Ectopic','6':'Non-developing',
    '7':'Medical abort','8':'Currently pregnant'
}

samples = []
all_events = []

for row in rows:
    sid = row[0].strip()
    if not sid: continue
    age_raw = row[7].strip()
    try: age = int(age_raw)
    except: age = None
    orig_label = row[1].strip()  # 0=case, 1=control in this spreadsheet

    events = []
    for slot in range(MAX_SLOTS):
        base = SLOT_START + slot * SLOT_SIZE
        if base + 4 >= len(row): break
        outcome = row[base + 3].strip()
        ga_raw = row[base + 4].strip()
        if not outcome or outcome == 'NA': continue
        try: ga = float(ga_raw) if ga_raw else None
        except: ga = None
        events.append({'outcome': outcome, 'ga': ga, 'slot': slot+1})
        all_events.append({'sid': sid, 'outcome': outcome, 'ga': ga})

    loss_strict = sum(1 for e in events if e['outcome'] in LOSS_CODES and e['ga'] is not None and e['ga'] < 20)
    loss_no_ga = sum(1 for e in events if e['outcome'] in LOSS_CODES and e['ga'] is None)
    loss_any = sum(1 for e in events if e['outcome'] in LOSS_CODES)
    loss_ge20 = sum(1 for e in events if e['outcome'] in LOSS_CODES and e['ga'] is not None and e['ga'] >= 20)
    live = sum(1 for e in events if e['outcome'] == '1')
    abortions = sum(1 for e in events if e['outcome'] == '4')
    ectopic = sum(1 for e in events if e['outcome'] == '5')
    medical = sum(1 for e in events if e['outcome'] == '7')
    preg_now = sum(1 for e in events if e['outcome'] == '8')
    n_preg = len(events)

    # Classification
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
    elif n_preg == 0:
        cls = 'UNCLASSIFIED'
    else:
        cls = 'UNCLASSIFIED'

    # Bin for gradient visualization
    if loss_strict >= 5: bin_label = '5+ losses'
    elif loss_strict == 4: bin_label = '4 losses'
    elif loss_strict == 3: bin_label = '3 losses'
    elif loss_strict == 2 and live == 0: bin_label = '2 losses, 0 live'
    elif loss_strict == 2 and live >= 1: bin_label = '2 losses + live'
    elif loss_any == 1 and loss_ge20 > 0: bin_label = '1 loss (GA≥20)'
    elif loss_any == 1 and live == 0: bin_label = '1 loss, 0 live'
    elif loss_any == 1 and live >= 1: bin_label = '1 loss + live'
    elif loss_any == 0 and live == 0 and (ectopic > 0 or medical > 0): bin_label = 'Ectopic/medical only'
    elif loss_any == 0 and live == 0: bin_label = 'No pregnancies'
    elif loss_any == 0 and live == 1: bin_label = '1 live birth only'
    elif loss_any == 0 and live == 2: bin_label = '2 live births'
    elif loss_any == 0 and live >= 3: bin_label = '3+ live births'
    else: bin_label = 'Other'

    samples.append({
        'sid': sid, 'age': age, 'orig_label': orig_label,
        'n_preg': n_preg, 'loss_strict': loss_strict, 'loss_any': loss_any,
        'loss_ge20': loss_ge20, 'loss_no_ga': loss_no_ga,
        'live': live, 'ectopic': ectopic, 'medical': medical,
        'abortions': abortions, 'preg_now': preg_now,
        'cls': cls, 'bin': bin_label,
        'code_seq': ''.join(e['outcome'] for e in events)
    })

# === Build output data ===
data = {}

# 1. Summary stats
cls_counts = Counter(s['cls'] for s in samples)
data['summary'] = {
    'total': len(samples),
    'age_18_45': sum(1 for s in samples if s['age'] and 18 <= s['age'] <= 45),
    'classification': {k: cls_counts.get(k, 0) for k in [
        'CASE','CASE_PROBABLE','CONTROL','CONTROL_PROBABLE',
        'GRAY','REVIEW','EXCLUDE_AGE','EXCLUDE_ECTOPIC','EXCLUDE_MEDICAL','UNCLASSIFIED'
    ]}
}

# 2. Outcome code distribution
outcome_dist = Counter(e['outcome'] for e in all_events if e['outcome'] not in ('NA',''))
data['outcome_codes'] = [
    {'code': c, 'name': CODE_NAMES.get(c, c), 'count': outcome_dist.get(c, 0)}
    for c in ['1','2','3','4','5','6','7','8']
]

# 3. GA distribution per loss code
ga_data = defaultdict(list)
for e in all_events:
    if e['outcome'] in LOSS_CODES and e['ga'] is not None:
        ga_data[e['outcome']].append(e['ga'])
data['ga_distributions'] = {}
for code in ['2','3','6']:
    vals = sorted(ga_data.get(code, []))
    # Bin into weeks
    bins = defaultdict(int)
    for v in vals:
        week = int(v)
        bins[week] += 1
    data['ga_distributions'][code] = {
        'name': CODE_NAMES[code],
        'bins': [{'week': w, 'count': bins[w]} for w in sorted(bins.keys())],
        'n': len(vals),
        'n_ge20': sum(1 for v in vals if v >= 20),
        'n_lt20': sum(1 for v in vals if v < 20)
    }

# 4. Gradient bins (case → control spectrum)
BIN_ORDER = [
    '5+ losses', '4 losses', '3 losses', '2 losses, 0 live', '2 losses + live',
    '1 loss (GA≥20)', '1 loss, 0 live', '1 loss + live',
    'Ectopic/medical only', 'No pregnancies',
    '1 live birth only', '2 live births', '3+ live births', 'Other'
]
bin_counts = Counter(s['bin'] for s in samples)
data['gradient_bins'] = []
for b in BIN_ORDER:
    if bin_counts.get(b, 0) == 0: continue
    ss = [s for s in samples if s['bin'] == b]
    age_vals = [s['age'] for s in ss if s['age'] is not None]
    data['gradient_bins'].append({
        'label': b,
        'n': bin_counts[b],
        'mean_age': round(sum(age_vals)/len(age_vals), 1) if age_vals else None,
        'cls_dist': dict(Counter(s['cls'] for s in ss)),
        'orig_label_dist': dict(Counter(s['orig_label'] for s in ss)),
        'mean_loss': round(sum(s['loss_any'] for s in ss) / len(ss), 2),
        'mean_live': round(sum(s['live'] for s in ss) / len(ss), 2)
    })

# 5. Loss count distribution
loss_dist = Counter(s['loss_strict'] for s in samples)
data['loss_distribution'] = [{'losses': k, 'count': loss_dist[k]} for k in sorted(loss_dist.keys())]

# 6. Label cross-tab
xtab = Counter()
for s in samples:
    orig = 'case(0)' if s['orig_label'] == '0' else 'ctrl(1)' if s['orig_label'] == '1' else 'other'
    comp = 'CASE' if s['cls'] in ('CASE','CASE_PROBABLE') else 'CTRL' if s['cls'] in ('CONTROL','CONTROL_PROBABLE') else 'OTHER'
    xtab[(orig, comp)] += 1
data['label_crosstab'] = [
    {'orig': o, 'computed': c, 'n': xtab[(o,c)]}
    for o, c in sorted(xtab.keys())
]

# 7. Top sequence patterns
pat_counts = Counter(s['code_seq'] for s in samples if s['code_seq'])
data['top_patterns'] = []
for pat, cnt in pat_counts.most_common(25):
    loss_in = sum(1 for c in pat if c in '236')
    live_in = sum(1 for c in pat if c == '1')
    abort_in = sum(1 for c in pat if c == '4')
    readable = ' → '.join(CODE_NAMES.get(c, c) for c in pat)
    data['top_patterns'].append({
        'pattern': pat, 'readable': readable, 'count': cnt,
        'n_loss': loss_in, 'n_live': live_in, 'n_abort': abort_in
    })

# 8. Age distribution by classification
for cls_name in ['CASE','CONTROL','GRAY','REVIEW']:
    ages = [s['age'] for s in samples if s['cls'] == cls_name and s['age'] is not None]
    age_bins = defaultdict(int)
    for a in ages:
        bucket = (a // 5) * 5  # 5-year bins
        age_bins[bucket] += 1
    data[f'age_dist_{cls_name.lower()}'] = [
        {'age_bin': f'{b}-{b+4}', 'count': age_bins[b]}
        for b in sorted(age_bins.keys())
    ]

with open('data/pregnancy_analysis.json', 'w', encoding='utf-8') as f:
    json.dump(data, f, indent=2, ensure_ascii=False)

print(f'Written data/pregnancy_analysis.json ({len(samples)} samples)')
print(f'Gradient bins: {len(data["gradient_bins"])}')
print(f'Top patterns: {len(data["top_patterns"])}')
for b in data['gradient_bins']:
    print(f'  {b["label"]:25s}: {b["n"]:4d}')
