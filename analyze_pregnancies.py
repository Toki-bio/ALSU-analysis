"""Comprehensive pregnancy outcome analysis for case/control classification."""
import csv, statistics
from collections import Counter, defaultdict

with open('GWAS от 27.08 - Sheet1.csv', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

SLOT_START = 28; SLOT_SIZE = 7; MAX_SLOTS = 20
LOSS_CODES = {'2','3','6'}  # intrauterine spontaneous losses

samples = []
for row in rows:
    sid = row[0].strip()
    if not sid: continue
    age = row[7].strip()
    try: age_int = int(age)
    except: age_int = None
    case_label = row[1].strip()
    currently_preg = row[26].strip()

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

    # Clean loss count: codes 2,3,6 with GA<20
    loss_strict = sum(1 for e in events if e['outcome'] in LOSS_CODES and e['ga'] is not None and e['ga'] < 20)
    loss_no_ga = sum(1 for e in events if e['outcome'] in LOSS_CODES and e['ga'] is None)
    loss_any = sum(1 for e in events if e['outcome'] in LOSS_CODES)
    loss_ge20 = sum(1 for e in events if e['outcome'] in LOSS_CODES and e['ga'] is not None and e['ga'] >= 20)

    live = sum(1 for e in events if e['outcome'] == '1')
    abortions = sum(1 for e in events if e['outcome'] == '4')
    ectopic = sum(1 for e in events if e['outcome'] == '5')
    medical = sum(1 for e in events if e['outcome'] == '7')
    preg_now = sum(1 for e in events if e['outcome'] == '8')

    code_seq = ''.join(e['outcome'] for e in events)

    samples.append({
        'sid': sid, 'age': age_int, 'label': case_label,
        'n_events': len(events), 'code_seq': code_seq,
        'loss_strict': loss_strict, 'loss_no_ga': loss_no_ga,
        'loss_any': loss_any, 'loss_ge20': loss_ge20,
        'live': live, 'ectopic': ectopic, 'medical': medical,
        'abortions': abortions, 'preg_now': preg_now
    })

print(f'Total samples parsed: {len(samples)}')
print(f'Age available: {sum(1 for s in samples if s["age"] is not None)}')
print(f'Age 18-45: {sum(1 for s in samples if s["age"] and 18 <= s["age"] <= 45)}')
print()

# === LOSS COUNT DISTRIBUTION ===
print('=== LOSS COUNT (strict: codes 2,3,6, GA<20) ===')
loss_dist = Counter(s['loss_strict'] for s in samples)
for k in sorted(loss_dist): print(f'  {k} losses: {loss_dist[k]} samples')

print()
print('=== LOSS COUNT (any: codes 2,3,6, ALL GA) ===')
loss_dist2 = Counter(s['loss_any'] for s in samples)
for k in sorted(loss_dist2): print(f'  {k} losses: {loss_dist2[k]} samples')

print()
print('=== LOSSES WITH GA>=20 (need review: stillbirth vs late miscarriage) ===')
ge20_events = []
for row in rows:
    sid = row[0].strip()
    if not sid: continue
    for slot in range(MAX_SLOTS):
        base = SLOT_START + slot * SLOT_SIZE
        if base + 4 >= len(row): break
        outcome = row[base + 3].strip()
        ga_raw = row[base + 4].strip()
        if outcome in LOSS_CODES:
            try:
                ga = float(ga_raw) if ga_raw else None
            except:
                ga = None
            if ga is not None and ga >= 20:
                ge20_events.append({'sid': sid, 'outcome': outcome, 'ga': ga, 'slot': slot+1})

print(f'  Total events with code 2/3/6 and GA>=20: {len(ge20_events)}')
code_dist = Counter(e['outcome'] for e in ge20_events)
for c in sorted(code_dist):
    names = {'2':'stillbirth','3':'miscarriage','6':'non-developing'}
    gas = [e['ga'] for e in ge20_events if e['outcome'] == c]
    print(f'  Code {c} ({names.get(c,"?")}): n={code_dist[c]}, GA range={min(gas):.0f}-{max(gas):.0f}, median={statistics.median(gas):.0f}')

print()
print('=== SAMPLES WITH MISSING GA ON LOSS EVENTS ===')
noga = Counter(s['loss_no_ga'] for s in samples)
for k in sorted(noga):
    if k > 0: print(f'  {k} loss events without GA: {noga[k]} samples')
total_noga = sum(s['loss_no_ga'] for s in samples)
print(f'  Total loss events without GA: {total_noga}')

print()
print('=== ORIGINAL CASE/CONTROL LABELS ===')
label_dist = Counter(s['label'] for s in samples)
for k in sorted(label_dist, key=lambda x: -label_dist[x]):
    print(f'  "{k}": {label_dist[k]}')

# === CASE/CONTROL CLASSIFICATION ===
print()
print('='*70)
print('=== PROPOSED CASE/CONTROL CLASSIFICATION ===')
print('='*70)

cases = []; controls = []; exclude = []; gray = []
for s in samples:
    age_ok = s['age'] is not None and 18 <= s['age'] <= 45
    # For strict: only count GA<20 losses
    strict_loss = s['loss_strict']
    # For inclusive: also count losses where GA is missing (assume <20 for codes 3,6)
    inclusive_loss = s['loss_strict'] + s['loss_no_ga']
    has_code8 = s['preg_now'] > 0
    only_ectopic = s['ectopic'] > 0 and s['loss_any'] == 0
    only_medical = s['medical'] > 0 and s['loss_any'] == 0

    classification = None
    reason = ''

    if only_ectopic:
        classification = 'EXCLUDE'
        reason = 'ectopic-only losses'
    elif only_medical and s['loss_any'] == 0:
        classification = 'EXCLUDE'
        reason = 'medical-abort-only'
    elif s['loss_ge20'] > 0 and strict_loss < 2:
        classification = 'REVIEW'
        reason = f'{s["loss_ge20"]} losses GA>=20wk (late stillbirth?)'
    elif strict_loss >= 2 and age_ok:
        classification = 'CASE'
        reason = f'{strict_loss} strict losses'
    elif inclusive_loss >= 2 and strict_loss < 2 and age_ok:
        classification = 'CASE_IF_GA'
        reason = f'{strict_loss} strict + {s["loss_no_ga"]} missing-GA losses'
    elif s['loss_any'] == 0 and s['live'] >= 1 and not has_code8 and age_ok:
        classification = 'CONTROL'
        reason = f'{s["live"]} live births, 0 losses'
    elif s['loss_any'] == 0 and s['live'] >= 1 and has_code8 and age_ok:
        classification = 'CONTROL_IF_DELIVERED'
        reason = f'{s["live"]} live + currently pregnant'
    elif s['loss_any'] == 1 and age_ok:
        classification = 'GRAY_1LOSS'
        reason = '1 loss (neither case nor control)'
    elif not age_ok:
        classification = 'EXCLUDE_AGE'
        reason = f'age={s["age"]}'
    else:
        classification = 'UNCLASSIFIED'
        reason = f'seq={s["code_seq"]}, live={s["live"]}, loss={s["loss_any"]}'

    s['classification'] = classification
    s['reason'] = reason

class_dist = Counter(s['classification'] for s in samples)
for cls in ['CASE','CASE_IF_GA','CONTROL','CONTROL_IF_DELIVERED','GRAY_1LOSS','REVIEW','EXCLUDE','EXCLUDE_AGE','UNCLASSIFIED']:
    if cls in class_dist:
        print(f'  {cls:25s}: {class_dist[cls]:5d}')

print()
print('=== CASE BREAKDOWN (by strict loss count) ===')
case_loss = Counter(s['loss_strict'] for s in samples if s['classification'] == 'CASE')
for k in sorted(case_loss): print(f'  {k} losses: {case_loss[k]} cases')

print()
print('=== REVIEW SAMPLES (losses at GA>=20) ===')
for s in samples:
    if s['classification'] == 'REVIEW':
        print(f'  {s["sid"]}: seq={s["code_seq"]}, strict={s["loss_strict"]}, ge20={s["loss_ge20"]}, live={s["live"]}')
    if s['classification'] == 'REVIEW' and sum(1 for _ in [0]) > 20:
        print('  ... (truncated)')
        break

# === PATTERN-BASED GROUPING ===
print()
print('='*70)
print('=== TOP 30 PREGNANCY SEQUENCE PATTERNS ===')
print('='*70)
print('  Code legend: 1=live, 2=stillbirth, 3=miscarriage, 4=abortion,')
print('  5=ectopic, 6=non-developing, 7=med-abort, 8=currently-pregnant')
print()
pattern_counts = Counter()
pattern_examples = defaultdict(list)
for s in samples:
    p = s['code_seq']
    pattern_counts[p] += 1
    if len(pattern_examples[p]) < 3:
        pattern_examples[p].append(s['sid'])

for pat, cnt in pattern_counts.most_common(30):
    # Classify pattern
    loss_in_pat = sum(1 for c in pat if c in '236')
    live_in_pat = sum(1 for c in pat if c == '1')
    print(f'  {pat:20s}  n={cnt:4d}  (live={live_in_pat}, loss={loss_in_pat})')

# === GRADIENT VIEW ===
print()
print('='*70)
print('=== CASE→CONTROL GRADIENT (aggregated pattern groups) ===')
print('='*70)

# Group by (loss_count, live_count, has_ectopic, has_medical, has_code8)
gradient = defaultdict(lambda: {'n':0, 'sids':[]})
for s in samples:
    key = (s['loss_any'], s['loss_ge20'], s['live'], s['ectopic']>0, s['medical']>0, s['preg_now']>0, s['abortions'])
    gradient[key]['n'] += 1
    if len(gradient[key]['sids']) < 2:
        gradient[key]['sids'].append(s['sid'])

print(f'  {"Loss":>4} {"≥20w":>4} {"Live":>4} {"Ect":>3} {"Med":>3} {"Prg":>3} {"Abt":>3} | {"N":>5}  Classification')
print(f'  {"-"*4} {"-"*4} {"-"*4} {"-"*3} {"-"*3} {"-"*3} {"-"*3} | {"-"*5}  {"-"*30}')
for key in sorted(gradient.keys(), key=lambda k: (-k[0], -k[1], k[2])):
    loss, ge20, live, ect, med, prg, abt = key
    n = gradient[key]['n']
    # Classify
    if loss >= 2 and ge20 == 0: cls = 'CASE'
    elif loss >= 2 and ge20 > 0: cls = 'CASE+REVIEW(late)'
    elif loss == 1: cls = 'GRAY (1 loss)'
    elif loss == 0 and live >= 1 and not prg: cls = 'CONTROL'
    elif loss == 0 and live >= 1 and prg: cls = 'CONTROL (if delivered)'
    elif loss == 0 and live == 0 and ect: cls = 'EXCLUDE (ectopic only)'
    elif loss == 0 and live == 0 and med: cls = 'EXCLUDE (med abort only)'
    else: cls = 'UNCLASSIFIED'
    ect_s = 'Y' if ect else '.'
    med_s = 'Y' if med else '.'
    prg_s = 'Y' if prg else '.'
    print(f'  {loss:4d} {ge20:4d} {live:4d} {ect_s:>3} {med_s:>3} {prg_s:>3} {abt:3d} | {n:5d}  {cls}')

# === SUMMARY STATS ===
print()
print('='*70)
print('=== FINAL SUMMARY ===')
print('='*70)
n_case = sum(1 for s in samples if s['classification'] == 'CASE')
n_case_if = sum(1 for s in samples if s['classification'] == 'CASE_IF_GA')
n_ctrl = sum(1 for s in samples if s['classification'] == 'CONTROL')
n_ctrl_if = sum(1 for s in samples if s['classification'] == 'CONTROL_IF_DELIVERED')
n_gray = sum(1 for s in samples if s['classification'] == 'GRAY_1LOSS')
n_review = sum(1 for s in samples if s['classification'] == 'REVIEW')
n_excl = sum(1 for s in samples if 'EXCLUDE' in s['classification'])
n_unc = sum(1 for s in samples if s['classification'] == 'UNCLASSIFIED')

print(f'  Definite CASES (>=2 losses GA<20, age 18-45):  {n_case}')
print(f'  Probable CASES (if missing GA confirmed <20):  {n_case_if}')
print(f'  Definite CONTROLS (0 loss, >=1 live, age 18-45): {n_ctrl}')
print(f'  Probable CONTROLS (if currently preg delivered): {n_ctrl_if}')
print(f'  GRAY ZONE (exactly 1 loss):                     {n_gray}')
print(f'  REVIEW NEEDED (losses at GA>=20):                {n_review}')
print(f'  EXCLUDED (age/ectopic-only/med-abort-only):      {n_excl}')
print(f'  UNCLASSIFIED:                                    {n_unc}')
print(f'  TOTAL:                                           {len(samples)}')
# Cross-tab original labels vs computed
print()
print('=== ORIGINAL LABEL vs COMPUTED ===')
from itertools import product
xtab = Counter()
for s in samples:
    orig = 'case' if s['label'] == '1' else 'ctrl' if s['label'] == '0' else 'unk'
    comp = 'CASE' if s['classification'] in ('CASE','CASE_IF_GA') else 'CTRL' if s['classification'] in ('CONTROL','CONTROL_IF_DELIVERED') else 'OTHER'
    xtab[(orig, comp)] += 1
print(f'                  Computed_CASE  Computed_CTRL  Computed_OTHER')
for orig in ['case','ctrl','unk']:
    c = xtab.get((orig,'CASE'),0)
    t = xtab.get((orig,'CTRL'),0)
    o = xtab.get((orig,'OTHER'),0)
    print(f'  Orig_{orig:4s}:      {c:5d}          {t:5d}          {o:5d}')