"""Extract per-sample detail lists for GRAY, REVIEW, EXCLUDED groups."""
import csv, json
from collections import defaultdict

with open('GWAS от 27.08 - Sheet1.csv', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

SLOT_START = 28; SLOT_SIZE = 7; MAX_SLOTS = 20
LOSS_CODES = {'2','3','6'}
CODE_NAMES = {
    '1':'Live','2':'Stillbirth','3':'Miscarriage',
    '4':'Abortion','5':'Ectopic','6':'Non-dev',
    '7':'Med-abort','8':'Pregnant'
}

details = {'GRAY': [], 'REVIEW': [], 'EXCLUDE_AGE': [], 'EXCLUDE_ECTOPIC': [], 'EXCLUDE_MEDICAL': [], 'UNCLASSIFIED': []}

for row in rows:
    sid = row[0].strip()
    if not sid: continue
    age_raw = row[7].strip()
    try: age = int(age_raw)
    except: age = None

    events = []
    for slot in range(MAX_SLOTS):
        base = SLOT_START + slot * SLOT_SIZE
        if base + 4 >= len(row): break
        outcome = row[base + 3].strip()
        ga_raw = row[base + 4].strip()
        if not outcome or outcome == 'NA': continue
        try: ga = float(ga_raw) if ga_raw else None
        except: ga = None
        events.append({'outcome': outcome, 'ga': ga})

    loss_strict = sum(1 for e in events if e['outcome'] in LOSS_CODES and e['ga'] is not None and e['ga'] < 20)
    loss_any = sum(1 for e in events if e['outcome'] in LOSS_CODES)
    loss_ge20 = sum(1 for e in events if e['outcome'] in LOSS_CODES and e['ga'] is not None and e['ga'] >= 20)
    loss_no_ga = sum(1 for e in events if e['outcome'] in LOSS_CODES and e['ga'] is None)
    live = sum(1 for e in events if e['outcome'] == '1')
    ectopic = sum(1 for e in events if e['outcome'] == '5')
    medical = sum(1 for e in events if e['outcome'] == '7')
    preg_now = sum(1 for e in events if e['outcome'] == '8')
    n_preg = len(events)
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

    if cls not in details:
        continue

    # Build compact event string with GA
    evt_strs = []
    for e in events:
        name = CODE_NAMES.get(e['outcome'], '?')
        ga_str = f" @{int(e['ga'])}w" if e['ga'] is not None else ""
        evt_strs.append(f"{name}{ga_str}")

    rec = {
        'id': sid,
        'age': age,
        'n_preg': n_preg,
        'live': live,
        'loss': loss_any,
        'loss_lt20': loss_strict,
        'loss_ge20': loss_ge20,
        'loss_noGA': loss_no_ga,
        'ectopic': ectopic,
        'events': ' → '.join(evt_strs) if evt_strs else '(none)',
        'seq': ''.join(e['outcome'] for e in events),
        'reason': cls
    }
    details[cls].append(rec)

# Sort each group
for k in details:
    if k == 'GRAY':
        details[k].sort(key=lambda x: (-x['loss'], -x['n_preg'], x['id']))
    elif k == 'REVIEW':
        details[k].sort(key=lambda x: (-x['loss_ge20'], -x['loss'], x['id']))
    else:
        details[k].sort(key=lambda x: x['id'])

# Merge EXCLUDE groups
excluded = details['EXCLUDE_AGE'] + details['EXCLUDE_ECTOPIC'] + details['EXCLUDE_MEDICAL'] + details['UNCLASSIFIED']
# Tag each with sub-reason
for r in details['EXCLUDE_AGE']:
    r['sub'] = f"Age {r['age']}" if r['age'] is not None else "Age unknown"
for r in details['EXCLUDE_ECTOPIC']:
    r['sub'] = 'Ectopic only'
for r in details['EXCLUDE_MEDICAL']:
    r['sub'] = 'Med-abort only'
for r in details['UNCLASSIFIED']:
    r['sub'] = 'No pregnancies' if r['n_preg'] == 0 else 'Unclassified'

output = {
    'gray': details['GRAY'],
    'review': details['REVIEW'],
    'excluded': excluded
}

print(f"GRAY: {len(details['GRAY'])}")
print(f"REVIEW: {len(details['REVIEW'])}")
print(f"EXCLUDED: {len(excluded)} (age:{len(details['EXCLUDE_AGE'])}, ectopic:{len(details['EXCLUDE_ECTOPIC'])}, medical:{len(details['EXCLUDE_MEDICAL'])}, unclass:{len(details['UNCLASSIFIED'])})")

with open('data/sample_details.json', 'w', encoding='utf-8') as f:
    json.dump(output, f, ensure_ascii=False)
print("Wrote data/sample_details.json")
