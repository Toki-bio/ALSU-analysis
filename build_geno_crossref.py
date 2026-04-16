"""Cross-reference pregnancy detail groups with genotyped sample list,
then re-inject into pre_gwas.html."""
import csv, json, re

# 1. Build genotyped ID set from fam
geno_ids = set()
with open('data/global/global_v2_admix.fam', 'r') as f:
    for line in f:
        parts = line.strip().split()
        if not parts: continue
        iid = parts[1]
        # Skip 1000G IDs (HG/NA prefix)
        if iid.startswith('HG') or iid.startswith('NA'):
            continue
        # Strip numeric prefix: "2_01-02" -> "01-02"
        base = re.sub(r'^\d+_', '', iid)
        geno_ids.add(base)

print(f"Genotyped ALSU IDs: {len(geno_ids)}")

# 2. Load detail data
with open('data/sample_details.json', 'r', encoding='utf-8') as f:
    details = json.load(f)

# 3. Tag each sample and compute stats
stats = {}
for key in ['gray', 'review', 'excluded']:
    n_geno = 0
    for s in details[key]:
        s['genotyped'] = s['id'] in geno_ids
        if s['genotyped']:
            n_geno += 1
    total = len(details[key])
    stats[key] = {'total': total, 'genotyped': n_geno, 'not_geno': total - n_geno}
    print(f"{key.upper()}: {total} total, {n_geno} genotyped, {total - n_geno} not genotyped")

# Also compute for CASE and CONTROL for summary
with open('GWAS от 27.08 - Sheet1.csv', encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    header = next(reader)
    rows = list(reader)

SLOT_START = 28; SLOT_SIZE = 7; MAX_SLOTS = 20
LOSS_CODES = {'2','3','6'}

case_geno = case_total = ctrl_geno = ctrl_total = 0
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

    if only_ectopic: cls = 'EXCLUDE_ECTOPIC'
    elif only_medical and loss_any == 0: cls = 'EXCLUDE_MEDICAL'
    elif not age_ok: cls = 'EXCLUDE_AGE'
    elif loss_ge20 > 0 and loss_strict < 2: cls = 'REVIEW'
    elif loss_strict >= 2: cls = 'CASE'
    elif loss_strict + loss_no_ga >= 2 and loss_strict < 2: cls = 'CASE_PROBABLE'
    elif loss_any == 0 and live >= 1 and not preg_now: cls = 'CONTROL'
    elif loss_any == 0 and live >= 1 and preg_now: cls = 'CONTROL_PROBABLE'
    elif loss_any == 1: cls = 'GRAY'
    elif n_preg == 0: cls = 'UNCLASSIFIED'
    else: cls = 'UNCLASSIFIED'

    is_geno = sid in geno_ids
    if cls == 'CASE':
        case_total += 1
        if is_geno: case_geno += 1
    elif cls == 'CONTROL':
        ctrl_total += 1
        if is_geno: ctrl_geno += 1

print(f"\nCASE: {case_total} total, {case_geno} genotyped")
print(f"CONTROL: {ctrl_total} total, {ctrl_geno} genotyped")
print(f"\nGWAS-ready genotyped: {case_geno} cases × {ctrl_geno} controls = {case_geno + ctrl_geno} samples")

# 4. Save updated details
with open('data/sample_details.json', 'w', encoding='utf-8') as f:
    json.dump(details, f, ensure_ascii=False)

# 5. Build genotype stats object for HTML
geno_stats = {
    'case': {'total': case_total, 'genotyped': case_geno},
    'control': {'total': ctrl_total, 'genotyped': ctrl_geno},
    'gray': stats['gray'],
    'review': stats['review'],
    'excluded': stats['excluded']
}
with open('data/geno_stats.json', 'w', encoding='utf-8') as f:
    json.dump(geno_stats, f, ensure_ascii=False, indent=2)

print(f"\nWrote data/sample_details.json (with genotyped flags)")
print(f"Wrote data/geno_stats.json")
print(json.dumps(geno_stats, indent=2))
