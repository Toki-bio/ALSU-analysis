#!/usr/bin/env python3
"""
Build comprehensive investigation data including chip mapping, sex check, and heterozygosity.
All data from real server extracts — no fabrication.
"""
import re
import json
import csv
from collections import defaultdict

DUMP = r"C:\work\ALSU-analysis\data\investigation_dump.txt"
CHECKS = r"C:\work\ALSU-analysis\data\checks_output.txt"
SAMPLESHEET = r"C:\work\ALSU-analysis\data\samplesheet_dump.txt"
IBD_FILE = r"C:\work\ALSU-analysis\data\ibd_results.txt"
OLD_DATA = r"C:\work\ALSU-analysis\data\investigation_data.json"

# ============================================================
# 1. Parse sample sheet for chip mapping
# ============================================================
def parse_sample_sheet(path):
    """Extract Sample_ID -> (SentrixBarcode, SentrixPosition) from the sample sheet dump."""
    mapping = {}
    in_full = False
    in_data = False
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('=== FULL SAMPLE SHEET ==='):
                in_full = True
                continue
            if not in_full:
                continue
            if line.startswith('Sample_ID,SentrixBarcode'):
                in_data = True
                continue
            if in_data and '===' in line:
                # Handle last line concatenated with === DONE ===
                if ',' in line and line.split(',')[0] not in ('', '==='):
                    data_part = line.split('===')[0].strip()
                    if data_part:
                        parts = data_part.split(',')
                        if len(parts) >= 3 and parts[0] and parts[1]:
                            mapping[parts[0].strip()] = {
                                'barcode': parts[1].strip(),
                                'position': parts[2].strip() if len(parts) > 2 else '',
                                'plate': parts[3].strip() if len(parts) > 3 else '',
                                'well': parts[4].strip() if len(parts) > 4 else ''
                            }
                break
            if in_data and line:
                parts = line.split(',')
                if len(parts) >= 3 and parts[0] and parts[1]:
                    sample_id = parts[0].strip()
                    barcode = parts[1].strip()
                    position = parts[2].strip() if len(parts) > 2 else ''
                    # Some later entries have more columns (plate, well, date)
                    plate = parts[3].strip() if len(parts) > 3 else ''
                    well = parts[4].strip() if len(parts) > 4 else ''
                    mapping[sample_id] = {
                        'barcode': barcode,
                        'position': position,
                        'plate': plate,
                        'well': well
                    }
    return mapping

# ============================================================
# 2. Parse sex check
# ============================================================
def parse_sexcheck(path):
    """Parse plink --check-sex output from the multi-section file."""
    result = {}
    in_section = False
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if '=== SEX CHECK ===' in line:
                in_section = True
                continue
            if '=== HET CHECK ===' in line:
                in_section = False
                continue
            if in_section and line:
                parts = line.split()
                if len(parts) >= 6 and parts[0] != 'FID':
                    try:
                        fid = parts[0]
                        iid = parts[1]
                        pedsex = int(parts[2])
                        snpsex = int(parts[3])
                        f_val = float(parts[5])
                        result[iid] = {
                            'fid': fid, 'iid': iid,
                            'pedsex': pedsex, 'snpsex': snpsex,
                            'f': f_val
                        }
                    except (ValueError, IndexError):
                        pass
    return result

# ============================================================
# 3. Parse heterozygosity
# ============================================================
def parse_het(path):
    """Parse plink --het output."""
    result = {}
    in_section = False
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if '=== HET CHECK ===' in line:
                in_section = True
                continue
            if '=== ALL DONE ===' in line:
                break
            if in_section and line:
                parts = line.split()
                if len(parts) >= 6 and parts[0] != 'FID':
                    try:
                        iid = parts[1]
                        o_hom = int(parts[2])
                        e_hom = float(parts[3])
                        n_nm = int(parts[4])
                        f_het = float(parts[5])
                        het_rate = (n_nm - o_hom) / n_nm if n_nm > 0 else 0
                        result[iid] = {
                            'iid': iid, 'o_hom': o_hom, 'e_hom': e_hom,
                            'n_nm': n_nm, 'f_het': f_het, 'het_rate': het_rate
                        }
                    except (ValueError, IndexError):
                        pass
    return result

# ============================================================
# 4. Parse .fam and .imiss from investigation dump
# ============================================================
def parse_sections(path):
    sections = {}
    current = None
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            m = re.match(r'^===(\w+.*?)===$', line)
            if m:
                current = m.group(1)
                sections[current] = []
                continue
            if current:
                sections[current].append(line)
    return sections

def strip_suffix(iid):
    m = re.match(r'^(.+?)[dt]$', iid)
    return m.group(1) if m else iid

def get_suffix(iid):
    if iid.endswith('d'): return 'd'
    if iid.endswith('t'): return 't'
    return ''

def get_prefix(iid):
    base = strip_suffix(iid)
    m = re.match(r'^(\d+)-', base)
    return m.group(1) if m else '??'

# Load old data
with open(OLD_DATA, encoding='utf-8') as f:
    old_data = json.load(f)

sections = parse_sections(DUMP)
chip_map = parse_sample_sheet(SAMPLESHEET)
sexcheck = parse_sexcheck(CHECKS)
het_data = parse_het(CHECKS)

print(f"Sample sheet entries: {len(chip_map)}")
print(f"Sex check entries: {len(sexcheck)}")
print(f"Het entries: {len(het_data)}")

# Build imiss lookup
imiss_raw = {}
for line in sections.get('IMISS_RAW', []):
    parts = line.split()
    if len(parts) >= 6 and parts[0] != 'FID':
        try:
            imiss_raw[parts[1]] = float(parts[5])
        except ValueError:
            pass

# Build fam lookup
fam_data = []
for line in sections.get('FAM', []):
    parts = line.split()
    if len(parts) >= 2:
        fam_data.append({'fid': parts[0], 'iid': parts[1], 'idx': int(parts[0])})

iid_to_fid = {s['iid']: s['fid'] for s in fam_data}

# ============================================================
# ANALYSIS: Chip-level quality
# ============================================================
print("\n" + "=" * 70)
print("CHIP-LEVEL QUALITY ANALYSIS")
print("=" * 70)

# For each sample in fam, find its chip barcode
# Note: d/t suffix samples won't be in the sheet under their d/t name
# They'll be under different IDs or duplicated entries

# Build chip -> [samples] mapping
chip_samples = defaultdict(list)  # barcode -> list of {iid, position, fmiss, het_rate}
unmapped_samples = []

for s in fam_data:
    iid = s['iid']
    base = strip_suffix(iid)
    suf = get_suffix(iid)
    
    fmiss = imiss_raw.get(iid, None)
    het = het_data.get(iid, {}).get('het_rate', None)
    sex = sexcheck.get(iid, {}).get('snpsex', None)
    sex_f = sexcheck.get(iid, {}).get('f', None)
    n_nm = het_data.get(iid, {}).get('n_nm', None)
    
    # Try to find in chip map: direct match first, then base
    chip_info = chip_map.get(iid) or chip_map.get(base)
    
    entry = {
        'iid': iid, 'idx': s['idx'], 'fmiss': fmiss,
        'het_rate': het, 'snpsex': sex, 'sex_f': sex_f,
        'n_nm': n_nm, 'suffix': suf
    }
    
    if chip_info:
        entry['barcode'] = chip_info['barcode']
        entry['position'] = chip_info['position']
        entry['plate'] = chip_info.get('plate', '')
        entry['well'] = chip_info.get('well', '')
        chip_samples[chip_info['barcode']].append(entry)
    else:
        entry['barcode'] = 'UNMAPPED'
        entry['position'] = ''
        unmapped_samples.append(entry)

print(f"\nMapped to chips: {sum(len(v) for v in chip_samples.values())}")
print(f"Unmapped: {len(unmapped_samples)}")
if unmapped_samples:
    print(f"  Unmapped IDs: {[s['iid'] for s in unmapped_samples[:20]]}")

# Per-chip statistics
chip_stats = []
for barcode, samples in sorted(chip_samples.items()):
    fmiss_vals = [s['fmiss'] for s in samples if s['fmiss'] is not None]
    het_vals = [s['het_rate'] for s in samples if s['het_rate'] is not None]
    n_samples = len(samples)
    mean_fmiss = sum(fmiss_vals) / len(fmiss_vals) if fmiss_vals else 0
    max_fmiss = max(fmiss_vals) if fmiss_vals else 0
    high_miss = sum(1 for f in fmiss_vals if f > 0.20)
    mean_het = sum(het_vals) / len(het_vals) if het_vals else 0
    
    chip_stats.append({
        'barcode': barcode,
        'n_samples': n_samples,
        'mean_fmiss': mean_fmiss,
        'max_fmiss': max_fmiss,
        'high_miss_count': high_miss,
        'mean_het': mean_het,
        'samples': samples
    })

chip_stats.sort(key=lambda x: -x['mean_fmiss'])

print(f"\n{'Barcode':<18} {'N':>3} {'Mean F_MISS':>11} {'Max F_MISS':>10} {'High(>20%)':>10} {'Mean Het':>9}")
print("-" * 65)
for cs in chip_stats:
    print(f"{cs['barcode']:<18} {cs['n_samples']:>3} {cs['mean_fmiss']:>11.4f} {cs['max_fmiss']:>10.4f} {cs['high_miss_count']:>10} {cs['mean_het']:>9.4f}")

# ============================================================
# ANALYSIS: Sex check for d/t pair cross-ref
# ============================================================
print("\n" + "=" * 70)
print("SEX CHECK: d/t PAIRS CROSS-REFERENCE")
print("=" * 70)

dt_samples = [s for s in fam_data if get_suffix(s['iid']) in ('d', 't')]
sex_mismatches = []

for dt in dt_samples:
    iid = dt['iid']
    base = strip_suffix(iid)
    sex_dt = sexcheck.get(iid, {}).get('snpsex', None)
    sex_base = sexcheck.get(base, {}).get('snpsex', None)
    
    if sex_dt is not None and sex_base is not None:
        match = 'MATCH' if sex_dt == sex_base else 'MISMATCH'
        if sex_dt != sex_base:
            sex_mismatches.append(iid)
        sex_str = lambda x: {0: 'unknown', 1: 'male', 2: 'female'}.get(x, '?')
        print(f"  {iid:>12} ({sex_str(sex_dt)}) vs {base:>12} ({sex_str(sex_base)}) => {match}")

print(f"\nSex mismatches between d/t and base: {len(sex_mismatches)}")

# ============================================================
# ANALYSIS: Hyper-connected samples heterozygosity
# ============================================================
print("\n" + "=" * 70)
print("HYPER-CONNECTED SAMPLES — HETEROZYGOSITY CHECK")
print("=" * 70)

hyper_ids = ['08-495', '08-123', '20-05', '08-25', '08-86', '07-04d', '08-701']
all_het_rates = [v['het_rate'] for v in het_data.values() if v['het_rate']]
mean_global_het = sum(all_het_rates) / len(all_het_rates)
std_het = (sum((h - mean_global_het)**2 for h in all_het_rates) / len(all_het_rates)) ** 0.5

print(f"Global mean het rate: {mean_global_het:.4f} +/- {std_het:.4f}")
for iid in hyper_ids:
    h = het_data.get(iid, {})
    if h:
        z = (h['het_rate'] - mean_global_het) / std_het if std_het > 0 else 0
        fmiss = imiss_raw.get(iid, 0)
        print(f"  {iid:>10}: het_rate={h['het_rate']:.4f} (z={z:+.2f}), F={h['f_het']:+.4f}, N_NM={h['n_nm']}, F_MISS={fmiss:.4f}")

# ============================================================
# ANALYSIS: Unexpected identical pairs — check chip positions
# ============================================================
print("\n" + "=" * 70)
print("IDENTICAL PAIRS: CHIP POSITION ANALYSIS")
print("=" * 70)

# Load pairs from old data
pairs_098 = old_data['pairs_098']

for p in pairs_098:
    iid1, iid2 = p['iid1'], p['iid2']
    base1, base2 = strip_suffix(iid1), strip_suffix(iid2)
    
    if base1 == base2:
        continue  # expected
    
    chip1 = chip_map.get(iid1) or chip_map.get(base1)
    chip2 = chip_map.get(iid2) or chip_map.get(base2)
    
    bc1 = chip1['barcode'] if chip1 else '?'
    pos1 = chip1['position'] if chip1 else '?'
    bc2 = chip2['barcode'] if chip2 else '?'
    pos2 = chip2['position'] if chip2 else '?'
    
    same_chip = 'SAME CHIP' if bc1 == bc2 else 'diff chip'
    
    print(f"  {iid1:>12} ({bc1} {pos1}) <=> {iid2:>12} ({bc2} {pos2})  {same_chip}  PI_HAT={p['pihat']:.4f}")

# ============================================================
# BUILD OUTPUT JSON for the HTML page
# ============================================================
print("\n\nBuilding output JSON...")

# Enrich fam data with all fields
enriched_samples = []
for s in fam_data:
    iid = s['iid']
    base = strip_suffix(iid)
    chip_info = chip_map.get(iid) or chip_map.get(base) or {}
    
    enriched_samples.append({
        'iid': iid,
        'idx': s['idx'],
        'fmiss': imiss_raw.get(iid),
        'het_rate': het_data.get(iid, {}).get('het_rate'),
        'f_het': het_data.get(iid, {}).get('f_het'),
        'n_nm': het_data.get(iid, {}).get('n_nm'),
        'snpsex': sexcheck.get(iid, {}).get('snpsex'),
        'sex_f': sexcheck.get(iid, {}).get('f'),
        'barcode': chip_info.get('barcode', ''),
        'position': chip_info.get('position', ''),
        'suffix': get_suffix(iid),
        'prefix': get_prefix(iid)
    })

# Chip stats for chart
chip_stats_out = [{
    'barcode': cs['barcode'],
    'n_samples': cs['n_samples'],
    'mean_fmiss': cs['mean_fmiss'],
    'max_fmiss': cs['max_fmiss'],
    'high_miss_count': cs['high_miss_count'],
    'mean_het': cs['mean_het']
} for cs in chip_stats]

# Identical pair chip analysis
pair_chip = []
for p in pairs_098:
    iid1, iid2 = p['iid1'], p['iid2']
    base1, base2 = strip_suffix(iid1), strip_suffix(iid2)
    
    chip1 = chip_map.get(iid1) or chip_map.get(base1) or {}
    chip2 = chip_map.get(iid2) or chip_map.get(base2) or {}
    
    pair_chip.append({
        'iid1': iid1, 'iid2': iid2, 'pihat': p['pihat'],
        'bc1': chip1.get('barcode', '?'), 'pos1': chip1.get('position', '?'),
        'bc2': chip2.get('barcode', '?'), 'pos2': chip2.get('position', '?'),
        'same_base': base1 == base2,
        'same_chip': chip1.get('barcode') == chip2.get('barcode')
    })

# Sex mismatch data for d/t
dt_sex = []
for dt in dt_samples:
    iid = dt['iid']
    base = strip_suffix(iid)
    sex_dt = sexcheck.get(iid, {})
    sex_base = sexcheck.get(base, {})
    dt_sex.append({
        'iid': iid, 'base': base,
        'snpsex_dt': sex_dt.get('snpsex'),
        'snpsex_base': sex_base.get('snpsex'),
        'f_dt': sex_dt.get('f'),
        'f_base': sex_base.get('f'),
        'match': sex_dt.get('snpsex') == sex_base.get('snpsex') if sex_dt.get('snpsex') and sex_base.get('snpsex') else None
    })

output = {
    **old_data,  # keep existing data
    'enriched_samples': enriched_samples,
    'chip_stats': chip_stats_out,
    'pair_chip': pair_chip,
    'dt_sex': dt_sex,
    'sex_mismatches': sex_mismatches,
    'global_het_mean': mean_global_het,
    'global_het_std': std_het,
    'hyper_connected': [{
        'iid': iid,
        'het_rate': het_data.get(iid, {}).get('het_rate'),
        'f_het': het_data.get(iid, {}).get('f_het'),
        'n_nm': het_data.get(iid, {}).get('n_nm'),
        'fmiss': imiss_raw.get(iid)
    } for iid in hyper_ids]
}

with open(r"C:\work\ALSU-analysis\data\investigation_data_v2.json", 'w') as f:
    json.dump(output, f)

print(f"Wrote investigation_data_v2.json")
print(f"  enriched_samples: {len(enriched_samples)}")
print(f"  chip_stats: {len(chip_stats_out)}")
print(f"  pair_chip: {len(pair_chip)}")
print(f"  dt_sex: {len(dt_sex)}")
