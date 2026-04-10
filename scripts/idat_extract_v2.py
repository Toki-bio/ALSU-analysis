#!/usr/bin/env python3
"""Extract control probe intensities from IDAT files — fixed format.
IDAT v3: no count prefix on arrays, use nSNPsRead.
Manifest: colon-separated addresses with leading zeros.
"""
import struct
import os
import sys
from collections import defaultdict

BASE = '/staging/ALSU-analysis/Conversion/OUTPUT/'
MANIFEST = '/staging/ALSU-analysis/Conversion/GSA-24v3-0_A2.csv'
BAD = ['208993030034', '208993030080', '208993030112', '208993030109']
GOOD = ['207591080076', '207585740001', '207591070002']

def avg(l):
    return sum(l)/len(l) if l else 0

# ========== Parse manifest controls ==========
controls = []
in_controls = False
with open(MANIFEST, 'r', errors='ignore') as f:
    for line in f:
        s = line.strip()
        if '[Controls]' in s:
            in_controls = True
            continue
        if in_controls:
            if s.startswith('['):
                break
            if not s:
                continue
            parts = s.split(',')
            # First field: colon-separated address like 0027630314:0027630314:...
            addr_field = parts[0].strip()
            addr_parts = addr_field.split(':')
            try:
                addr = int(addr_parts[0].lstrip('0') or '0')
            except:
                continue
            # Remaining fields: Type, Color, Extended_Type
            ctype = parts[1].strip() if len(parts) > 1 else ''
            color = parts[2].strip() if len(parts) > 2 else ''
            ext = parts[3].strip() if len(parts) > 3 else ''
            controls.append({'address': addr, 'type': ctype, 'color': color, 'ext': ext})

print(f"Control probes: {len(controls)}")
for c in controls:
    print(f"  addr={c['address']:>10d}  type={c['type']:<20} color={c['color']:<8} ext={c['ext']}")

# ========== Parse IDAT v3 ==========
def read_idat(path):
    """Parse IDAT v3 file. No count prefix on arrays."""
    with open(path, 'rb') as f:
        magic = f.read(4)
        if magic != b'IDAT':
            return None
        version = struct.unpack('<q', f.read(8))[0]
        nf = struct.unpack('<i', f.read(4))[0]
        ft = {}
        for _ in range(nf):
            c = struct.unpack('<H', f.read(2))[0]
            o = struct.unpack('<q', f.read(8))[0]
            ft[c] = o
        # nSNPsRead
        f.seek(ft[1000])
        n = struct.unpack('<i', f.read(4))[0]
        # IlluminaID — no count prefix
        f.seek(ft[102])
        ids = struct.unpack(f'<{n}i', f.read(4*n))
        # Mean — no count prefix
        f.seek(ft[104])
        means = struct.unpack(f'<{n}H', f.read(2*n))
        return dict(zip(ids, means))

# ========== Check control addresses exist in IDAT ==========
test_path = os.path.join(BASE, GOOD[0], f'{GOOD[0]}_R01C01_Grn.idat')
test_data = read_idat(test_path)
if test_data:
    ctrl_addrs = set(c['address'] for c in controls)
    found = ctrl_addrs.intersection(test_data.keys())
    print(f"\nTest IDAT: {len(test_data)} probes, control match: {len(found)}/{len(ctrl_addrs)}")
    if len(found) < len(ctrl_addrs):
        missing = sorted(ctrl_addrs - set(test_data.keys()))
        print(f"  Missing: {missing}")
        # Show range of IDAT IDs
        all_ids = sorted(test_data.keys())
        print(f"  IDAT ID range: {all_ids[0]} to {all_ids[-1]}")
        print(f"  Control addr range: {min(ctrl_addrs)} to {max(ctrl_addrs)}")

# ========== Parse all chips ==========
print("\n" + "=" * 100)
print("PARSING ALL CHIPS (position R01C01)")
print("=" * 100)
chip_data = {}
for bc in BAD + GOOD:
    chip_dir = os.path.join(BASE, bc)
    if not os.path.isdir(chip_dir):
        print(f"  SKIP {bc}: no dir")
        continue
    grn = os.path.join(chip_dir, f'{bc}_R01C01_Grn.idat')
    red = os.path.join(chip_dir, f'{bc}_R01C01_Red.idat')
    if not os.path.exists(grn) or not os.path.exists(red):
        print(f"  SKIP {bc}: no R01C01 IDATs")
        continue
    try:
        gd = read_idat(grn)
        rd = read_idat(red)
        if gd and rd:
            chip_data[bc] = {'grn': gd, 'red': rd}
            tag = 'BAD' if bc in BAD else 'GOOD'
            print(f"  OK {bc} ({tag:4}): {len(gd)} probes")
        else:
            print(f"  FAIL {bc}")
    except Exception as e:
        print(f"  ERROR {bc}: {e}")

# ========== Control Probe Table ==========
chips = [bc for bc in BAD + GOOD if bc in chip_data]
if not chips:
    print("No chips parsed. Exiting.")
    sys.exit(1)

print("\n" + "=" * 140)
print("CONTROL PROBE INTENSITIES (R01C01, 1 sample per chip)")
print("=" * 140)
hdr = f"{'Type':<20} {'Ext':<22} {'Col':<7}"
for bc in chips:
    tag = 'B' if bc in BAD else 'G'
    hdr += f" | {bc[-4:]}{tag}:Grn  Red"
print(hdr)
print("-" * len(hdr))

for cp in controls:
    a = cp['address']
    row = f"{cp['type']:<20} {cp['ext']:<22} {cp['color']:<7}"
    for bc in chips:
        g = chip_data[bc]['grn'].get(a, -1)
        r = chip_data[bc]['red'].get(a, -1)
        row += f" | {g:6d} {r:5d}"
    print(row)

# ========== Summary by Type ==========
print("\n" + "=" * 110)
print("MEAN INTENSITY BY CONTROL TYPE: BAD vs GOOD")
print("=" * 110)

type_agg = defaultdict(lambda: defaultdict(lambda: {'grn': [], 'red': []}))
for cp in controls:
    a = cp['address']
    ct = cp['type']
    for bc in BAD:
        if bc in chip_data:
            g = chip_data[bc]['grn'].get(a)
            r = chip_data[bc]['red'].get(a)
            if g is not None: type_agg[ct]['BAD']['grn'].append(g)
            if r is not None: type_agg[ct]['BAD']['red'].append(r)
    for bc in GOOD:
        if bc in chip_data:
            g = chip_data[bc]['grn'].get(a)
            r = chip_data[bc]['red'].get(a)
            if g is not None: type_agg[ct]['GOOD']['grn'].append(g)
            if r is not None: type_agg[ct]['GOOD']['red'].append(r)

print(f"{'Control Type':<25} | {'BAD_Grn':>9} {'BAD_Red':>9} | {'GOOD_Grn':>9} {'GOOD_Red':>9} | {'B/G Grn':>8} {'B/G Red':>8}")
print("-" * 110)
for ct in sorted(type_agg):
    bg = avg(type_agg[ct]['BAD']['grn']); br = avg(type_agg[ct]['BAD']['red'])
    gg = avg(type_agg[ct]['GOOD']['grn']); gr = avg(type_agg[ct]['GOOD']['red'])
    rg = bg/gg if gg > 0 else 0; rr = br/gr if gr > 0 else 0
    print(f"{ct:<25} | {bg:9.0f} {br:9.0f} | {gg:9.0f} {gr:9.0f} | {rg:8.3f} {rr:8.3f}")

# ========== Staining Detail (most diagnostic) ==========
print("\n" + "=" * 80)
print("STAINING CONTROLS — DNP = Green channel dye, Biotin = Red channel dye")
print("=" * 80)
staining = [cp for cp in controls if cp['type'] == 'Staining']
for cp in staining:
    a = cp['address']
    print(f"\n  {cp['ext']} (color={cp['color']}, addr={a}):")
    for bc in BAD + GOOD:
        if bc not in chip_data: continue
        tag = 'BAD' if bc in BAD else 'GOOD'
        g = chip_data[bc]['grn'].get(a, -1)
        r = chip_data[bc]['red'].get(a, -1)
        print(f"    {bc} ({tag:4}): Grn={g:6d}  Red={r:6d}")

# ========== Per-chip diagnostics ==========
print("\n" + "=" * 80)
print("PER-CHIP DIAGNOSTIC SUMMARY")
print("=" * 80)
for dtype in ['Staining', 'Hybridization', 'Extension', 'Non-Polymorphic',
              'Target Removal', 'Stringency', 'Non-Specific Binding']:
    probes = [cp for cp in controls if cp['type'] == dtype]
    if not probes: continue
    print(f"\n--- {dtype} ({len(probes)} probes: {', '.join(cp['ext'] for cp in probes)}) ---")
    for bc in BAD + GOOD:
        if bc not in chip_data: continue
        tag = 'BAD' if bc in BAD else 'GOOD'
        gv = [chip_data[bc]['grn'].get(cp['address'], 0) for cp in probes]
        rv = [chip_data[bc]['red'].get(cp['address'], 0) for cp in probes]
        print(f"  {bc} ({tag:4}): Grn avg={avg(gv):7.0f} max={max(gv):6d} min={min(gv):5d}  "
              f"Red avg={avg(rv):7.0f} max={max(rv):6d} min={min(rv):5d}")

# ========== Also: overall probe intensity distribution (top 100 probes) ==========
print("\n" + "=" * 80)
print("OVERALL PROBE INTENSITY (median/mean/p95 of ALL 704K probes)")
print("=" * 80)
for bc in BAD + GOOD:
    if bc not in chip_data: continue
    tag = 'BAD' if bc in BAD else 'GOOD'
    grn_vals = sorted(chip_data[bc]['grn'].values())
    red_vals = sorted(chip_data[bc]['red'].values())
    n = len(grn_vals)
    g_med = grn_vals[n//2]; g_mean = avg(grn_vals); g_p95 = grn_vals[int(n*0.95)]
    r_med = red_vals[n//2]; r_mean = avg(red_vals); r_p95 = red_vals[int(n*0.95)]
    print(f"  {bc} ({tag:4}): Grn median={g_med:5d} mean={g_mean:7.0f} p95={g_p95:6d}  "
          f"Red median={r_med:5d} mean={r_mean:7.0f} p95={r_p95:6d}")

print("\nDONE")
