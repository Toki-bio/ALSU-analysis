#!/usr/bin/env python3
"""Extract control probe intensities from IDAT files for chip QC comparison.
Parses binary IDAT v3 format directly using struct. No external libraries needed.
Compares bad chips vs good reference chips.
"""
import struct
import os
import sys
import re
from collections import defaultdict

def read_idat(path):
    """Parse Illumina IDAT genotyping file. Returns {bead_id: mean_intensity} and metadata dict."""
    with open(path, 'rb') as f:
        magic = f.read(4)
        if magic != b'IDAT':
            return None, {'error': 'not IDAT'}
        version = struct.unpack('<q', f.read(8))[0]
        nfields = struct.unpack('<i', f.read(4))[0]
        ftable = {}
        for _ in range(nfields):
            code = struct.unpack('<H', f.read(2))[0]
            offset = struct.unpack('<q', f.read(8))[0]
            ftable[code] = offset

        def read_int(code):
            f.seek(ftable[code])
            return struct.unpack('<i', f.read(4))[0]

        def read_int_array(code):
            f.seek(ftable[code])
            n = struct.unpack('<i', f.read(4))[0]
            return list(struct.unpack(f'<{n}i', f.read(4 * n)))

        def read_uint16_array(code):
            f.seek(ftable[code])
            n = struct.unpack('<i', f.read(4))[0]
            return list(struct.unpack(f'<{n}H', f.read(2 * n)))

        def read_string(code):
            if code not in ftable:
                return ''
            f.seek(ftable[code])
            # IDAT strings: 1 byte length-type, then data
            # Actually: variable-length encoding
            b = struct.unpack('B', f.read(1))[0]
            if b < 128:
                length = b
            else:
                b2 = struct.unpack('B', f.read(1))[0]
                length = (b & 0x7F) | (b2 << 7)
            return f.read(length).decode('utf-8', errors='replace')

        meta = {'version': version, 'fields': sorted(ftable.keys())}

        # Field 102 = IlluminaID (bead pool addresses), 104 = Mean intensity
        if 102 not in ftable or 104 not in ftable:
            return None, {'error': f'missing fields, available: {sorted(ftable.keys())}'}

        ids = read_int_array(102)
        means = read_uint16_array(104)

        # Optional metadata
        for code, name in [(402, 'barcode'), (403, 'chiptype'), (406, 'sampleid')]:
            try:
                meta[name] = read_string(code)
            except:
                pass

        return dict(zip(ids, means)), meta


def parse_manifest_controls(path):
    """Parse [Controls] section from Illumina manifest CSV."""
    controls = []
    in_controls = False
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            stripped = line.strip()
            if '[Controls]' in stripped:
                in_controls = True
                continue
            if in_controls:
                if stripped.startswith('['):
                    break
                if not stripped:
                    continue
                parts = stripped.split(',')
                # Find address (large integer, typically > 1000000)
                addr = None
                for p in parts:
                    p = p.strip()
                    if p.isdigit() and int(p) > 100000:
                        addr = int(p)
                        break
                if addr is None:
                    continue
                # Find control type
                known_types = ['Staining', 'Extension', 'Target Removal', 'Hybridization',
                               'Stringency', 'Non-Specific Binding', 'Non-Polymorphic',
                               'Restoration', 'Negative']
                ctype = ''
                ext_type = ''
                for i, p in enumerate(parts):
                    p = p.strip()
                    if p in known_types:
                        ctype = p
                        if i + 1 < len(parts):
                            ext_type = parts[i + 1].strip()
                        break
                controls.append({'address': addr, 'type': ctype, 'ext': ext_type})
    return controls


def avg(lst):
    return sum(lst) / len(lst) if lst else 0


# ========================= MAIN =========================
BASE = '/staging/ALSU-analysis/Conversion/OUTPUT/'
BAD = ['208993030034', '208993030080', '208993030112', '208993030109']
GOOD = ['207591080076', '207585740001', '207591070002']

# --- Find manifest ---
manifest = '/staging/ALSU-analysis/Conversion/GSA-24v3-0_A2.csv'
if not os.path.exists(manifest):
    # Search broader
    for d in [BASE, '/staging/ALSU-analysis/Conversion/', os.path.join(BASE, 'ConvSK')]:
        if not os.path.isdir(d):
            continue
        for fn in os.listdir(d):
            if 'GSA-24v3' in fn and fn.endswith('.csv'):
                manifest = os.path.join(d, fn)
                break

if not os.path.exists(manifest):
    print(f"ERROR: Manifest not found at {manifest}")
    sys.exit(1)

print(f"Manifest: {manifest}")
controls = parse_manifest_controls(manifest)
print(f"Control probes: {len(controls)}")
types_found = sorted(set(c['type'] for c in controls if c['type']))
print(f"Types: {types_found}")
if controls:
    print(f"Sample entry: {controls[0]}")
print()

# --- Discover IDAT format from one file ---
test_dir = os.path.join(BASE, GOOD[0])
test_idats = sorted(f for f in os.listdir(test_dir) if f.endswith('.idat'))
print(f"Test chip {GOOD[0]}: {len(test_idats)} IDATs")
print(f"  Names: {test_idats[:4]}")
test_data, test_meta = read_idat(os.path.join(test_dir, test_idats[0]))
if test_data:
    print(f"  Format OK: {len(test_data)} probes, fields={test_meta.get('fields')}")
    ctrl_addrs = set(c['address'] for c in controls)
    found_in_idat = ctrl_addrs.intersection(set(test_data.keys()))
    print(f"  Control addresses matched: {len(found_in_idat)}/{len(ctrl_addrs)}")
    if len(found_in_idat) < len(ctrl_addrs):
        missing = ctrl_addrs - set(test_data.keys())
        print(f"  Missing addresses: {sorted(missing)[:10]}")
else:
    print(f"  PARSE FAILED: {test_meta}")
print()

# --- Parse all chips (position R01C01) ---
print("=" * 100)
print("PARSING CHIPS")
print("=" * 100)
chip_data = {}
for bc in BAD + GOOD:
    chip_dir = os.path.join(BASE, bc)
    if not os.path.isdir(chip_dir):
        print(f"  SKIP {bc}: no directory")
        continue
    idats = sorted(f for f in os.listdir(chip_dir) if f.endswith('.idat'))
    # Find R01C01
    grn_path = red_path = None
    for idat in idats:
        if 'R01C01' in idat:
            if 'Grn' in idat:
                grn_path = os.path.join(chip_dir, idat)
            elif 'Red' in idat:
                red_path = os.path.join(chip_dir, idat)
    if not grn_path or not red_path:
        # Fallback: first available position
        positions = sorted(set(re.search(r'(R\d+C\d+)', f).group(1) for f in idats if re.search(r'(R\d+C\d+)', f)))
        if positions:
            pos = positions[0]
            for idat in idats:
                if pos in idat and 'Grn' in idat:
                    grn_path = os.path.join(chip_dir, idat)
                elif pos in idat and 'Red' in idat:
                    red_path = os.path.join(chip_dir, idat)
            print(f"  {bc}: using {pos} (R01C01 not found)")
    if not grn_path or not red_path:
        print(f"  SKIP {bc}: no IDAT pair")
        continue

    gd, gm = read_idat(grn_path)
    rd, rm = read_idat(red_path)
    if gd and rd:
        chip_data[bc] = {'grn': gd, 'red': rd}
        tag = 'BAD' if bc in BAD else 'GOOD'
        sid = gm.get('sampleid', '?')
        print(f"  {bc} ({tag:4}): Grn={len(gd)} Red={len(rd)} probes  sample={sid}")
    else:
        print(f"  FAIL {bc}: Grn={'ok' if gd else 'fail'} Red={'ok' if rd else 'fail'}")
print()

# --- Control Probe Intensities Table ---
print("=" * 120)
print("CONTROL PROBE INTENSITIES (position R01C01)")
print("=" * 120)
chips_present = [bc for bc in BAD + GOOD if bc in chip_data]
header = f"{'Type':<18} {'Detail':<22}"
for bc in chips_present:
    tag = 'BAD' if bc in BAD else 'OK'
    header += f" | {bc[-6:]}_{tag}  Grn    Red"
print(header)
print("-" * len(header))

for cp in controls:
    addr = cp['address']
    row = f"{cp['type']:<18} {cp['ext']:<22}"
    for bc in chips_present:
        g = chip_data[bc]['grn'].get(addr, -1)
        r = chip_data[bc]['red'].get(addr, -1)
        row += f" | {g:8d} {r:6d}"
    print(row)

# --- Summary by Control Type ---
print("\n" + "=" * 110)
print("MEAN CONTROL INTENSITY BY TYPE")
print("=" * 110)

type_agg = defaultdict(lambda: defaultdict(lambda: {'grn': [], 'red': []}))
for cp in controls:
    addr = cp['address']
    ct = cp['type'] or 'Unknown'
    for bc in BAD:
        if bc in chip_data:
            g = chip_data[bc]['grn'].get(addr)
            r = chip_data[bc]['red'].get(addr)
            if g is not None: type_agg[ct]['BAD']['grn'].append(g)
            if r is not None: type_agg[ct]['BAD']['red'].append(r)
    for bc in GOOD:
        if bc in chip_data:
            g = chip_data[bc]['grn'].get(addr)
            r = chip_data[bc]['red'].get(addr)
            if g is not None: type_agg[ct]['GOOD']['grn'].append(g)
            if r is not None: type_agg[ct]['GOOD']['red'].append(r)

print(f"{'Control Type':<25} | {'BAD_Grn':>8} {'BAD_Red':>8} | {'GOOD_Grn':>8} {'GOOD_Red':>8} | {'B/G Grn':>8} {'B/G Red':>8}")
print("-" * 110)
for ct in sorted(type_agg.keys()):
    bg = avg(type_agg[ct]['BAD']['grn'])
    br = avg(type_agg[ct]['BAD']['red'])
    gg = avg(type_agg[ct]['GOOD']['grn'])
    gr = avg(type_agg[ct]['GOOD']['red'])
    rg = bg / gg if gg > 0 else 0
    rr = br / gr if gr > 0 else 0
    print(f"{ct:<25} | {bg:8.0f} {br:8.0f} | {gg:8.0f} {gr:8.0f} | {rg:8.3f} {rr:8.3f}")

# --- Per-chip detail for key diagnostics ---
print("\n" + "=" * 80)
print("KEY DIAGNOSTIC CONTROLS (per-chip)")
print("=" * 80)
diag_types = ['Staining', 'Hybridization', 'Extension', 'Non-Polymorphic', 'Target Removal', 'Stringency', 'Non-Specific Binding']
for dtype in diag_types:
    probes = [cp for cp in controls if cp['type'] == dtype]
    if not probes:
        continue
    print(f"\n--- {dtype} ({len(probes)} probes) ---")
    # Show extended types
    exts = sorted(set(cp['ext'] for cp in probes))
    print(f"  Sub-types: {exts}")
    for bc in BAD + GOOD:
        if bc not in chip_data:
            continue
        tag = 'BAD' if bc in BAD else 'GOOD'
        gvals = [chip_data[bc]['grn'].get(cp['address'], 0) for cp in probes]
        rvals = [chip_data[bc]['red'].get(cp['address'], 0) for cp in probes]
        print(f"  {bc} ({tag:4}): Grn mean={avg(gvals):7.0f} max={max(gvals):6d} min={min(gvals):5d}  "
              f"Red mean={avg(rvals):7.0f} max={max(rvals):6d} min={min(rvals):5d}")

# --- Staining detail (most diagnostic) ---
print("\n" + "=" * 80)
print("STAINING CONTROL DETAIL (DNP=Green-channel, Biotin=Red-channel)")
print("=" * 80)
staining = [cp for cp in controls if cp['type'] == 'Staining']
for cp in staining:
    addr = cp['address']
    print(f"\n  {cp['ext']} (addr={addr}):")
    for bc in BAD + GOOD:
        if bc not in chip_data:
            continue
        tag = 'BAD' if bc in BAD else 'GOOD'
        g = chip_data[bc]['grn'].get(addr, -1)
        r = chip_data[bc]['red'].get(addr, -1)
        print(f"    {bc} ({tag:4}): Grn={g:6d}  Red={r:6d}")

print("\nDONE")
