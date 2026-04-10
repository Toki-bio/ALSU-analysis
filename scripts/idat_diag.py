#!/usr/bin/env python3
"""Diagnostic + extraction: understand IDAT format and manifest controls, then extract."""
import struct
import os
import sys
from collections import defaultdict

BASE = '/staging/ALSU-analysis/Conversion/OUTPUT/'
MANIFEST = '/staging/ALSU-analysis/Conversion/GSA-24v3-0_A2.csv'
BAD = ['208993030034', '208993030080', '208993030112', '208993030109']
GOOD = ['207591080076', '207585740001', '207591070002']

# ===================== MANIFEST =====================
print("=" * 80)
print("MANIFEST [Controls] SECTION")
print("=" * 80)
in_controls = False
ctrl_lines = []
with open(MANIFEST, 'r', errors='ignore') as f:
    for line in f:
        if '[Controls]' in line:
            in_controls = True
            print(f"HEADER: {line.rstrip()}")
            continue
        if in_controls:
            if line.strip().startswith('['):
                print(f"END at: {line.rstrip()}")
                break
            if line.strip():
                ctrl_lines.append(line.strip())

print(f"Total control lines: {len(ctrl_lines)}")
for i, ln in enumerate(ctrl_lines[:8]):
    print(f"  [{i}] {ln}")
if len(ctrl_lines) > 8:
    print(f"  ... ({len(ctrl_lines) - 8} more)")

# Parse controls - try different column formats
controls = []
for ln in ctrl_lines:
    parts = ln.split(',')
    # typical Illumina format: Index, Address, Type, Extended_Type, Color, ...
    # or: Index, Type, Extended_Type, Address, ...
    # Find the large integer (address) and the type string
    addr = None
    ctype = ''
    ext = ''
    known_types = ['Staining', 'Extension', 'Target Removal', 'Hybridization',
                   'Stringency', 'Non-Specific Binding', 'Non-Polymorphic',
                   'Restoration', 'Negative']
    for i, p in enumerate(parts):
        p = p.strip()
        if p.isdigit() and len(p) >= 5:  # addresses are typically 7-8 digits
            if addr is None:
                addr = int(p)
        if p in known_types:
            ctype = p
            if i + 1 < len(parts):
                ext = parts[i + 1].strip()
    if addr:
        controls.append({'address': addr, 'type': ctype, 'ext': ext, 'raw': ln[:60]})

print(f"\nParsed {len(controls)} control probes")
if controls:
    print(f"  Types: {sorted(set(c['type'] for c in controls if c['type']))}")
    print(f"  First 3:")
    for c in controls[:3]:
        print(f"    addr={c['address']} type={c['type']} ext={c['ext']}")
        print(f"    raw: {c['raw']}")

# ===================== IDAT FORMAT =====================
print("\n" + "=" * 80)
print("IDAT FORMAT DISCOVERY")
print("=" * 80)

chip = GOOD[0]
idat_path = os.path.join(BASE, chip, f'{chip}_R01C01_Grn.idat')
print(f"Test file: {idat_path}")
print(f"File size: {os.path.getsize(idat_path)} bytes ({os.path.getsize(idat_path)/1024/1024:.1f} MB)")

with open(idat_path, 'rb') as f:
    magic = f.read(4)
    version = struct.unpack('<q', f.read(8))[0]
    nfields = struct.unpack('<i', f.read(4))[0]
    print(f"Magic: {magic}, Version: {version}, Fields: {nfields}")

    ftable = {}
    for _ in range(nfields):
        code = struct.unpack('<H', f.read(2))[0]
        offset = struct.unpack('<q', f.read(8))[0]
        ftable[code] = offset

    print("Field table:")
    for code in sorted(ftable.keys()):
        print(f"  Code {code:5d} -> offset {ftable[code]:12d}")

    # nSNPsRead
    if 1000 in ftable:
        f.seek(ftable[1000])
        nsnps = struct.unpack('<i', f.read(4))[0]
        print(f"\nnSNPsRead = {nsnps}")

    # Probe field 102 (IlluminaID)
    if 102 in ftable:
        f.seek(ftable[102])
        raw20 = f.read(20)
        print(f"\nField 102 at {ftable[102]}: hex={raw20.hex()}")
        # Interpret as: count_prefix (int32) + int32_array
        cnt = struct.unpack_from('<i', raw20, 0)[0]
        print(f"  As [count_prefix={cnt}] + int32[]: {struct.unpack_from('<4i', raw20, 4)}")
        # Interpret as: direct int32[] (no prefix, use nsnps)
        print(f"  As direct int32[] (no prefix): {struct.unpack_from('<5i', raw20, 0)}")

    # Mean field 104
    if 104 in ftable:
        f.seek(ftable[104])
        raw20 = f.read(20)
        print(f"\nField 104 at {ftable[104]}: hex={raw20.hex()}")
        cnt = struct.unpack_from('<i', raw20, 0)[0]
        print(f"  As [count_prefix={cnt}] + uint16[]: {struct.unpack_from('<8H', raw20, 4)}")
        print(f"  As direct uint16[] (no prefix): {struct.unpack_from('<10H', raw20, 0)}")

    # SD field 103
    if 103 in ftable:
        f.seek(ftable[103])
        raw20 = f.read(20)
        print(f"\nField 103 at {ftable[103]}: hex={raw20.hex()}")

    # NBeads field 107
    if 107 in ftable:
        f.seek(ftable[107])
        raw20 = f.read(20)
        print(f"\nField 107 at {ftable[107]}: hex={raw20.hex()}")

    # String fields (barcode etc)
    for code, name in [(300, 'RunInfo'), (400, 'RedGreen'), (402, 'Barcode'), (403, 'ChipType'), (406, 'SampleID')]:
        if code in ftable:
            f.seek(ftable[code])
            raw = f.read(40)
            print(f"\nField {code} ({name}) at {ftable[code]}: hex={raw[:20].hex()} ascii={raw[:40]}")


# ===================== PARSE IDATs WITH CORRECT FORMAT =====================
def read_idat_v3(path, nsnps_expected):
    """Parse IDAT using nSNPsRead for array lengths (no count prefix)."""
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

        # Read nSNPsRead
        f.seek(ft[1000])
        nsnps = struct.unpack('<i', f.read(4))[0]

        # Try both formats: with count prefix and without
        # IlluminaID (code 102)
        f.seek(ft[102])
        first4 = struct.unpack('<i', f.read(4))[0]

        if first4 == nsnps:
            # Has count prefix
            ids = list(struct.unpack(f'<{nsnps}i', f.read(4 * nsnps)))
            f.seek(ft[104])
            n2 = struct.unpack('<i', f.read(4))[0]
            means = list(struct.unpack(f'<{n2}H', f.read(2 * n2)))
        else:
            # No count prefix - read nsnps values starting from offset
            f.seek(ft[102])
            ids = list(struct.unpack(f'<{nsnps}i', f.read(4 * nsnps)))
            f.seek(ft[104])
            means = list(struct.unpack(f'<{nsnps}H', f.read(2 * nsnps)))

        return dict(zip(ids, means))


# Get nsnps from test file
print("\n\n" + "=" * 80)
print("EXTRACTING CONTROL PROBE DATA")
print("=" * 80)

if not controls:
    print("No control probes parsed from manifest. Cannot proceed.")
    print("\nLet me try: dump first 5 lines after [Controls] as raw bytes:")
    for i, ln in enumerate(ctrl_lines[:5]):
        print(f"  {i}: {repr(ln)}")
    sys.exit(0)

# Parse all chips
chip_data = {}
for bc in BAD + GOOD:
    chip_dir = os.path.join(BASE, bc)
    if not os.path.isdir(chip_dir):
        continue
    grn = os.path.join(chip_dir, f'{bc}_R01C01_Grn.idat')
    red = os.path.join(chip_dir, f'{bc}_R01C01_Red.idat')
    if not os.path.exists(grn) or not os.path.exists(red):
        print(f"  SKIP {bc}: R01C01 IDATs missing")
        continue
    try:
        gd = read_idat_v3(grn, nsnps)
        rd = read_idat_v3(red, nsnps)
        if gd and rd:
            chip_data[bc] = {'grn': gd, 'red': rd}
            tag = 'BAD' if bc in BAD else 'GOOD'
            # Check how many control addresses are present
            ctrl_addrs = set(c['address'] for c in controls)
            found_g = len(ctrl_addrs.intersection(gd.keys()))
            found_r = len(ctrl_addrs.intersection(rd.keys()))
            print(f"  {bc} ({tag:4}): {len(gd)} probes, ctrl_match: Grn={found_g} Red={found_r}")
        else:
            print(f"  FAIL {bc}")
    except Exception as e:
        print(f"  ERROR {bc}: {e}")

if not chip_data:
    print("No chips parsed successfully.")
    sys.exit(0)

# Output control probe table
print("\n" + "=" * 120)
print("CONTROL PROBE INTENSITIES (R01C01)")
print("=" * 120)

def avg(l):
    return sum(l)/len(l) if l else 0

chips = [bc for bc in BAD + GOOD if bc in chip_data]
hdr = f"{'Type':<18} {'Detail':<22}"
for bc in chips:
    tag = 'B' if bc in BAD else 'G'
    hdr += f" |{bc[-4:]}{tag} Grn  Red"
print(hdr)
print("-" * len(hdr))

for cp in controls:
    a = cp['address']
    row = f"{cp['type']:<18} {cp['ext']:<22}"
    for bc in chips:
        g = chip_data[bc]['grn'].get(a, -1)
        r = chip_data[bc]['red'].get(a, -1)
        row += f" |{g:5d} {r:5d}"
    print(row)

# Summary by type
print("\n" + "=" * 110)
print("MEAN INTENSITY BY CONTROL TYPE: BAD vs GOOD")
print("=" * 110)
type_agg = defaultdict(lambda: defaultdict(lambda: {'grn': [], 'red': []}))
for cp in controls:
    a = cp['address']
    ct = cp['type'] or 'Unknown'
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

print(f"{'Control Type':<25} | {'BAD_Grn':>8} {'BAD_Red':>8} | {'GOOD_Grn':>8} {'GOOD_Red':>8} | {'B/G_Grn':>8} {'B/G_Red':>8}")
print("-" * 110)
for ct in sorted(type_agg):
    bg = avg(type_agg[ct]['BAD']['grn']); br = avg(type_agg[ct]['BAD']['red'])
    gg = avg(type_agg[ct]['GOOD']['grn']); gr = avg(type_agg[ct]['GOOD']['red'])
    rg = bg/gg if gg > 0 else 0; rr = br/gr if gr > 0 else 0
    print(f"{ct:<25} | {bg:8.0f} {br:8.0f} | {gg:8.0f} {gr:8.0f} | {rg:8.3f} {rr:8.3f}")

# Staining detail (most diagnostic)
print("\n" + "=" * 80)
print("STAINING CONTROL DETAIL")
print("=" * 80)
staining = [cp for cp in controls if cp['type'] == 'Staining']
for cp in staining:
    a = cp['address']
    print(f"\n  {cp['ext']} (addr={a}):")
    for bc in BAD + GOOD:
        if bc not in chip_data: continue
        tag = 'BAD' if bc in BAD else 'GOOD'
        g = chip_data[bc]['grn'].get(a, -1)
        r = chip_data[bc]['red'].get(a, -1)
        print(f"    {bc} ({tag:4}): Grn={g:6d}  Red={r:6d}")

# Per-chip summary for key types
print("\n" + "=" * 80)
print("PER-CHIP KEY DIAGNOSTICS")
print("=" * 80)
for dtype in ['Staining', 'Hybridization', 'Extension', 'Non-Polymorphic', 'Target Removal']:
    probes = [cp for cp in controls if cp['type'] == dtype]
    if not probes: continue
    print(f"\n--- {dtype} ({len(probes)} probes) ---")
    for bc in BAD + GOOD:
        if bc not in chip_data: continue
        tag = 'BAD' if bc in BAD else 'GOOD'
        gv = [chip_data[bc]['grn'].get(cp['address'], 0) for cp in probes]
        rv = [chip_data[bc]['red'].get(cp['address'], 0) for cp in probes]
        print(f"  {bc} ({tag:4}): Grn mean={avg(gv):7.0f} max={max(gv):6d}  Red mean={avg(rv):7.0f} max={max(rv):6d}")

print("\nDONE")
