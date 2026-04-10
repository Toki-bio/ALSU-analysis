"""Parse iScan QC files and compare bad vs good chips quantitatively."""
import re, statistics

text = open('data/deep_chip_qc_utf8.txt', encoding='utf-8').read()

# Extract QC blocks per chip
# Pattern: "--- BARCODE_qc.txt ---" followed by data until next "---" or "==="
chip_sections = re.split(r'--- (\d+)_qc\.txt ---', text)

chips = {}
for i in range(1, len(chip_sections), 2):
    barcode = chip_sections[i]
    data = chip_sections[i+1]
    
    # Parse rows: only Swath lines have all the metrics
    rows = []
    for line in data.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) >= 8 and 'Swath' in parts[0]:
            try:
                rows.append({
                    'section': parts[0],
                    'reg_score': float(parts[1]),
                    'mean_on': float(parts[2]),
                    'mean_off': float(parts[3]),
                    'sep_metric': float(parts[4]),
                    'num_on': float(parts[5]),
                    'snr': float(parts[6]),
                    'cutoff': float(parts[7]),
                    'focus': float(parts[8]) if len(parts) > 8 and parts[8] else None
                })
            except (ValueError, IndexError):
                pass
        if line.startswith('==='): 
            break
    
    if rows:
        chips[barcode] = rows

BAD = ['208993030034', '208993030080', '208993030112', '208993030109']
GOOD = ['207591080076', '207585740001', '207591070002']

print("=" * 80)
print("CHIP SCANNER QC COMPARISON: BAD vs GOOD")
print("=" * 80)
print(f"\n{'Barcode':<15} {'N_rows':>6} {'Mean_ON':>10} {'Mean_OFF':>10} {'S/N':>8} {'Focus':>8} {'Reg':>8} {'Sep':>10} {'Status'}")
print("-" * 95)

for barcode in BAD + GOOD:
    if barcode not in chips:
        print(f"{barcode:<15} NOT FOUND")
        continue
    rows = chips[barcode]
    # Channel1 = Green, Channel2 = Red
    ch1 = [r for r in rows if 'Channel1' in r['section']]
    ch2 = [r for r in rows if 'Channel2' in r['section']]
    all_rows = rows
    
    mean_on = statistics.mean(r['mean_on'] for r in all_rows)
    mean_off = statistics.mean(r['mean_off'] for r in all_rows)
    snr = statistics.mean(r['snr'] for r in all_rows)
    foc_vals = [r['focus'] for r in all_rows if r['focus'] is not None]
    focus = statistics.mean(foc_vals) if foc_vals else 0
    reg = statistics.mean(abs(r['reg_score']) for r in all_rows)
    sep = statistics.mean(r['sep_metric'] for r in all_rows)
    status = 'BAD' if barcode in BAD else 'GOOD'
    
    print(f"{barcode:<15} {len(all_rows):>6} {mean_on:>10.1f} {mean_off:>10.1f} {snr:>8.4f} {focus:>8.4f} {reg:>8.4f} {sep:>10.1f} {status}")

# Detailed position-by-position analysis for worst bad chip
print("\n" + "=" * 80)
print("POSITION-LEVEL DETAIL: 208993030034 (worst chip)")
print("=" * 80)
if '208993030034' in chips:
    rows = chips['208993030034']
    # Group by R position (R01-R12)
    by_row = {}
    for r in rows:
        m = re.match(r'R(\d+)', r['section'])
        if m:
            row_num = int(m.group(1))
            by_row.setdefault(row_num, []).append(r)
    
    print(f"\n{'Row':>4} {'N':>4} {'Mean_ON':>10} {'S/N':>8} {'Focus':>8} {'Reg':>8}")
    print("-" * 50)
    for row_num in sorted(by_row.keys()):
        rrows = by_row[row_num]
        mon = statistics.mean(r['mean_on'] for r in rrows)
        snr = statistics.mean(r['snr'] for r in rrows)
        foc = statistics.mean(r['focus'] for r in rrows if r['focus'] is not None) if any(r['focus'] is not None for r in rrows) else 0
        reg = statistics.mean(r['reg_score'] for r in rrows)
        print(f"R{row_num:02d}  {len(rrows):>4} {mon:>10.1f} {snr:>8.4f} {foc:>8.4f} {reg:>8.4f}")

# Compare control probe regions (positions R12 usually have staining controls)
print("\n" + "=" * 80)
print("KEY DIAGNOSTIC: Mean ON signal comparison")
print("=" * 80)
for barcode in BAD + GOOD:
    if barcode not in chips:
        continue
    rows = chips[barcode]
    ch1 = [r for r in rows if 'Channel1' in r['section']]
    ch2 = [r for r in rows if 'Channel2' in r['section']]
    
    grn_mean = statistics.mean(r['mean_on'] for r in ch1) if ch1 else 0
    red_mean = statistics.mean(r['mean_on'] for r in ch2) if ch2 else 0
    grn_min = min(r['mean_on'] for r in ch1) if ch1 else 0
    red_min = min(r['mean_on'] for r in ch2) if ch2 else 0
    
    status = 'BAD' if barcode in BAD else 'GOOD'
    print(f"{barcode} ({status}): Grn mean={grn_mean:.0f} min={grn_min:.0f}  Red mean={red_mean:.0f} min={red_min:.0f}")

# Scan dates from XML
print("\n" + "=" * 80)
print("SCAN DATES (from XML metadata)")
print("=" * 80)
dates = re.findall(r'<SentrixMatrixID>(\d+)</SentrixMatrixID>.*?<ImageDate>([^<]+)</ImageDate>', text, re.DOTALL)
for barcode, date in dates:
    status = 'BAD' if barcode in BAD else 'GOOD' if barcode in GOOD else '???'
    print(f"{barcode} ({status}): Scanned {date}")

# P95 values from XML  
print("\n" + "=" * 80)
print("P95 INTENSITY (from XML — overall chip brightness)")
print("=" * 80)
p95s = re.findall(r'<SentrixMatrixID>(\d+)</SentrixMatrixID>.*?<P95>([^<]+)</P95>', text, re.DOTALL)
for barcode, p95 in p95s:
    status = 'BAD' if barcode in BAD else 'GOOD' if barcode in GOOD else '???'
    print(f"{barcode} ({status}): P95={float(p95):.1f}")
