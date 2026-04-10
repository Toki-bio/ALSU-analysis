#!/usr/bin/env python3
import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

with open('sample_investigation_v2.html', 'r', encoding='utf-8') as f:
    content = f.read()

open_br = content.count('{')
close_br = content.count('}')

print(f'Braces: {open_br} open, {close_br} close')
if open_br != close_br:
    print('ERROR: Braces unbalanced')
else:
    print('OK: Braces balanced')

print()
print('Checks:')
if 'const ES = D.enriched_samples' in content:
    print('  render() function: OK')
print(f'  Plotly charts: {content.count("Plotly.newPlot")} instances')
print(f'  Section 2 mentions "FID ordering": {"yes" if "FID ordering" in content else "NO"}')
print(f'  Section 3 horizontal bars: {"yes" if "orientation: \'h\'" in content else "NO"}')
print(f'  Section 3 height 1000: {"yes" if "height: 1000" in content else "NO"}')
