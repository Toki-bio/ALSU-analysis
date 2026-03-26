#!/usr/bin/env python3
"""Update step10.html: Replace ADMIX_DATA with K=2-8, enable buttons, complete CV table, remove banners."""
import re

with open("steps/step10.html", encoding="utf-8") as f:
    html = f.read()

original_len = len(html)

# 1. Replace ADMIX_DATA line (line starting with 'const ADMIX_DATA =')
with open("data/admix_v2_jsline_full.txt", encoding="utf-8") as f:
    new_js_line = f.read().strip()

html = re.sub(
    r'^const ADMIX_DATA = \{.*?\};$',
    new_js_line,
    html,
    count=1,
    flags=re.MULTILINE
)
print("1. Replaced ADMIX_DATA blob")

# 2. Update comment above ADMIX_DATA
html = html.replace(
    "// ADMIXTURE Q-matrix data (1,047 samples x K=2-5; K=6-8 pending)",
    "// ADMIXTURE Q-matrix data (1,047 samples × K=2-8)"
)
print("2. Updated ADMIX_DATA comment")

# 3. Remove the "V2 re-analysis in progress" info-box
html = re.sub(
    r'<div class="info-box"[^>]*>\s*<strong>&#9888; V2 re-analysis in progress\.</strong>.*?</div>',
    '',
    html,
    count=1,
    flags=re.DOTALL
)
print("3. Removed V2 in-progress banner")

# 4. Update K values line
html = html.replace(
    "K=2 through K=8 (K=2&ndash;5 complete; K=6&ndash;8 running)",
    "K=2 through K=8"
)
print("4. Updated K values description")

# 5. Enable K=6, K=7, K=8 buttons
html = html.replace(
    '<button disabled style="opacity:0.4" title="Running on server">K = 6</button>',
    '<button onclick="drawAdmixture(6)">K = 6</button>'
)
html = html.replace(
    '<button disabled style="opacity:0.4" title="Pending">K = 7</button>',
    '<button onclick="drawAdmixture(7)">K = 7</button>'
)
html = html.replace(
    '<button disabled style="opacity:0.4" title="Pending">K = 8</button>',
    '<button onclick="drawAdmixture(8)">K = 8</button>'
)
print("5. Enabled K=6,7,8 buttons")

# 6. Replace K=6 Running row
html = html.replace(
    '''    <tr style="color: #999;">
        <td>K=6</td>
        <td colspan="3"><em>Running&hellip;</em></td>
    </tr>''',
    '''    <tr>
        <td><strong>K=6</strong></td>
        <td>0.31266</td>
        <td>+0.00504</td>
        <td style="color:#c62828;">Worse</td>
    </tr>'''
)
print("6. Updated K=6 CV row")

# 7. Replace K=7 Pending row
html = html.replace(
    '''    <tr style="color: #999;">
        <td>K=7</td>
        <td colspan="3"><em>Pending</em></td>
    </tr>''',
    '''    <tr>
        <td><strong>K=7</strong></td>
        <td>0.31445</td>
        <td>+0.00683</td>
        <td style="color:#c62828;">Worse</td>
    </tr>'''
)
print("7. Updated K=7 CV row")

# 8. Replace K=8 Pending row
html = html.replace(
    '''    <tr style="color: #999;">
        <td>K=8</td>
        <td colspan="3"><em>Pending</em></td>
    </tr>''',
    '''    <tr>
        <td><strong>K=8</strong></td>
        <td>0.31627</td>
        <td>+0.00865</td>
        <td style="color:#c62828;">Worse</td>
    </tr>'''
)
print("8. Updated K=8 CV row")

# 9. Update the interpretation paragraph
html = html.replace(
    "CV error increases <strong>monotonically</strong> from K=2 (0.30762) through K=5 (0.31103), and the remaining K=6&ndash;8 runs are expected to continue upward.",
    "CV error increases <strong>monotonically</strong> from K=2 (0.30762) through K=8 (0.31627)."
)
print("9. Updated interpretation paragraph")

with open("steps/step10.html", "w", encoding="utf-8") as f:
    f.write(html)

new_len = len(html)
print(f"\nDone! {original_len} -> {new_len} chars ({new_len - original_len:+d})")
