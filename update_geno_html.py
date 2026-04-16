"""Update pre_gwas.html with genotyping status indicators.

Adds:
1. Genotyped counts to summary cards
2. Genotyped counts to expand-panel headers
3. Genotyped column in detail tables
4. Updated GWAS-ready info-box with genotyped counts
5. Updated DETAILS blob with genotyped boolean field
6. CSS for genotyped indicators
"""
import json, re

HTML_PATH = 'pre_gwas.html'
DETAILS_PATH = 'data/sample_details.json'
GENO_PATH = 'data/geno_stats.json'

with open(HTML_PATH, 'r', encoding='utf-8') as f:
    html = f.read()

with open(DETAILS_PATH, 'r', encoding='utf-8-sig') as f:
    details = json.load(f)

with open(GENO_PATH, 'r', encoding='utf-8-sig') as f:
    geno = json.load(f)

# === 1. Add CSS for genotyped indicators ===
geno_css = """
        .geno-yes { color:#2e7d32; font-weight:600; }
        .geno-no { color:#bdbdbd; }
        .geno-badge { display:inline-block; font-size:0.78em; padding:2px 7px; border-radius:10px;
                       margin-left:6px; font-weight:600; }
        .geno-badge-green { background:#e8f5e9; color:#2e7d32; }
        .geno-badge-amber { background:#fff3e0; color:#e65100; }
"""
html = html.replace(
    '        .footer {',
    geno_css + '        .footer {'
)

# === 2. Update summary cards to show genotyped counts ===
# Replace renderSummaryCards function
old_render = """function renderSummaryCards() {
    const c = DATA.summary.classification;
    const cards = [
        {v: DATA.summary.total, l: 'Total Subjects'},
        {v: c.CASE, l: 'Definite Cases'},
        {v: c.CONTROL, l: 'Definite Controls'},
        {v: c.GRAY, l: 'Gray Zone (1 loss)'},
        {v: c.REVIEW, l: 'Review Needed'},
        {v: c.EXCLUDE_AGE + c.EXCLUDE_ECTOPIC + c.EXCLUDE_MEDICAL + c.UNCLASSIFIED, l: 'Excluded'},
    ];
    document.getElementById('summaryCards').innerHTML = cards.map(c =>
        `<div class="stat-card"><div class="stat-value">${c.v.toLocaleString()}</div><div class="stat-label">${c.l}</div></div>`
    ).join('');
}"""

new_render = """function renderSummaryCards() {
    const c = DATA.summary.classification;
    const G = DATA.geno_stats;
    const cards = [
        {v: DATA.summary.total, l: 'Total Subjects', g: G.case.genotyped + G.control.genotyped + G.gray.genotyped + G.review.genotyped + G.excluded.genotyped, gt: DATA.summary.total},
        {v: c.CASE, l: 'Definite Cases', g: G.case.genotyped, gt: G.case.total},
        {v: c.CONTROL, l: 'Definite Controls', g: G.control.genotyped, gt: G.control.total},
        {v: c.GRAY, l: 'Gray Zone (1 loss)', g: G.gray.genotyped, gt: G.gray.total},
        {v: c.REVIEW, l: 'Review Needed', g: G.review.genotyped, gt: G.review.total},
        {v: c.EXCLUDE_AGE + c.EXCLUDE_ECTOPIC + c.EXCLUDE_MEDICAL + c.UNCLASSIFIED, l: 'Excluded', g: G.excluded.genotyped, gt: G.excluded.total},
    ];
    document.getElementById('summaryCards').innerHTML = cards.map(c =>
        `<div class="stat-card"><div class="stat-value">${c.v.toLocaleString()}</div><div class="stat-label">${c.l}</div>` +
        `<div style="font-size:0.75em;color:#666;margin-top:4px;">${c.g} genotyped (${Math.round(c.g/c.gt*100)}%)</div></div>`
    ).join('');
}"""
html = html.replace(old_render, new_render)

# === 3. Update expand-panel headers with genotyped counts ===

# GRAY
html = html.replace(
    '<strong>313 samples</strong>\n            </div>\n            <span class="arrow" id="arrow-gray">&#9654;</span>',
    '<strong>313 samples</strong><span class="geno-badge geno-badge-amber">119 genotyped (38%)</span>\n            </div>\n            <span class="arrow" id="arrow-gray">&#9654;</span>'
)

# REVIEW  
html = html.replace(
    '<strong>116 samples</strong>\n            </div>\n            <span class="arrow" id="arrow-review">&#9654;</span>',
    '<strong>116 samples</strong><span class="geno-badge geno-badge-amber">62 genotyped (53%)</span>\n            </div>\n            <span class="arrow" id="arrow-review">&#9654;</span>'
)

# EXCLUDED
html = html.replace(
    '<strong>94 samples</strong>\n            </div>\n            <span class="arrow" id="arrow-excluded">&#9654;</span>',
    '<strong>94 samples</strong><span class="geno-badge geno-badge-amber">40 genotyped (43%)</span>\n            </div>\n            <span class="arrow" id="arrow-excluded">&#9654;</span>'
)

# === 4. Add "Geno" column to detail table headers ===
html = html.replace(
    '<thead><tr><th>ID</th><th>Age</th><th>Preg</th><th>Live</th><th>Loss</th><th>Pregnancy History</th></tr></thead>',
    '<thead><tr><th>ID</th><th>Geno</th><th>Age</th><th>Preg</th><th>Live</th><th>Loss</th><th>Pregnancy History</th></tr></thead>'
)
html = html.replace(
    '<thead><tr><th>ID</th><th>Age</th><th>Preg</th><th>Live</th><th>Loss</th><th>GA&ge;20</th><th>Pregnancy History</th></tr></thead>',
    '<thead><tr><th>ID</th><th>Geno</th><th>Age</th><th>Preg</th><th>Live</th><th>Loss</th><th>GA&ge;20</th><th>Pregnancy History</th></tr></thead>'
)
html = html.replace(
    '<thead><tr><th>ID</th><th>Age</th><th>Reason</th><th>Preg</th><th>Live</th><th>Loss</th><th>Pregnancy History</th></tr></thead>',
    '<thead><tr><th>ID</th><th>Geno</th><th>Age</th><th>Reason</th><th>Preg</th><th>Live</th><th>Loss</th><th>Pregnancy History</th></tr></thead>'
)

# === 5. Update renderDetailTable to include genotyped column ===
old_render_detail = """function renderDetailTable(key) {
    const data = DETAILS[key];
    if (!data) return;
    const tbody = document.querySelector('#dtable-' + key + ' tbody');
    if (tbody.children.length > 0) return; // already rendered
    const frag = document.createDocumentFragment();
    data.forEach(s => {
        const tr = document.createElement('tr');
        tr.dataset.searchable = (s.id + ' ' + s.events + ' ' + (s.sub || '')).toLowerCase();
        if (key === 'gray') {
            tr.innerHTML = `<td><strong>${s.id}</strong></td><td>${s.age ?? '?'}</td><td>${s.n_preg}</td><td>${s.live}</td><td>${s.loss}</td><td>${colorizeEvents(s.events)}</td>`;
        } else if (key === 'review') {
            tr.innerHTML = `<td><strong>${s.id}</strong></td><td>${s.age ?? '?'}</td><td>${s.n_preg}</td><td>${s.live}</td><td>${s.loss}</td><td style="color:#c62828;font-weight:600;">${s.loss_ge20}</td><td>${colorizeEvents(s.events)}</td>`;
        } else {
            tr.innerHTML = `<td><strong>${s.id}</strong></td><td>${s.age ?? '?'}</td><td><span class="sub-tag">${s.sub || s.reason}</span></td><td>${s.n_preg}</td><td>${s.live}</td><td>${s.loss}</td><td>${colorizeEvents(s.events)}</td>`;
        }
        frag.appendChild(tr);
    });
    tbody.appendChild(frag);
}"""

new_render_detail = """function renderDetailTable(key) {
    const data = DETAILS[key];
    if (!data) return;
    const tbody = document.querySelector('#dtable-' + key + ' tbody');
    if (tbody.children.length > 0) return; // already rendered
    const frag = document.createDocumentFragment();
    const genoCell = s => s.genotyped ? '<td class="geno-yes" title="In genotype array">&#10003;</td>' : '<td class="geno-no" title="Not genotyped">&#10007;</td>';
    data.forEach(s => {
        const tr = document.createElement('tr');
        tr.dataset.searchable = (s.id + ' ' + s.events + ' ' + (s.sub || '') + (s.genotyped ? ' genotyped' : ' not-genotyped')).toLowerCase();
        if (key === 'gray') {
            tr.innerHTML = `<td><strong>${s.id}</strong></td>${genoCell(s)}<td>${s.age ?? '?'}</td><td>${s.n_preg}</td><td>${s.live}</td><td>${s.loss}</td><td>${colorizeEvents(s.events)}</td>`;
        } else if (key === 'review') {
            tr.innerHTML = `<td><strong>${s.id}</strong></td>${genoCell(s)}<td>${s.age ?? '?'}</td><td>${s.n_preg}</td><td>${s.live}</td><td>${s.loss}</td><td style="color:#c62828;font-weight:600;">${s.loss_ge20}</td><td>${colorizeEvents(s.events)}</td>`;
        } else {
            tr.innerHTML = `<td><strong>${s.id}</strong></td>${genoCell(s)}<td>${s.age ?? '?'}</td><td><span class="sub-tag">${s.sub || s.reason}</span></td><td>${s.n_preg}</td><td>${s.live}</td><td>${s.loss}</td><td>${colorizeEvents(s.events)}</td>`;
        }
        frag.appendChild(tr);
    });
    tbody.appendChild(frag);
}"""
html = html.replace(old_render_detail, new_render_detail)

# === 6. Update GWAS-ready info-box ===
html = html.replace(
    '<strong>GWAS-ready summary:</strong> With strict criteria, we have <strong>398 cases &times; 861 controls</strong>\n'
    '    (1,259 samples), excluding 556 subjects (313 gray + 116 review + 94 excluded + 18 unclassified + 15 probable).\n'
    '    After GA review, cases could increase to ~450&ndash;500.',
    
    '<strong>GWAS-ready summary:</strong> With strict criteria, we have <strong>398 cases &times; 861 controls</strong>\n'
    '    (1,259 samples), excluding 556 subjects (313 gray + 116 review + 94 excluded + 18 unclassified + 15 probable).\n'
    '    After GA review, cases could increase to ~450&ndash;500.\n'
    '    <br><br>\n'
    '    <strong>Genotyped subset:</strong> <span class="geno-badge geno-badge-green">270 cases &times; 499 controls = 769 GWAS-ready samples</span>\n'
    '    <br><span style="font-size:0.9em;color:#666;">Of 1,047 genotyped ALSU samples: 270 cases (68%), 499 controls (58%), 119 gray, 62 review, 40 excluded, remainder unclassified/probable.</span>'
)

# === 7. Replace DETAILS blob with updated version (including genotyped field) ===
details_compact = json.dumps(details, ensure_ascii=False, separators=(',', ':'))
old_details_pattern = re.compile(r'const DETAILS = \{.*?\};\n', re.DOTALL)
match = old_details_pattern.search(html)
if match:
    html = html[:match.start()] + 'const DETAILS = ' + details_compact + ';\n' + html[match.end():]
    print(f"Replaced DETAILS blob: {len(details_compact)} chars")
else:
    print("WARNING: Could not find DETAILS blob to replace!")

# === 8. Add geno_stats to inline DATA ===
# Insert geno_stats into the DATA object
geno_insert = json.dumps(geno, separators=(',', ':'))
html = html.replace(
    'const DATA = {',
    'const DATA = {"geno_stats":' + geno_insert + ','
)

with open(HTML_PATH, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"Updated {HTML_PATH}: {len(html)} chars, {html.count(chr(10))+1} lines")

# Verify changes
checks = [
    ('geno_stats in DATA', 'geno_stats' in html),
    ('genotyped in DETAILS', '"genotyped":true' in html or '"genotyped":false' in html),
    ('geno-badge CSS', '.geno-badge' in html),
    ('Geno column header', '<th>Geno</th>' in html),
    ('genoCell function', 'genoCell' in html),
    ('769 GWAS-ready', '769' in html),
    ('geno-badge in gray panel', 'geno-badge-amber">119' in html),
]
for label, ok in checks:
    print(f"  {'✓' if ok else '✗'} {label}")
