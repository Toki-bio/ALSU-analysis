"""Inject sample detail data + JS functions into pre_gwas.html."""
import json

# Read the compact detail data
with open('data/sample_details.json', 'r', encoding='utf-8') as f:
    details = json.load(f)

detail_json = json.dumps(details, ensure_ascii=False, separators=(',',':'))

js_block = f"""
// ===== SAMPLE DETAIL DATA =====
const DETAILS = {detail_json};

// ===== Expandable panels =====
function togglePanel(key) {{
    const panel = document.getElementById('panel-' + key);
    const arrow = document.getElementById('arrow-' + key);
    const isOpen = panel.classList.contains('open');
    // Close all first
    ['gray','review','excluded'].forEach(k => {{
        document.getElementById('panel-' + k).classList.remove('open');
        document.getElementById('arrow-' + k).classList.remove('open');
    }});
    if (!isOpen) {{
        panel.classList.add('open');
        arrow.classList.add('open');
        renderDetailTable(key);
    }}
}}

function colorizeEvents(evtStr) {{
    return evtStr.replace(/(Live[^→]*)/g, '<span class="evt-live">$1</span>')
        .replace(/(Stillbirth[^→]*|Miscarriage[^→]*|Non-dev[^→]*)/g, '<span class="evt-loss">$1</span>')
        .replace(/(Abortion[^→]*|Med-abort[^→]*)/g, '<span class="evt-abort">$1</span>')
        .replace(/(Ectopic[^→]*|Pregnant[^→]*)/g, '<span class="evt-other">$1</span>');
}}

function renderDetailTable(key) {{
    const data = DETAILS[key];
    if (!data) return;
    const tbody = document.querySelector('#dtable-' + key + ' tbody');
    if (tbody.children.length > 0) return; // already rendered
    const frag = document.createDocumentFragment();
    data.forEach(s => {{
        const tr = document.createElement('tr');
        tr.dataset.searchable = (s.id + ' ' + s.events + ' ' + (s.sub || '')).toLowerCase();
        if (key === 'gray') {{
            tr.innerHTML = `<td><strong>${{s.id}}</strong></td><td>${{s.age ?? '?'}}</td><td>${{s.n_preg}}</td><td>${{s.live}}</td><td>${{s.loss}}</td><td>${{colorizeEvents(s.events)}}</td>`;
        }} else if (key === 'review') {{
            tr.innerHTML = `<td><strong>${{s.id}}</strong></td><td>${{s.age ?? '?'}}</td><td>${{s.n_preg}}</td><td>${{s.live}}</td><td>${{s.loss}}</td><td style="color:#c62828;font-weight:600;">${{s.loss_ge20}}</td><td>${{colorizeEvents(s.events)}}</td>`;
        }} else {{
            tr.innerHTML = `<td><strong>${{s.id}}</strong></td><td>${{s.age ?? '?'}}</td><td><span class="sub-tag">${{s.sub || s.reason}}</span></td><td>${{s.n_preg}}</td><td>${{s.live}}</td><td>${{s.loss}}</td><td>${{colorizeEvents(s.events)}}</td>`;
        }}
        frag.appendChild(tr);
    }});
    tbody.appendChild(frag);
}}

function filterDetail(key, query) {{
    const q = query.toLowerCase().trim();
    const rows = document.querySelectorAll('#dtable-' + key + ' tbody tr');
    rows.forEach(tr => {{
        tr.style.display = (!q || tr.dataset.searchable.includes(q)) ? '' : 'none';
    }});
}}
"""

# Read existing HTML
with open('pre_gwas.html', 'r', encoding='utf-8') as f:
    html = f.read()

# Inject before </script>
html = html.replace('</script>\n</body>', js_block + '</script>\n</body>')

with open('pre_gwas.html', 'w', encoding='utf-8') as f:
    f.write(html)

print(f"Injected {len(detail_json):,} chars of detail data + JS functions")
print(f"Final HTML size: {len(html):,} chars")
