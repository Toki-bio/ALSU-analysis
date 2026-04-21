"""Build interactive chr18 locus plot HTML from downloaded data."""
from __future__ import annotations
import json, math
from pathlib import Path

DATA = Path('data/chr18_locus')
OUT = Path('steps/chr18_locus.html')

# --- Load Firth-hybrid top SNPs ---
firth_rows = []
with open(DATA / 'snps_top.tsv') as fh:
    header = fh.readline().lstrip('#').split()
    for line in fh:
        parts = line.rstrip('\n').split('\t')
        r = dict(zip(header, parts))
        firth_rows.append(r)

# --- Load SAIGE results (region) ---
saige_by_id = {}
with open(DATA / 'snps_saige.tsv') as fh:
    sh = fh.readline().split()
    for line in fh:
        parts = line.rstrip('\n').split('\t')
        r = dict(zip(sh, parts))
        rid = r.get('MarkerID')
        if rid:
            saige_by_id[rid] = r

# --- Load LD matrix ---
ld_mat = []
with open(DATA / 'ld_r2_ours.phased.vcor2') as fh:
    for line in fh:
        ld_mat.append([float(x) for x in line.split()])
ld_vars = open(DATA / 'ld_r2_ours.phased.vcor2.vars').read().strip().split('\n')

# --- Genes ---
genes = json.load(open(DATA / 'genes.json'))

# --- Build SNP table ---
def f(x, default='—'):
    try:
        return float(x)
    except Exception:
        return default

snps = []
for r in firth_rows:
    rid = r['ID']
    pos = int(r['POS'])
    p = float(r['P'])
    saige = saige_by_id.get(rid, {})
    snps.append({
        'id': rid,
        'pos': pos,
        'ref': r['REF'],
        'alt': r['ALT'],
        'a1': r['A1'],
        'af': float(r['A1_FREQ']),
        'n': int(r['OBS_CT']),
        'or': float(r['OR']),
        'l95': float(r['L95']),
        'u95': float(r['U95']),
        'se': float(r['LOG(OR)_SE']),
        'z': float(r['Z_STAT']),
        'p': p,
        'logp': -math.log10(p),
        'saige_p': f(saige.get('p.value')),
        'saige_beta': f(saige.get('BETA')),
        'saige_af_case': f(saige.get('AF_case')),
        'saige_af_ctrl': f(saige.get('AF_ctrl')),
    })

# identify lead SNP (min P; tie -> smallest pos)
snps_sorted = sorted(snps, key=lambda s: (s['p'], s['pos']))
lead = snps_sorted[0]

# LD to lead (r²) if lead is in ld_vars, else max r² column
if lead['id'] in ld_vars:
    lead_idx = ld_vars.index(lead['id'])
    r2_by_id = {v: ld_mat[lead_idx][i] for i, v in enumerate(ld_vars)}
else:
    r2_by_id = {}
for s in snps:
    s['r2_lead'] = r2_by_id.get(s['id'])

# --- Prepare plot data ---
pos_min = min(s['pos'] for s in snps) - 15000
pos_max = max(s['pos'] for s in snps) + 15000

plot_data = {
    'snps': snps,
    'lead_id': lead['id'],
    'genes': genes,
    'ld_vars': ld_vars,
    'ld_mat': ld_mat,
    'window': [pos_min, pos_max],
}

TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>chr18:70.47 Mb locus — ALSU RPL GWAS</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root {
  --bg:#0e1117; --panel:#161b22; --fg:#c9d1d9; --muted:#8b949e;
  --accent:#58a6ff; --warn:#f0883e; --good:#3fb950; --bad:#f85149;
}
* { box-sizing: border-box; }
body {
  background: var(--bg); color: var(--fg);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  margin: 0; padding: 20px; line-height: 1.5;
}
h1 { margin:0 0 4px 0; font-size: 22px; color: var(--accent); }
h2 { margin-top: 28px; border-bottom: 1px solid #30363d; padding-bottom: 4px; }
.sub { color: var(--muted); font-size: 13px; }
.card { background: var(--panel); border: 1px solid #30363d; border-radius: 8px; padding: 14px 18px; margin: 12px 0; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; }
.stat { background:#0d1117; border:1px solid #30363d; border-radius:6px; padding:8px 12px; }
.stat .k { color: var(--muted); font-size: 11px; text-transform: uppercase; }
.stat .v { font-size: 18px; font-weight: 600; color: var(--accent); }
table { width: 100%; border-collapse: collapse; font-size: 12px; font-family: 'Consolas','Menlo',monospace; }
th, td { padding: 5px 8px; border-bottom: 1px solid #30363d; text-align: right; }
th:first-child, td:first-child { text-align: left; }
th { background:#21262d; cursor: pointer; user-select: none; position: sticky; top:0; }
th:hover { background:#30363d; }
tr:hover td { background:#1c2128; }
tr.lead td { background:#1f2937 !important; font-weight: 600; }
tr.selected td { background:#264653 !important; }
.r2-high { color: #ff7b72; }
.r2-mid { color: #ffa657; }
.r2-low { color: #8b949e; }
a { color: var(--accent); }
#figure { width: 100%; min-height: 800px; }
.note { color: var(--muted); font-size: 12px; margin-top: 4px; }
.tag { display:inline-block; padding: 2px 6px; border-radius: 4px; font-size: 11px; background:#21262d; color:#8b949e; margin-right: 4px;}
</style>
</head>
<body>

<h1>chr18 locus — RPL GWAS lead signal</h1>
<div class="sub">
  Region chr18:70,469,147–70,495,261 (GRCh38) · 47 SNPs at p&lt;1×10⁻⁵ · Firth-hybrid (PLINK2), cross-referenced with SAIGE
</div>

<div class="card">
  <div class="grid">
    <div class="stat"><div class="k">Lead SNP</div><div class="v" id="lead_id"></div></div>
    <div class="stat"><div class="k">Lead P (Firth)</div><div class="v" id="lead_p"></div></div>
    <div class="stat"><div class="k">Lead OR (95% CI)</div><div class="v" id="lead_or"></div></div>
    <div class="stat"><div class="k">Lead MAF (cohort)</div><div class="v" id="lead_af"></div></div>
    <div class="stat"><div class="k">Variants in peak</div><div class="v" id="n_snps"></div></div>
    <div class="stat"><div class="k">Span</div><div class="v" id="span"></div></div>
  </div>
  <div class="note" style="margin-top:10px;">
    <b>Genomic context:</b> intergenic region between <b>SOCS6</b> (~140 kb upstream, suppressor of cytokine signaling 6) and distal <b>NETO1</b> (~2.2 Mb downstream).
    Lead variants overlap lncRNA <b>ENSG00000288828</b> (70,464,905–70,592,638). All lead alleles are <b>protective</b> (OR &lt; 0.5) and <b>low-frequency</b> (MAF 11–14 %).
    The 47 peak SNPs form a single tight LD block (minimum pairwise r² ≈ 0.67 in our cohort).
  </div>
</div>

<h2>Interactive figure</h2>
<div class="card">
  <div id="figure"></div>
  <div class="note">Top: association panel (–log₁₀ P). Middle: gene/feature track. Bottom: pairwise r² (our cohort, 1090 post-QC samples). Hover for full annotation. Click a variant to highlight in the table below.</div>
</div>

<h2>All 47 peak variants</h2>
<div class="card" style="max-height:560px; overflow:auto;">
<table id="snptable">
<thead>
  <tr>
    <th data-k="id">rsID</th>
    <th data-k="pos">Position (GRCh38)</th>
    <th data-k="ref">REF</th>
    <th data-k="alt">ALT</th>
    <th data-k="af">A1 AF</th>
    <th data-k="or">OR</th>
    <th data-k="l95">L95</th>
    <th data-k="u95">U95</th>
    <th data-k="p">Firth P</th>
    <th data-k="saige_p">SAIGE P</th>
    <th data-k="saige_af_case">AF cases (SAIGE)</th>
    <th data-k="saige_af_ctrl">AF ctrls (SAIGE)</th>
    <th data-k="r2_lead">r² (lead)</th>
  </tr>
</thead>
<tbody></tbody>
</table>
</div>

<div class="card note">
<b>Columns:</b> Firth P from PLINK2 <code>--glm firth-fallback</code> (n=764 at this locus).
SAIGE P from mixed-model Firth-corrected SPA (n=764, saddlepoint approximation).
A1 = effect allele (minor); OR relative to A1.
AF cases / AF controls from SAIGE (case n=256, ctrl n=508).
r² = pairwise LD to lead SNP in our cohort (PLINK2 <code>--r2-phased square</code> on imputed HQ data, 1090 samples).
</div>

<h2>Interpretation notes</h2>
<div class="card">
<ul>
  <li><b>Single causal signal.</b> All 47 peak variants in tight LD (min pairwise r² ≈ 0.67) over ~26 kb — single underlying haplotype, fine-mapping within cohort impossible without external reference.</li>
  <li><b>Allele direction.</b> All leads protective (OR 0.34–0.43) — minor allele reduces RPL risk.</li>
  <li><b>Nearest protein-coding gene: SOCS6</b> (~140 kb upstream). SOCS6 is a negative regulator of JAK/STAT signaling, with roles in cytokine response — biologically plausible for pregnancy/immune phenotype. Not yet reported in RPL GWAS.</li>
  <li><b>Novel non-coding context.</b> Peak overlaps a lincRNA transcript (ENSG00000288828 / novel transcript) — possible regulatory role for SOCS6 or distal NETO1.</li>
  <li><b>Not genome-wide significant.</b> Firth p ≈ 3×10⁻⁷; below 5×10⁻⁸ GWS threshold. Requires replication in an independent RPL cohort.</li>
  <li><b>SAIGE confirms.</b> Best SAIGE p ≈ 4×10⁻⁶ at the same lead position — concordant direction; slightly attenuated by mixed-model shrinkage.</li>
</ul>
</div>

<script id="payload" type="application/json">__PAYLOAD__</script>
<script>
const PAYLOAD = JSON.parse(document.getElementById('payload').textContent);
const snps = PAYLOAD.snps;
const leadId = PAYLOAD.lead_id;
const genes = PAYLOAD.genes;
const ldVars = PAYLOAD.ld_vars;
const ldMat = PAYLOAD.ld_mat;
const [x0, x1] = PAYLOAD.window;
const lead = snps.find(s => s.id === leadId);

// --- Top-of-page stats ---
const fmt = (x, d=3) => (x==null || x==='—') ? '—' : (typeof x==='number' ? x.toExponential(d) : x);
const fmtP = (x) => (x==null || x==='—') ? '—' : Number(x).toExponential(2);
const fmtF = (x, d=3) => (x==null || x==='—') ? '—' : Number(x).toFixed(d);

document.getElementById('lead_id').textContent = lead.id;
document.getElementById('lead_p').textContent = fmtP(lead.p);
document.getElementById('lead_or').textContent = `${fmtF(lead.or,2)} (${fmtF(lead.l95,2)}–${fmtF(lead.u95,2)})`;
document.getElementById('lead_af').textContent = fmtF(lead.af,3);
document.getElementById('n_snps').textContent = snps.length;
document.getElementById('span').textContent = `${(Math.min(...snps.map(s=>s.pos))/1e6).toFixed(3)}–${(Math.max(...snps.map(s=>s.pos))/1e6).toFixed(3)} Mb`;

// --- Color by LD to lead ---
function r2Color(r2) {
  if (r2 == null) return '#8b949e';
  if (r2 >= 0.8) return '#f85149';
  if (r2 >= 0.6) return '#f0883e';
  if (r2 >= 0.4) return '#d2a8ff';
  if (r2 >= 0.2) return '#79c0ff';
  return '#6e7681';
}

// --- Build traces ---
// Panel 1 (y domain 0.65-1.0): -log10 P — one trace per method
const hoverFirth = snps.map(s => 
  `<b>${s.id}</b><br>chr18:${s.pos.toLocaleString()}<br>` +
  `${s.ref}>${s.alt} (A1=${s.a1})<br>` +
  `AF=${fmtF(s.af,3)} | N=${s.n}<br>` +
  `OR=${fmtF(s.or,3)} [${fmtF(s.l95,3)}, ${fmtF(s.u95,3)}]<br>` +
  `SE(logOR)=${fmtF(s.se,3)} | Z=${fmtF(s.z,2)}<br>` +
  `<b>Firth P = ${fmtP(s.p)}</b><br>` +
  `SAIGE P = ${fmtP(s.saige_p)}<br>` +
  `r² to lead = ${s.r2_lead==null ? '—' : fmtF(s.r2_lead,3)}`
);
const firthTrace = {
  x: snps.map(s => s.pos), y: snps.map(s => s.logp),
  mode: 'markers', type: 'scatter', name: 'Firth',
  marker: {
    size: snps.map(s => s.id === leadId ? 16 : 10),
    color: snps.map(s => r2Color(s.r2_lead)),
    symbol: snps.map(s => s.id === leadId ? 'diamond' : 'circle'),
    line: { width: snps.map(s => s.id === leadId ? 2 : 1), color: '#0e1117' }
  },
  text: hoverFirth, hovertemplate: '%{text}<extra></extra>',
  xaxis: 'x', yaxis: 'y1', customdata: snps.map(s => s.id),
};
const saigeSnps = snps.filter(s => typeof s.saige_p === 'number');
const saigeTrace = {
  x: saigeSnps.map(s => s.pos),
  y: saigeSnps.map(s => -Math.log10(s.saige_p)),
  mode: 'markers', type: 'scatter', name: 'SAIGE',
  marker: { size: 6, color: '#8b949e', symbol: 'triangle-down', opacity: 0.7 },
  text: saigeSnps.map(s => `${s.id}<br>SAIGE P = ${fmtP(s.saige_p)}`),
  hovertemplate: '%{text}<extra></extra>',
  xaxis: 'x', yaxis: 'y1',
};
// Significance lines as shapes (GWS and suggestive)
const shapes = [
  { type:'line', xref:'x', yref:'y1', x0:x0, x1:x1, y0:-Math.log10(5e-8), y1:-Math.log10(5e-8),
    line:{ color:'#f85149', width:1, dash:'dash' } },
  { type:'line', xref:'x', yref:'y1', x0:x0, x1:x1, y0:-Math.log10(1e-5), y1:-Math.log10(1e-5),
    line:{ color:'#f0883e', width:1, dash:'dot' } },
];

// Panel 2: genes — horizontal bars
const geneTraces = [];
const laneByGene = {};
// Assign lanes greedy
const sortedGenes = [...genes].sort((a,b)=>a.start-b.start);
const lanes = [];
sortedGenes.forEach(g => {
  let placed = false;
  for (let i=0;i<lanes.length;i++) {
    if (lanes[i] < g.start - 2000) { lanes[i] = g.end; laneByGene[g.id] = i; placed = true; break; }
  }
  if (!placed) { laneByGene[g.id] = lanes.length; lanes.push(g.end); }
});
const nLanes = Math.max(1, lanes.length);

genes.forEach(g => {
  const lane = laneByGene[g.id];
  const name = g.external_name || g.id;
  const btype = g.biotype;
  const color = btype === 'protein_coding' ? '#3fb950'
              : btype === 'lncRNA' ? '#79c0ff'
              : btype === 'processed_pseudogene' ? '#6e7681'
              : '#d2a8ff';
  geneTraces.push({
    x: [g.start, g.end], y: [lane, lane], mode:'lines',
    line: { color: color, width: 8 },
    hoverinfo: 'text',
    text: `<b>${name}</b><br>${btype}<br>${g.start.toLocaleString()}–${g.end.toLocaleString()}<br>strand ${g.strand>0?'+':'−'}`,
    showlegend: false, xaxis:'x', yaxis:'y2',
  });
});
const geneAnnots = genes.map(g => ({
  x: (Math.max(g.start,x0)+Math.min(g.end,x1))/2,
  y: laneByGene[g.id] + 0.35, xref:'x', yref:'y2',
  text: (g.external_name || g.id).replace('ENSG00000','ENSG..'),
  showarrow: false, font: { size: 10, color:'#c9d1d9' },
}));

// Panel 3: LD triangle heatmap
// Sort variants by position
const varPos = {};
snps.forEach(s => { varPos[s.id] = s.pos; });
const ldOrder = ldVars.map((v,i) => ({v, i, pos: varPos[v] || 0})).sort((a,b)=>a.pos-b.pos);
const ldLabels = ldOrder.map(o => o.v);
const ldPositions = ldOrder.map(o => o.pos);
// Reorder matrix — only show upper triangle
const n = ldOrder.length;
const Z = Array.from({length:n}, () => Array(n).fill(null));
for (let i=0;i<n;i++){
  for (let j=0;j<n;j++){
    const r2 = ldMat[ldOrder[i].i][ldOrder[j].i];
    if (j >= i) Z[i][j] = r2;
  }
}
const heat = {
  z: Z, x: ldLabels, y: ldLabels,
  type:'heatmap',
  colorscale: [[0,'#0d1117'],[0.2,'#0d47a1'],[0.4,'#1976d2'],[0.6,'#ffa726'],[0.8,'#ef5350'],[1,'#b71c1c']],
  zmin: 0, zmax: 1,
  colorbar: { title:'r²', len:0.25, y:0.12, x:1.02, thickness:10, tickfont:{color:'#c9d1d9'} },
  hovertemplate: '%{y} ↔ %{x}<br>r² = %{z:.3f}<extra></extra>',
  xaxis:'x2', yaxis:'y3',
};

const layout = {
  paper_bgcolor:'#161b22', plot_bgcolor:'#0d1117',
  font: { color:'#c9d1d9', size: 12 },
  margin: { t: 30, l: 70, r: 80, b: 60 },
  showlegend: true,
  legend: { orientation:'h', x:0, y:1.06, bgcolor: 'rgba(0,0,0,0)' },
  grid: { rows: 3, columns: 1, pattern:'independent' },
  xaxis:  { domain:[0,1], range:[x0,x1], showgrid:false, zeroline:false,
            anchor:'y2', title:'', tickformat:',' },
  xaxis2: { domain:[0,1], showgrid:false, zeroline:false, anchor:'y3',
            tickangle:-60, tickfont:{size:9} },
  yaxis:  { domain:[0.62, 1.0], title:'−log₁₀ P', zeroline:false, gridcolor:'#21262d' },
  yaxis2: { domain:[0.50, 0.60], range:[-0.5, nLanes-0.2], showticklabels:false, zeroline:false, showgrid:false, title:'genes' },
  yaxis3: { domain:[0.00, 0.48], autorange:'reversed', tickfont:{size:9}, scaleanchor:'x2', scaleratio:1 },
  shapes: shapes,
  annotations: [
    ...geneAnnots,
    { x: x1, y: -Math.log10(5e-8), xref:'x', yref:'y1', text: 'GWS 5×10⁻⁸', showarrow:false, font:{ size:10, color:'#f85149'}, xanchor:'right', yanchor:'bottom' },
    { x: x1, y: -Math.log10(1e-5), xref:'x', yref:'y1', text: 'suggestive 10⁻⁵', showarrow:false, font:{ size:10, color:'#f0883e'}, xanchor:'right', yanchor:'bottom' },
  ],
  height: 820,
};

Plotly.newPlot('figure', [firthTrace, saigeTrace, ...geneTraces, heat], layout, { responsive: true, displaylogo: false });

// --- Table ---
const tbody = document.querySelector('#snptable tbody');
let sortKey = 'p', sortDir = 1;
function fmtCell(v, d=3) {
  if (v==null || v==='—') return '—';
  if (typeof v !== 'number') return v;
  return Math.abs(v) < 1e-3 && v !== 0 ? v.toExponential(2) : v.toFixed(d);
}
function r2cls(r2) { if (r2==null) return ''; if (r2>=0.8) return 'r2-high'; if (r2>=0.5) return 'r2-mid'; return 'r2-low'; }
function renderTable() {
  const sorted = [...snps].sort((a,b)=>{
    const av = a[sortKey], bv = b[sortKey];
    if (av == null || av === '—') return 1;
    if (bv == null || bv === '—') return -1;
    if (av < bv) return -1*sortDir;
    if (av > bv) return 1*sortDir;
    return 0;
  });
  tbody.innerHTML = sorted.map(s => {
    const cls = (s.id === leadId ? 'lead' : '');
    return `<tr class="${cls}" data-id="${s.id}">
      <td><a href="https://www.ensembl.org/Homo_sapiens/Variation/Explore?v=${s.id}" target="_blank">${s.id}</a></td>
      <td>${s.pos.toLocaleString()}</td>
      <td>${s.ref}</td><td>${s.alt}</td>
      <td>${fmtCell(s.af,3)}</td>
      <td>${fmtCell(s.or,3)}</td>
      <td>${fmtCell(s.l95,3)}</td>
      <td>${fmtCell(s.u95,3)}</td>
      <td>${fmtCell(s.p,3)}</td>
      <td>${fmtCell(s.saige_p,3)}</td>
      <td>${fmtCell(s.saige_af_case,3)}</td>
      <td>${fmtCell(s.saige_af_ctrl,3)}</td>
      <td class="${r2cls(s.r2_lead)}">${fmtCell(s.r2_lead,3)}</td>
    </tr>`;
  }).join('');
}
renderTable();

document.querySelectorAll('#snptable th').forEach(th => {
  th.addEventListener('click', () => {
    const k = th.dataset.k;
    if (sortKey === k) sortDir *= -1; else { sortKey = k; sortDir = 1; }
    renderTable();
  });
});

// Click in plot -> highlight row
document.getElementById('figure').on('plotly_click', (ev) => {
  const p = ev.points && ev.points[0];
  if (!p || !p.customdata) return;
  const id = p.customdata;
  document.querySelectorAll('#snptable tr').forEach(tr => tr.classList.remove('selected'));
  const tr = document.querySelector(`#snptable tr[data-id="${id}"]`);
  if (tr) { tr.classList.add('selected'); tr.scrollIntoView({block:'center', behavior:'smooth'}); }
});

</script>
</body>
</html>
"""

html = TEMPLATE.replace('__PAYLOAD__', json.dumps(plot_data, separators=(',', ':')))
OUT.parent.mkdir(exist_ok=True, parents=True)
OUT.write_text(html, encoding='utf-8')
print(f'Wrote {OUT} ({OUT.stat().st_size} bytes)')
print(f'Lead SNP: {lead["id"]} at chr18:{lead["pos"]:,}  p={lead["p"]:.2e}  OR={lead["or"]:.3f} [{lead["l95"]:.3f}-{lead["u95"]:.3f}]  AF={lead["af"]:.3f}')
print(f'Number of SNPs in plot: {len(snps)}')
print(f'Genes in window: {len(genes)}')
print(f'LD matrix: {len(ld_mat)}x{len(ld_mat[0])}')
