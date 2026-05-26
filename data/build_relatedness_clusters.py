#!/usr/bin/env python3
"""
Build steps/relatedness_clusters.html: D3.js force graph + cluster table.

Thresholds offered:
  T1: PI_HAT >= 0.354  (1st degree, PI_HAT-based)
  T2: PI_HAT >= 0.177  (2nd degree, PI_HAT-based; inflated by admixture)
  T3: KING >= 0.177    (1st degree, KING-only)
  T4: KING >= 0.0884   (2nd degree, KING-only)
"""
import os, json
from collections import defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

NOTABLE = {
    "ALSU_0832": "ALSU_0832 cluster (PCRel x4)",
    "ALSU_0407": "PCRel confirmed (0832 cluster)",
    "ALSU_0808": "PCRel confirmed (0832 cluster)",
    "ALSU_0787": "PCRel confirmed (0832 cluster)",
    "ALSU_0855": "PCRel confirmed (0832 cluster)",
    "ALSU_0218": "PCRel+RelAdmix confirmed",
    "ALSU_0242": "PCRel+RelAdmix confirmed",
}

PC_RELATE_PAIRS = {
    tuple(sorted(["ALSU_0407", "ALSU_0832"])): 0.7955,
    tuple(sorted(["ALSU_0808", "ALSU_0832"])): 0.5753,
    tuple(sorted(["ALSU_0787", "ALSU_0832"])): 0.5435,
    tuple(sorted(["ALSU_0832", "ALSU_0855"])): 0.1570,
    tuple(sorted(["ALSU_0218", "ALSU_0242"])): 0.1532,
}

# ── 1. Sample map ──────────────────────────────────────────────────────────────
sample_map = {}
with open(os.path.join(BASE, 'dragen_qc_data/cross_gwas2026_alsu_winter/sample_map.tsv'), encoding='utf-8-sig') as f:
    hdr = f.readline().strip().split('\t')
    ci = hdr.index('cohort'); si = hdr.index('safe_iid'); oi = hdr.index('original_iid')
    for line in f:
        p = line.strip().split('\t')
        if len(p) > max(ci, si, oi):
            sample_map[p[si]] = {'cohort': p[ci], 'original': p[oi]}

# ── 2. PI_HAT pairs ────────────────────────────────────────────────────────────
all_pihat = {}
with open(os.path.join(BASE, 'dragen_qc_data/cross_gwas2026_alsu_winter/cross_ibd_pihat_0125.genome'), encoding='utf-8-sig') as f:
    hdr = f.readline().split()
    pi_i = hdr.index('PI_HAT'); i1_i = hdr.index('IID1'); i2_i = hdr.index('IID2')
    for line in f:
        p = line.split()
        if len(p) <= pi_i: continue
        try: pi = float(p[pi_i])
        except: continue
        all_pihat[tuple(sorted([p[i1_i], p[i2_i]]))] = pi

# ── 3. KING pairs ──────────────────────────────────────────────────────────────
all_king = {}
with open(os.path.join(BASE, 'dragen_qc_data/cross_gwas2026_alsu_winter/cross_king_00884.kin0'), encoding='utf-8') as f:
    hdr = f.readline().split(); hdr[0] = hdr[0].lstrip('#')
    ki_i = hdr.index('KINSHIP'); i1_i = hdr.index('IID1'); i2_i = hdr.index('IID2')
    for line in f:
        p = line.split()
        if len(p) <= ki_i: continue
        try: k = float(p[ki_i])
        except: continue
        all_king[tuple(sorted([p[i1_i], p[i2_i]]))] = k

# ── 4. RelateAdmix top-100 ────────────────────────────────────────────────────
all_ra = {}
with open(os.path.join(BASE, 'data/relateadmix_top_pairs.tsv'), encoding='utf-8') as f:
    hdr = f.readline().strip().split('\t')
    i1_i = hdr.index('iid1'); i2_i = hdr.index('iid2'); kin_i = hdr.index('relateadmix_kin_equiv')
    for line in f:
        p = line.strip().split('\t')
        if len(p) <= kin_i: continue
        try: kin = float(p[kin_i])
        except: continue
        all_ra[tuple(sorted([p[i1_i], p[i2_i]]))] = kin

print(f"Loaded: {len(all_pihat)} PI_HAT | {len(all_king)} KING | {len(all_ra)} RelAdmix | {len(sample_map)} samples")

# ── Connected components ───────────────────────────────────────────────────────
def connected_components(pairs):
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[rb] = ra
    for a, b in pairs: union(a, b)
    groups = defaultdict(set)
    for n in parent: groups[find(n)].add(n)
    return [sorted(g) for g in groups.values() if len(g) >= 2]

# ── Build dataset ─────────────────────────────────────────────────────────────
def build_dataset(pihat_min, king_min, label, warn=''):
    edge_map = {}
    if pihat_min < 1.0:
        for key, pi in all_pihat.items():
            if pi >= pihat_min:
                edge_map[key] = {'s': key[0], 't': key[1], 'pi': pi, 'ki': None, 'ra': None, 'pc': None, 'm': 'pi'}
    for key, k in all_king.items():
        if k >= king_min:
            if key in edge_map:
                edge_map[key]['ki'] = k; edge_map[key]['m'] = 'both'
            else:
                edge_map[key] = {'s': key[0], 't': key[1], 'pi': None, 'ki': k, 'ra': None, 'pc': None, 'm': 'ki'}
    for key in list(edge_map.keys()):
        edge_map[key]['ra'] = all_ra.get(key)
        edge_map[key]['pc'] = PC_RELATE_PAIRS.get(key)

    edges = list(edge_map.values())
    clusters_raw = connected_components([(e['s'], e['t']) for e in edges])

    node_to_cluster = {}
    for ci, cluster in enumerate(sorted(clusters_raw, key=len, reverse=True)):
        for n in cluster: node_to_cluster[n] = ci

    all_nodes_ids = set(e['s'] for e in edges) | set(e['t'] for e in edges)
    nodes = []
    for nid in sorted(all_nodes_ids):
        info = sample_map.get(nid, {'cohort': 'unknown', 'original': nid})
        nodes.append({'id': nid, 'co': info['cohort'], 'or': info['original'],
                      'cl': node_to_cluster.get(nid, -1), 'no': NOTABLE.get(nid, '')})

    cluster_edges = defaultdict(list)
    for e in edges:
        cluster_edges[node_to_cluster.get(e['s'], -1)].append(e)

    cluster_members = defaultdict(list)
    for n in nodes: cluster_members[n['cl']].append(n['id'])

    clusters = []
    for ci in sorted(cluster_members.keys()):
        members = sorted(cluster_members[ci])
        ce = cluster_edges[ci]
        max_pi = max((e['pi'] or 0) for e in ce) if ce else 0
        max_ki = max((e['ki'] or 0) for e in ce) if ce else 0
        notable_m = [n for n in members if n in NOTABLE]
        clusters.append({
            'id': ci,
            'sz': len(members),
            'mb': members,
            'mp': round(max_pi, 4),
            'mk': round(max_ki, 4),
            'ra': any(e['ra'] is not None for e in ce),
            'pc': any(e['pc'] is not None for e in ce),
            'ne': len(ce),
            'no': '; '.join(NOTABLE[n] for n in notable_m[:3]) if notable_m else ''
        })

    # Round edge values
    for e in edges:
        if e['pi'] is not None: e['pi'] = round(e['pi'], 4)
        if e['ki'] is not None: e['ki'] = round(e['ki'], 4)
        if e['ra'] is not None: e['ra'] = round(e['ra'], 4)
        if e['pc'] is not None: e['pc'] = round(e['pc'], 4)

    print(f"  {label}: {len(nodes)} nodes, {len(edges)} edges, {len(clusters)} clusters, "
          f"largest={max((c['sz'] for c in clusters), default=0)}")
    return {
        'label': label, 'warn': warn,
        'nn': len(nodes), 'ne': len(edges), 'nc': len(clusters),
        'nodes': nodes, 'edges': edges, 'clusters': clusters
    }

print("Building datasets...")
datasets = {
    't1': build_dataset(0.354, 999, '1st degree — PI_HAT ≥ 0.354'),
    't2': build_dataset(0.177, 999, '2nd degree — PI_HAT ≥ 0.177',
                        warn='⚠ PI_HAT 2nd-degree pairs are inflated by admixture background. Most pairs in this view are NOT true relatives. Use KING-based thresholds for reliable clusters.'),
    't3': build_dataset(999, 0.177, '1st degree — KING ≥ 0.177'),
    't4': build_dataset(999, 0.0884, '2nd degree — KING ≥ 0.0884'),
}

# ── Build HTML ─────────────────────────────────────────────────────────────────
data_json = json.dumps(datasets, separators=(',', ':'))

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Relatedness Clusters — ALSU+GWAS2026</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:'Segoe UI',Arial,sans-serif;background:#f5f7fa;color:#222}
.page-header{background:linear-gradient(135deg,#1565c0,#0d47a1);color:#fff;padding:18px 28px}
.page-header h1{font-size:1.4rem;font-weight:700}
.page-header p{font-size:.85rem;opacity:.85;margin-top:4px}
.controls{display:flex;gap:10px;flex-wrap:wrap;padding:14px 20px;background:#fff;border-bottom:1px solid #dde}
.btn-thresh{padding:7px 16px;border:2px solid #c5cae9;border-radius:20px;background:#fff;cursor:pointer;font-size:.82rem;font-weight:600;color:#3949ab;transition:.15s}
.btn-thresh.active{background:#3949ab;color:#fff;border-color:#3949ab}
.btn-thresh:hover:not(.active){background:#e8eaf6}
.warn-box{background:#fff3e0;border-left:4px solid #ff9800;padding:8px 14px;font-size:.82rem;color:#6d4c00;margin:8px 20px;border-radius:4px;display:none}
.stats-bar{display:flex;gap:24px;padding:10px 22px;background:#f9f9fb;border-bottom:1px solid #e5e8f0;font-size:.83rem;color:#444}
.stats-bar span{font-weight:700;color:#1565c0}
.main-area{display:flex;gap:0;height:calc(100vh - 200px);min-height:480px}
#graph-area{flex:1;position:relative;background:#fff;border-right:1px solid #dde}
#graph-svg{width:100%;height:100%}
.graph-overlay{position:absolute;top:8px;right:12px;font-size:.75rem;color:#888;background:rgba(255,255,255,.8);padding:4px 8px;border-radius:4px;pointer-events:none}
.panel-msg{position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);text-align:center;color:#777;font-size:.9rem;max-width:340px;line-height:1.5}
.legend{position:absolute;bottom:10px;left:10px;background:rgba(255,255,255,.92);padding:8px 12px;border-radius:6px;font-size:.72rem;box-shadow:0 1px 4px rgba(0,0,0,.1)}
.legend-row{display:flex;align-items:center;gap:6px;margin-bottom:3px}
.legend-dot{width:10px;height:10px;border-radius:50%;flex-shrink:0}
.legend-line{width:22px;height:3px;flex-shrink:0}
.table-panel{width:400px;overflow-y:auto;background:#fff}
.table-panel table{width:100%;border-collapse:collapse;font-size:.78rem}
.table-panel th{position:sticky;top:0;background:#3949ab;color:#fff;padding:7px 8px;text-align:left;cursor:pointer;user-select:none;white-space:nowrap}
.table-panel th:hover{background:#303f9f}
.table-panel td{padding:5px 8px;border-bottom:1px solid #eee;vertical-align:top}
.table-panel tr.highlighted td{background:#e8eaf6}
.table-panel tr:hover td{background:#f3f4fc;cursor:pointer}
.badge-pc{background:#b71c1c;color:#fff;font-size:.65rem;padding:1px 5px;border-radius:8px;white-space:nowrap}
.badge-ra{background:#1b5e20;color:#fff;font-size:.65rem;padding:1px 5px;border-radius:8px;white-space:nowrap}
.badge-no{background:#e65100;color:#fff;font-size:.65rem;padding:1px 5px;border-radius:8px;white-space:nowrap;display:block;margin-top:2px}
.member-list{max-height:60px;overflow-y:auto;line-height:1.4}
#tooltip{position:fixed;background:rgba(20,30,60,.9);color:#fff;padding:8px 12px;border-radius:6px;font-size:.75rem;pointer-events:none;display:none;max-width:280px;z-index:999;line-height:1.5}
</style>
</head>
<body>
<div class="page-header">
  <h1>Relatedness Clusters — ALSU + GWAS2026 (N=1,193)</h1>
  <p>Connected components by PI_HAT / KING thresholds. Edges coloured by method: <b style="color:#a5d6a7">both</b> · <b style="color:#81d4fa">KING only</b> · <b style="color:#ffab91">PI_HAT only</b>. Nodes: <b style="color:#90caf9">ALSU</b> · <b style="color:#ffcc80">GWAS2026</b>.</p>
</div>

<div class="controls">
  <button class="btn-thresh" data-key="t1">1st degree — PI_HAT ≥ 0.354</button>
  <button class="btn-thresh" data-key="t2">2nd degree — PI_HAT ≥ 0.177 ⚠</button>
  <button class="btn-thresh active" data-key="t3">1st degree — KING ≥ 0.177</button>
  <button class="btn-thresh" data-key="t4">2nd degree — KING ≥ 0.0884</button>
</div>
<div class="warn-box" id="warn-box"></div>

<div class="stats-bar" id="stats-bar">
  <div>Nodes: <span id="s-nodes">—</span></div>
  <div>Edges: <span id="s-edges">—</span></div>
  <div>Clusters: <span id="s-clusters">—</span></div>
  <div>Largest: <span id="s-largest">—</span></div>
</div>

<div class="main-area">
  <div id="graph-area">
    <svg id="graph-svg"></svg>
    <div class="graph-overlay" id="graph-overlay"></div>
    <div class="panel-msg" id="panel-msg" style="display:none"></div>
    <div class="legend">
      <div class="legend-row"><div class="legend-dot" style="background:#1565c0"></div>ALSU</div>
      <div class="legend-row"><div class="legend-dot" style="background:#e65100"></div>GWAS2026</div>
      <div class="legend-row"><div class="legend-dot" style="background:gold;border:2px solid #b71c1c"></div>Notable (PCRel/RelAdmix)</div>
      <div style="margin-top:5px;border-top:1px solid #ddd;padding-top:5px">
      <div class="legend-row"><div class="legend-line" style="background:#2e7d32"></div>Both KING+PI_HAT</div>
      <div class="legend-row"><div class="legend-line" style="background:#0277bd"></div>KING only</div>
      <div class="legend-row"><div class="legend-line" style="background:#e53935"></div>PI_HAT only</div>
      </div>
    </div>
  </div>
  <div class="table-panel">
    <table id="cluster-table">
      <thead><tr>
        <th data-col="id" title="Click to sort">#</th>
        <th data-col="sz">Size</th>
        <th data-col="mp">Max PI_HAT</th>
        <th data-col="mk">Max KING</th>
        <th>Members</th>
        <th>Flags</th>
      </tr></thead>
      <tbody id="cluster-tbody"></tbody>
    </table>
  </div>
</div>
<div id="tooltip"></div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
const ALL_DATA = """ + data_json + r""";

// ── State ──────────────────────────────────────────────────────────────────────
let currentKey = 't3';
let sortCol = 'sz', sortDir = -1;
let highlightedCluster = null;
let simulation = null;
const MAX_GRAPH_EDGES = 2000;

// ── Colour helpers ─────────────────────────────────────────────────────────────
const nodeColor  = d => d.co === 'ALSU' ? '#1565c0' : '#e65100';
const edgeColor  = d => d.m === 'both' ? '#2e7d32' : d.m === 'ki' ? '#0277bd' : '#e53935';
const edgeWidth  = d => {
  const v = Math.max(d.pi || 0, (d.ki || 0) * 2);
  return Math.max(1, Math.min(5, v * 6));
};

// ── D3 graph ──────────────────────────────────────────────────────────────────
function renderGraph(ds) {
  const svg = d3.select('#graph-svg');
  svg.selectAll('*').remove();
  if (simulation) { simulation.stop(); simulation = null; }

  const panelMsg = document.getElementById('panel-msg');
  panelMsg.style.display = 'none';

  if (ds.ne > MAX_GRAPH_EDGES) {
    panelMsg.innerHTML = `<b>Graph view disabled</b><br>${ds.ne.toLocaleString()} edges exceed the ${MAX_GRAPH_EDGES.toLocaleString()}-edge render limit.<br><br>Use the table on the right to explore clusters, or switch to a stricter threshold.`;
    panelMsg.style.display = '';
    document.getElementById('graph-overlay').textContent = '';
    return;
  }

  const area = document.getElementById('graph-area');
  const W = area.clientWidth, H = area.clientHeight;

  const nodeMap = {};
  ds.nodes.forEach(n => nodeMap[n.id] = {...n, x: W/2 + (Math.random()-.5)*200, y: H/2 + (Math.random()-.5)*200});
  const nodes = Object.values(nodeMap);
  const edges = ds.edges.map(e => ({...e, source: e.s, target: e.t}));

  // Group nodes by cluster for initial placement
  const clusterCenters = {};
  nodes.forEach(n => {
    if (!(n.cl in clusterCenters)) clusterCenters[n.cl] = {x: W/2+(Math.random()-.5)*W*.6, y: H/2+(Math.random()-.5)*H*.6, count: 0};
    clusterCenters[n.cl].count++;
  });
  nodes.forEach(n => {
    const c = clusterCenters[n.cl];
    n.x = c.x + (Math.random()-.5)*80;
    n.y = c.y + (Math.random()-.5)*80;
  });

  const g = svg.append('g');

  // Zoom
  svg.call(d3.zoom().scaleExtent([.1, 8]).on('zoom', ({transform}) => g.attr('transform', transform)));

  const linkSel = g.append('g').selectAll('line').data(edges).join('line')
    .attr('stroke', edgeColor)
    .attr('stroke-width', edgeWidth)
    .attr('stroke-opacity', .7);

  const nodeSel = g.append('g').selectAll('circle').data(nodes).join('circle')
    .attr('r', d => d.no ? 9 : 6)
    .attr('fill', d => d.no ? 'gold' : nodeColor(d))
    .attr('stroke', d => d.no ? '#b71c1c' : '#fff')
    .attr('stroke-width', d => d.no ? 2.5 : 1)
    .attr('cursor', 'pointer')
    .on('mouseover', showTooltip)
    .on('mousemove', moveTooltip)
    .on('mouseout', hideTooltip)
    .on('click', (event, d) => { highlightCluster(d.cl); event.stopPropagation(); });

  svg.on('click', () => { highlightCluster(null); });

  simulation = d3.forceSimulation(nodes)
    .force('link', d3.forceLink(edges).id(d => d.id).distance(55).strength(.8))
    .force('charge', d3.forceManyBody().strength(-120))
    .force('center', d3.forceCenter(W/2, H/2))
    .force('collide', d3.forceCollide(10))
    .on('tick', () => {
      linkSel.attr('x1', d => d.source.x).attr('y1', d => d.source.y)
             .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
      nodeSel.attr('cx', d => d.x).attr('cy', d => d.y);
    });

  // Store selections for highlighting
  window._nodeSel = nodeSel;
  window._linkSel = linkSel;
  window._nodes = nodes;
  window._edges = edges;
  document.getElementById('graph-overlay').textContent = `Scroll to zoom · Drag to pan · Click node to highlight cluster`;
}

// ── Tooltip ────────────────────────────────────────────────────────────────────
const tip = document.getElementById('tooltip');
function showTooltip(event, d) {
  const ds = ALL_DATA[currentKey];
  const cl = ds.clusters[d.cl] || {};
  let html = `<b>${d.id}</b><br>Original: ${d.or}<br>Cohort: ${d.co}`;
  html += `<br>Cluster #${d.cl+1} (${cl.sz || '?'} members)`;
  if (cl.mp) html += `<br>Max PI_HAT: ${cl.mp}`;
  if (cl.mk) html += `<br>Max KING: ${cl.mk}`;
  if (d.no) html += `<br><span style="color:#ffd54f">★ ${d.no}</span>`;
  tip.innerHTML = html;
  tip.style.display = 'block';
}
function moveTooltip(event) {
  tip.style.left = (event.clientX + 14) + 'px';
  tip.style.top  = (event.clientY - 28) + 'px';
}
function hideTooltip() { tip.style.display = 'none'; }

// ── Cluster highlight ─────────────────────────────────────────────────────────
function highlightCluster(ci) {
  highlightedCluster = ci;
  if (window._nodeSel) {
    window._nodeSel
      .attr('opacity', ci === null ? 1 : d => d.cl === ci ? 1 : .12)
      .attr('r', d => d.no ? 9 : 6);
    window._linkSel
      .attr('opacity', ci === null ? .7 : d => {
        const sn = window._nodes.find(n => n.id === (d.source.id||d.source));
        return (sn && sn.cl === ci) ? .9 : .05;
      });
  }
  // Highlight table row
  document.querySelectorAll('#cluster-tbody tr').forEach(tr => {
    tr.classList.toggle('highlighted', ci !== null && parseInt(tr.dataset.ci) === ci);
  });
  if (ci !== null) {
    const row = document.querySelector(`#cluster-tbody tr[data-ci="${ci}"]`);
    if (row) row.scrollIntoView({block:'nearest', behavior:'smooth'});
  }
}

// ── Table ──────────────────────────────────────────────────────────────────────
function renderTable(ds) {
  const tbody = document.getElementById('cluster-tbody');
  const sorted = [...ds.clusters].sort((a, b) => sortDir * (a[sortCol] > b[sortCol] ? 1 : a[sortCol] < b[sortCol] ? -1 : 0));
  tbody.innerHTML = '';
  sorted.forEach(c => {
    const tr = document.createElement('tr');
    tr.dataset.ci = c.id;
    tr.addEventListener('click', () => highlightCluster(highlightedCluster === c.id ? null : c.id));
    const members = c.mb.slice(0, 8).join('<br>') + (c.mb.length > 8 ? `<br><i>+${c.mb.length-8} more</i>` : '');
    const flags = [
      c.pc ? '<span class="badge-pc">PCRel</span>' : '',
      c.ra ? '<span class="badge-ra">RelAdmix</span>' : '',
      c.no ? `<span class="badge-no">${c.no}</span>` : ''
    ].filter(Boolean).join(' ');
    tr.innerHTML = `<td>${c.id+1}</td><td><b>${c.sz}</b></td>
      <td>${c.mp || '—'}</td><td>${c.mk || '—'}</td>
      <td><div class="member-list">${members}</div></td>
      <td>${flags}</td>`;
    tbody.appendChild(tr);
  });
}

// ── Load threshold ────────────────────────────────────────────────────────────
function loadThreshold(key) {
  currentKey = key;
  highlightedCluster = null;
  const ds = ALL_DATA[key];

  // Warn box
  const wb = document.getElementById('warn-box');
  if (ds.warn) { wb.textContent = ds.warn; wb.style.display = ''; }
  else { wb.style.display = 'none'; }

  // Stats
  document.getElementById('s-nodes').textContent = ds.nn.toLocaleString();
  document.getElementById('s-edges').textContent = ds.ne.toLocaleString();
  document.getElementById('s-clusters').textContent = ds.nc.toLocaleString();
  const largest = ds.clusters.length ? Math.max(...ds.clusters.map(c => c.sz)) : 0;
  document.getElementById('s-largest').textContent = largest;

  renderGraph(ds);
  renderTable(ds);

  document.querySelectorAll('.btn-thresh').forEach(b => b.classList.toggle('active', b.dataset.key === key));
}

// ── Table sort ────────────────────────────────────────────────────────────────
document.querySelectorAll('#cluster-table th[data-col]').forEach(th => {
  th.addEventListener('click', () => {
    const col = th.dataset.col;
    if (sortCol === col) sortDir *= -1; else { sortCol = col; sortDir = -1; }
    renderTable(ALL_DATA[currentKey]);
  });
});

// ── Button clicks ─────────────────────────────────────────────────────────────
document.querySelectorAll('.btn-thresh').forEach(b => {
  b.addEventListener('click', () => loadThreshold(b.dataset.key));
});

// ── Init ──────────────────────────────────────────────────────────────────────
loadThreshold('t3');
</script>
</body>
</html>
"""

out_path = os.path.join(BASE, 'steps', 'relatedness_clusters.html')
with open(out_path, 'w', encoding='utf-8') as f:
    f.write(HTML)
print(f"\nWritten: {out_path}")
