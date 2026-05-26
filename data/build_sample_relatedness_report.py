"""
build_sample_relatedness_report.py
Generates steps/sample_relatedness_report.html:
  - Per-sample table: closest match by PI_HAT, KING, RelateAdmix, PC-Relate
  - Status flags and prominent-case highlighting
  - Sortable/filterable with inline JS

Data sources (all local):
  dragen_qc_data/cross_gwas2026_alsu_winter/cross_ibd_pihat_0125.genome
  dragen_qc_data/cross_gwas2026_alsu_winter/cross_king_00884.kin0
  dragen_qc_data/cross_gwas2026_alsu_winter/sample_map.tsv
  data/relateadmix_top_pairs.tsv

PC-Relate per-pair data is not locally available; the 5 PI_HAT-only PC-Relate-supported
pairs and their kinship values are hardcoded from the admixture_aware_relatedness.html report.
The 43 KING+PI_HAT concordant pairs with PC-Relate >= 0.0884 are noted but individual
pair identities are unknown (counted in the summary only).
"""

import os, sys, math, json, html as html_lib
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

GENOME_FILE  = os.path.join(ROOT, "dragen_qc_data", "cross_gwas2026_alsu_winter", "cross_ibd_pihat_0125.genome")
KING_FILE    = os.path.join(ROOT, "dragen_qc_data", "cross_gwas2026_alsu_winter", "cross_king_00884.kin0")
SMAP_FILE    = os.path.join(ROOT, "dragen_qc_data", "cross_gwas2026_alsu_winter", "sample_map.tsv")
RLA_FILE     = os.path.join(ROOT, "data", "relateadmix_top_pairs.tsv")
OUT_FILE     = os.path.join(ROOT, "steps", "sample_relatedness_report.html")

# ── PC-Relate known pairs (from admixture_aware_relatedness.html) ──────────────
# These 5 are PI_HAT-only discordant but confirmed by PC-Relate >= 0.0884.
# The 43 KING+PI_HAT concordant pairs with PC-Relate >= 0.0884 are not listed
# individually in the local data.
PC_RELATE_PIHAT_ONLY = [
    ("ALSU_0407", "ALSU_0832", 0.5000, 0.795545, "PI_HAT-only confirmed by PC-Relate"),
    ("ALSU_0808", "ALSU_0832", 0.5000, 0.575325, "PI_HAT-only confirmed by PC-Relate"),
    ("ALSU_0787", "ALSU_0832", 0.5000, 0.543498, "PI_HAT-only confirmed by PC-Relate"),
    ("ALSU_0832", "ALSU_0855", 0.5000, 0.157009, "PI_HAT-only confirmed by PC-Relate"),
    ("ALSU_0218", "ALSU_0242", 0.5361, 0.153170, "PI_HAT-only confirmed by PC-Relate AND RelateAdmix"),
]
# Build lookup: frozenset -> pc_relate_kin
PCREL_LOOKUP = {}
for i1, i2, pi, kin, note in PC_RELATE_PIHAT_ONLY:
    PCREL_LOOKUP[frozenset([i1, i2])] = (i1, i2, pi, kin, note)

# ── 1. Sample map ─────────────────────────────────────────────────────────────
smap = {}  # safe_iid -> {cohort, original_iid}
with open(SMAP_FILE, encoding="utf-8-sig") as f:
    hdr = f.readline().strip().split("\t")
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < len(hdr):
            continue
        row = dict(zip(hdr, parts))
        smap[row["safe_iid"]] = {"cohort": row["cohort"], "original_iid": row["original_iid"]}

def orig(iid):
    return smap.get(iid, {}).get("original_iid", iid)

def cohort(iid):
    return smap.get(iid, {}).get("cohort", "?")

# ── 2. PI_HAT ─────────────────────────────────────────────────────────────────
# per_sample_pihat[iid] = (best_match_iid, pi_hat_value)
per_sample_pihat = {}  # best max PI_HAT per sample
all_pihat_pairs = []   # all pairs >= 0.185 for concordance analysis

with open(GENOME_FILE, encoding="utf-8-sig") as f:
    hdr = f.readline().split()
    pi_idx = hdr.index("PI_HAT")
    i1_idx = hdr.index("IID1")
    i2_idx = hdr.index("IID2")
    for line in f:
        parts = line.split()
        if len(parts) <= pi_idx:
            continue
        try:
            pi = float(parts[pi_idx])
        except ValueError:
            continue
        i1, i2 = parts[i1_idx], parts[i2_idx]
        # update best per-sample
        for samp, mate in [(i1, i2), (i2, i1)]:
            if samp not in per_sample_pihat or pi > per_sample_pihat[samp][1]:
                per_sample_pihat[samp] = (mate, pi)
        if pi >= 0.185:
            all_pihat_pairs.append((i1, i2, pi))

# ── 3. KING ───────────────────────────────────────────────────────────────────
per_sample_king = {}   # best max KING kinship per sample
king_pair_set = set()  # frozensets of pairs in KING output

with open(KING_FILE, encoding="utf-8-sig") as f:
    raw_hdr = f.readline().strip().lstrip("#")
    hdr = raw_hdr.split()
    # handle if first col is '#FID1'
    if hdr[0].startswith("#"):
        hdr[0] = hdr[0].lstrip("#")
    ki_idx = hdr.index("KINSHIP")
    i1_idx = hdr.index("IID1")
    i2_idx = hdr.index("IID2")
    for line in f:
        parts = line.split()
        if len(parts) <= ki_idx:
            continue
        try:
            kin = float(parts[ki_idx])
        except ValueError:
            continue
        i1, i2 = parts[i1_idx], parts[i2_idx]
        king_pair_set.add(frozenset([i1, i2]))
        for samp, mate in [(i1, i2), (i2, i1)]:
            if samp not in per_sample_king or kin > per_sample_king[samp][1]:
                per_sample_king[samp] = (mate, kin)

# ── 4. RelateAdmix top pairs ──────────────────────────────────────────────────
per_sample_rla = {}  # best RelateAdmix kinship-equiv per sample

with open(RLA_FILE, encoding="utf-8-sig") as f:
    hdr = f.readline().strip().split("\t")
    kin_idx = hdr.index("relateadmix_kin_equiv")
    i1_idx = hdr.index("iid1")
    i2_idx = hdr.index("iid2")
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) <= kin_idx:
            continue
        try:
            kin = float(parts[kin_idx])
        except ValueError:
            continue
        i1, i2 = parts[i1_idx], parts[i2_idx]
        for samp, mate in [(i1, i2), (i2, i1)]:
            if samp not in per_sample_rla or kin > per_sample_rla[samp][1]:
                per_sample_rla[samp] = (mate, kin)

# ── 5. Build per-sample rows ──────────────────────────────────────────────────
# Union of samples from PI_HAT (>=0.185 threshold for concordance) and KING
pihat_samps = set(i for i, j, p in all_pihat_pairs) | set(j for i, j, p in all_pihat_pairs)
king_samps   = set(per_sample_king.keys())
rla_samps    = set(per_sample_rla.keys())
pcrel_samps  = set()
for i1, i2, pi, kin, note in PC_RELATE_PIHAT_ONLY:
    pcrel_samps.add(i1); pcrel_samps.add(i2)

all_samps = pihat_samps | king_samps | rla_samps | pcrel_samps

rows = []

for samp in sorted(all_samps):
    c = cohort(samp)
    orig_id = orig(samp)

    # PI_HAT best
    if samp in per_sample_pihat:
        ph_mate, ph_val = per_sample_pihat[samp]
        ph_str   = f"{ph_val:.4f}"
        ph_match = f"{orig(ph_mate)} ({ph_mate})"
    else:
        ph_mate, ph_val, ph_str, ph_match = None, 0.0, "—", "—"

    # KING best
    if samp in per_sample_king:
        ki_mate, ki_val = per_sample_king[samp]
        ki_str   = f"{ki_val:.4f}"
        ki_match = f"{orig(ki_mate)} ({ki_mate})"
    else:
        ki_mate, ki_val, ki_str, ki_match = None, 0.0, "—", "—"

    # RelateAdmix best (top-100 only)
    if samp in per_sample_rla:
        ra_mate, ra_val = per_sample_rla[samp]
        ra_str   = f"{ra_val:.4f}"
        ra_match = f"{orig(ra_mate)} ({ra_mate})"
    else:
        ra_mate, ra_val, ra_str, ra_match = None, 0.0, "—", "—"

    # PC-Relate annotation — collect ALL confirmed pairs for this sample
    pcrel_pairs = []  # list of (mate, kin, note_str)
    for i1, i2, pi, kin, note in PC_RELATE_PIHAT_ONLY:
        if samp == i1:
            pcrel_pairs.append((i2, kin, f"kin={kin:.4f} vs {orig(i2)} ({i2})"))
        elif samp == i2:
            pcrel_pairs.append((i1, kin, f"kin={kin:.4f} vs {orig(i1)} ({i1})"))
    pcrel_val  = max((k for _, k, _ in pcrel_pairs), default=0.0)
    pcrel_note = "; ".join(n for _, _, n in pcrel_pairs)

    # ── Concordance flags ─────────────────────────────────────────────────────
    # Is PI_HAT best match pair in KING output?
    ph_pair_in_king = (ph_mate is not None) and (frozenset([samp, ph_mate]) in king_pair_set)
    # PI_HAT best >= 0.185?
    ph_above_185 = ph_val >= 0.185
    # KING best >= 0.0884?
    ki_above_0884 = ki_val >= 0.0884

    # Determine status — PC-Relate confirmed takes priority over best-pair concordance
    if pcrel_val > 0:
        if samp in {"ALSU_0218", "ALSU_0242"}:
            status = "PCRel+RelAdmix confirmed"
        else:
            status = "PI_HAT-only, PCRel confirmed"
    elif ph_pair_in_king and ph_above_185:
        if ph_val >= 0.354 and ki_val >= 0.354:
            status = "1st-degree confirmed"
        elif ph_val >= 0.177 and ki_val >= 0.177:
            status = "2nd-degree confirmed"
        else:
            status = "Concordant (KING+PI_HAT)"
    elif ph_above_185 and not ph_pair_in_king:
        status = "PI_HAT-only"
    elif ki_above_0884:
        status = "KING only"
    else:
        status = "No significant pair"

    # Prominent flag
    prominent = ""
    if samp == "ALSU_0832":
        prominent = "ALSU_0832 cluster (4 PCRel pairs)"
    elif samp in {"ALSU_0218", "ALSU_0242"} and pcrel_val > 0:
        prominent = "Both PCRel+RelAdmix confirmed"
    elif samp in {"ALSU_0407", "ALSU_0808", "ALSU_0787", "ALSU_0855"} and pcrel_val > 0:
        prominent = "PCRel confirmed PI_HAT-only"
    elif ph_val >= 0.5 and ki_val >= 0.354:
        prominent = "1st-degree (all methods)"

    rows.append({
        "safe_iid": samp,
        "orig_id":  orig_id,
        "cohort":   c,
        "ph_val":   ph_val,
        "ph_str":   ph_str,
        "ph_match": ph_match,
        "ki_val":   ki_val,
        "ki_str":   ki_str,
        "ki_match": ki_match,
        "ra_val":   ra_val,
        "ra_str":   ra_str,
        "ra_match": ra_match,
        "pcrel_note": pcrel_note,
        "status":   status,
        "prominent": prominent,
    })

# Sort: prominent first, then by PI_HAT desc
def sort_key(r):
    order = {
        "PCRel+RelAdmix confirmed": 0,
        "PI_HAT-only, PCRel confirmed": 1,
        "1st-degree confirmed": 2,
        "2nd-degree confirmed": 3,
        "Concordant (KING+PI_HAT)": 4,
        "PI_HAT-only": 5,
        "KING only": 6,
        "No significant pair": 7,
    }
    return (order.get(r["status"], 9), -r["ph_val"])

rows.sort(key=sort_key)

# ── 6. Summary statistics ─────────────────────────────────────────────────────
n_concordant = sum(1 for r in rows if "Concordant" in r["status"] or "confirmed" in r["status"].lower())
n_pihat_only = sum(1 for r in rows if r["status"] == "PI_HAT-only")
n_1st_degree = sum(1 for r in rows if r["status"] == "1st-degree confirmed")
n_pcrel      = sum(1 for r in rows if "PCRel" in r["status"])

# ── 7. HTML generation ────────────────────────────────────────────────────────
STATUS_COLORS = {
    "PCRel+RelAdmix confirmed":     ("#b71c1c", "#ffebee"),   # dark red / light red
    "PI_HAT-only, PCRel confirmed": ("#e65100", "#fff3e0"),   # orange / light orange
    "1st-degree confirmed":         ("#1b5e20", "#e8f5e9"),   # dark green / light green
    "2nd-degree confirmed":         ("#2e7d32", "#f1f8e9"),   # green
    "Concordant (KING+PI_HAT)":     ("#0d47a1", "#e3f2fd"),   # blue
    "PI_HAT-only":                  ("#f57f17", "#fffde7"),   # amber / light amber
    "KING only":                    ("#4a148c", "#f3e5f5"),   # purple
    "No significant pair":          ("#546e7a", "#f5f5f5"),   # grey
}

def status_badge(status):
    fg, bg = STATUS_COLORS.get(status, ("#333", "#eee"))
    return (f'<span class="badge" style="background:{bg};color:{fg};border:1px solid {fg};">'
            f'{html_lib.escape(status)}</span>')

def val_cell(val, strval):
    if val == 0.0:
        return '<td style="color:#aaa">—</td>'
    if val >= 0.354:
        color = "#b71c1c"
    elif val >= 0.177:
        color = "#e65100"
    elif val >= 0.0884:
        color = "#1565c0"
    else:
        color = "#555"
    return f'<td style="color:{color};font-weight:600">{strval}</td>'

# Build table rows HTML
table_rows_html = []
for r in rows:
    prom_td = (f'<td><span style="background:#fce4ec;color:#c62828;border:1px solid #ef9a9a;'
               f'border-radius:3px;padding:2px 6px;font-size:0.78em;white-space:nowrap;">'
               f'{html_lib.escape(r["prominent"])}</span></td>'
               if r["prominent"] else '<td></td>')

    pcrel_td = (f'<td style="font-size:0.85em;color:#6a1b9a">{html_lib.escape(r["pcrel_note"])}</td>'
                if r["pcrel_note"] else '<td style="color:#aaa">—</td>')

    ra_td_val = val_cell(r["ra_val"], r["ra_str"] + (" *" if r["ra_val"] > 0 else ""))
    # * means "in top-100 only"

    table_rows_html.append(
        f'<tr data-cohort="{r["cohort"]}" data-status="{html_lib.escape(r["status"])}">'
        f'<td><strong>{html_lib.escape(r["orig_id"])}</strong>'
        f'<br><small style="color:#888">{r["safe_iid"]}</small></td>'
        f'<td>{r["cohort"]}</td>'
        + val_cell(r["ph_val"], r["ph_str"])
        + f'<td style="font-size:0.85em">{html_lib.escape(r["ph_match"])}</td>'
        + val_cell(r["ki_val"], r["ki_str"])
        + f'<td style="font-size:0.85em">{html_lib.escape(r["ki_match"])}</td>'
        + ra_td_val
        + f'<td style="font-size:0.85em">{html_lib.escape(r["ra_match"])}</td>'
        + pcrel_td
        + f'<td>{status_badge(r["status"])}</td>'
        + prom_td
        + '</tr>\n'
    )

table_body = "".join(table_rows_html)

# Legend HTML
legend_items = "".join(
    f'<span class="badge" style="background:{bg};color:{fg};border:1px solid {fg};margin:3px 4px;display:inline-block">'
    f'{html_lib.escape(s)}</span>'
    for s, (fg, bg) in STATUS_COLORS.items()
)

html_out = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Per-Sample Relatedness Cross-Method Report — ALSU+GWAS2026</title>
<style>
  *, *::before, *::after {{ box-sizing: border-box; }}
  body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #f5f7fa; color: #222; margin: 0; padding: 0; }}
  .hero {{ background: linear-gradient(135deg,#1a237e,#283593); color: #fff; padding: 28px 32px 22px; }}
  .hero a {{ color: #90caf9; text-decoration: none; font-size: 0.85em; }}
  .hero h1 {{ margin: 8px 0 4px; font-size: 1.7em; }}
  .hero p {{ margin: 0; opacity: 0.85; font-size: 0.95em; }}
  .cards {{ display: flex; flex-wrap: wrap; gap: 14px; padding: 20px 28px 0; }}
  .card {{ background: #fff; border-radius: 6px; box-shadow: 0 1px 4px rgba(0,0,0,.12);
           padding: 14px 20px; min-width: 130px; text-align: center; }}
  .card .num {{ font-size: 1.9em; font-weight: 700; }}
  .card .lbl {{ font-size: 0.8em; color: #666; margin-top: 2px; }}
  .callout {{ background: #fff8e1; border-left: 4px solid #f9a825; padding: 12px 18px;
              margin: 18px 28px; border-radius: 4px; font-size: 0.93em; }}
  .section {{ padding: 10px 28px 24px; }}
  h2 {{ color: #1a237e; border-bottom: 2px solid #e8eaf6; padding-bottom: 6px; margin-top: 24px; }}
  .controls {{ display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-bottom: 14px; }}
  .controls input, .controls select {{
    border: 1px solid #ccc; border-radius: 4px; padding: 6px 10px; font-size: 0.9em; }}
  .controls input {{ width: 240px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; background: #fff;
           box-shadow: 0 1px 4px rgba(0,0,0,.08); }}
  th {{ background: #283593; color: #fff; padding: 9px 10px; text-align: left;
        cursor: pointer; user-select: none; white-space: nowrap; }}
  th:hover {{ background: #3949ab; }}
  td {{ padding: 8px 10px; border-bottom: 1px solid #eceff1; vertical-align: top; }}
  tr:hover td {{ background: #f3f4f6; }}
  .badge {{ border-radius: 3px; padding: 2px 7px; font-size: 0.82em; white-space: nowrap; font-weight: 600; }}
  .note {{ font-size: 0.8em; color: #666; margin: 10px 0 0; }}
  .legend {{ background: #fff; border: 1px solid #e0e0e0; border-radius: 6px;
             padding: 12px 16px; margin: 0 0 18px; }}
  #rowcount {{ font-size: 0.85em; color: #666; margin-bottom: 6px; }}
</style>
</head>
<body>

<div class="hero">
  <a href="../index.html">← Back to Home</a> &nbsp;|&nbsp;
  <a href="admixture_aware_relatedness.html">Admixture-Aware Relatedness Report</a>
  <h1>Per-Sample Relatedness — Cross-Method Comparison</h1>
  <p>All 1,193 merged ALSU+GWAS2026 samples · Closest match by PI_HAT, KING, RelateAdmix, and PC-Relate</p>
</div>

<div class="cards">
  <div class="card"><div class="num" style="color:#1565c0">{len(rows)}</div><div class="lbl">Samples shown</div></div>
  <div class="card"><div class="num" style="color:#c62828">{n_pihat_only}</div><div class="lbl">PI_HAT-only samples</div></div>
  <div class="card"><div class="num" style="color:#1b5e20">{n_concordant}</div><div class="lbl">Method-concordant</div></div>
  <div class="card"><div class="num" style="color:#6a1b9a">5</div><div class="lbl">PCRel-confirmed PI_HAT-only</div></div>
  <div class="card"><div class="num" style="color:#6a1b9a">1</div><div class="lbl">Both PCRel+RelAdmix confirmed</div></div>
</div>

<div class="callout" style="margin-top:18px">
  <strong>How to read this table:</strong>
  For each sample, the <em>closest match</em> (highest score) found by each method is shown.
  <strong>PI_HAT</strong> uses all pairs ≥ 0.125 in the merged ALSU+GWAS2026 dataset.
  <strong>KING</strong> shows only pairs with kinship ≥ 0.0884.
  <strong>RelateAdmix</strong> shows results from the top 100 pairs only (marked with *); — means the sample did not appear in the top 100.
  <strong>PC-Relate</strong> shows the 5 PI_HAT-only pairs confirmed by PC-Relate at kinship ≥ 0.0884 (hardcoded from the detailed report; 43 KING+PI_HAT concordant pairs also have PC-Relate ≥ 0.0884 but are not individually listed here).
  The <em>Status</em> column reflects the concordance between PI_HAT and KING for the sample's best PI_HAT pair.
</div>

<div class="section">
  <h2>Status Legend</h2>
  <div class="legend">{legend_items}</div>

  <h2>Per-Sample Table</h2>
  <div class="controls">
    <input type="text" id="srchInput" placeholder="Search sample ID…" oninput="filterTable()">
    <select id="cohortFilter" onchange="filterTable()">
      <option value="">All cohorts</option>
      <option value="ALSU">ALSU</option>
      <option value="GWAS2026">GWAS2026</option>
    </select>
    <select id="statusFilter" onchange="filterTable()">
      <option value="">All statuses</option>
      {''.join(f'<option value="{s}">{s}</option>' for s in STATUS_COLORS)}
    </select>
    <label><input type="checkbox" id="prominentOnly" onchange="filterTable()"> Prominent cases only</label>
  </div>
  <div id="rowcount"></div>

  <table id="mainTable">
    <thead>
      <tr>
        <th onclick="sortTable(0)">Sample ▲▼</th>
        <th onclick="sortTable(1)">Cohort ▲▼</th>
        <th onclick="sortTable(2)" title="Max PI_HAT from all pairs ≥ 0.125">PI_HAT ▲▼</th>
        <th>PI_HAT Best Match</th>
        <th onclick="sortTable(4)" title="KING kinship from pairs ≥ 0.0884">KING ▲▼</th>
        <th>KING Best Match</th>
        <th onclick="sortTable(6)" title="RelateAdmix kinship-equiv, top-100 only (*)">RelAdmix* ▲▼</th>
        <th>RelAdmix Best Match</th>
        <th>PC-Relate Note</th>
        <th onclick="sortTable(9)">Status ▲▼</th>
        <th>Notable</th>
      </tr>
    </thead>
    <tbody id="tableBody">
{table_body}
    </tbody>
  </table>
  <p class="note">* RelateAdmix shows best match from the top-100 pairs only (711,028 all-pair estimates computed on DRAGEN; see <a href="admixture_aware_relatedness.html">full report</a>).</p>
  <p class="note">PC-Relate column: only 5 PI_HAT-only discordant pairs confirmed at kinship ≥ 0.0884 are listed here (from <a href="admixture_aware_relatedness.html">admixture-aware report</a>). 43 additional KING+PI_HAT concordant samples also had PC-Relate kinship ≥ 0.0884.</p>
</div>

<script>
let sortDir = {{}};

function sortTable(col) {{
  const tbody = document.getElementById('tableBody');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  const asc = !sortDir[col];
  sortDir[col] = asc;
  rows.sort((a, b) => {{
    const ca = a.cells[col]?.textContent.trim() || '';
    const cb = b.cells[col]?.textContent.trim() || '';
    const na = parseFloat(ca), nb = parseFloat(cb);
    if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
    return asc ? ca.localeCompare(cb) : cb.localeCompare(ca);
  }});
  rows.forEach(r => tbody.appendChild(r));
  filterTable();
}}

function filterTable() {{
  const srch = document.getElementById('srchInput').value.toLowerCase();
  const cohort = document.getElementById('cohortFilter').value;
  const status = document.getElementById('statusFilter').value;
  const prominentOnly = document.getElementById('prominentOnly').checked;
  const tbody = document.getElementById('tableBody');
  let visible = 0;
  tbody.querySelectorAll('tr').forEach(tr => {{
    const txt = tr.textContent.toLowerCase();
    const dc = tr.dataset.cohort || '';
    const ds = tr.dataset.status || '';
    const prominent = tr.cells[10]?.textContent.trim();
    let show = true;
    if (srch && !txt.includes(srch)) show = false;
    if (cohort && dc !== cohort) show = false;
    if (status && ds !== status) show = false;
    if (prominentOnly && !prominent) show = false;
    tr.style.display = show ? '' : 'none';
    if (show) visible++;
  }});
  document.getElementById('rowcount').textContent = `Showing ${{visible}} of {len(rows)} samples`;
}}

// Initial count
filterTable();
</script>
</body>
</html>
"""

with open(OUT_FILE, "w", encoding="utf-8") as f:
    f.write(html_out)

print(f"Written: {OUT_FILE}")
print(f"Total samples in table: {len(rows)}")
print(f"  PI_HAT-only: {n_pihat_only}")
print(f"  Concordant/confirmed: {n_concordant}")
print(f"  PCRel-confirmed: {n_pcrel}")
