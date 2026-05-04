#!/usr/bin/env python3
"""Build chr9/chr15 dense multi-population LD locus pages.

The inputs are regional GWAS TSVs and PLINK2 --r2-phased square matrices
created from the UZB post-imputation cohort. 1000G EUR/SAS matrices and
Ensembl gene tracks are fetched once and cached under data/regional_ld/.
"""
from __future__ import annotations

import json
import math
import statistics
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "regional_ld"
OUT_DIR = ROOT / "steps"

PEAK_P = 1e-5
EXTENSION_R2 = 0.60
REFERENCE_WINDOW_KB = 500


@dataclass(frozen=True)
class LocusConfig:
    chrom: int
    lead_id: str
    lead_pos: int
    start: int
    end: int
    display_start: int
    display_end: int
    label: str
    subtitle: str
    context: str
    functional_page: str
    canonical_page: str
    legacy_page: str


LOCI = [
    LocusConfig(
        chrom=9,
        lead_id="rs28446251",
        lead_pos=115_907_201,
        start=115_800_000,
        end=116_050_000,
        display_start=115_750_000,
        display_end=116_420_000,
        label="chr9:115.80-116.05 Mb",
        subtitle="LINC00474/PAPPA region",
        context="Suggestive RPL association overlapping LINC00474 and upstream of PAPPA",
        functional_page="chr9_locus_functional.html",
        canonical_page="chr9_ld_eur_sas_uzb.html",
        legacy_page="chr9_ld_eur_uzb.html",
    ),
    LocusConfig(
        chrom=15,
        lead_id="rs8027539",
        lead_pos=31_340_603,
        start=31_200_000,
        end=31_480_000,
        display_start=31_150_000,
        display_end=31_490_000,
        label="chr15:31.20-31.48 Mb",
        subtitle="KLF13 / lncRNA-rich 31.3 Mb region",
        context="Suggestive RPL association across a compact KLF13-proximal LD block",
        functional_page="chr15_locus_functional.html",
        canonical_page="chr15_ld_eur_sas_uzb.html",
        legacy_page="chr15_ld_eur_uzb.html",
    ),
]


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root {
  --bg:#f5f7fb; --panel:#ffffff; --ink:#172033; --muted:#5f6f84;
  --border:#d4dce8; --soft:#edf2f7; --blue:#1769aa; --teal:#0f766e;
  --orange:#c85f14; --red:#b8323f; --green:#2f7d32; --violet:#6d4aa3;
  --mono:"Courier New",monospace;
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:"Segoe UI",Arial,sans-serif;line-height:1.48;padding:22px}
.wrap{max-width:1260px;margin:0 auto}
a{color:var(--blue);text-decoration:none}
a:hover{text-decoration:underline}
h1{font-size:25px;line-height:1.25;margin:14px 0 7px;letter-spacing:0}
h2{font-size:18px;margin:0 0 12px;letter-spacing:0}
h3{font-size:12px;margin:0 0 6px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
p{margin:8px 0}.sub{color:var(--muted);font-size:14px;margin-bottom:12px}
code{font-family:var(--mono);font-size:12px;background:var(--soft);border-radius:3px;padding:2px 5px}
.status-row{display:flex;flex-wrap:wrap;gap:8px;margin:9px 0 18px}
.pill{font-size:12px;border:1px solid var(--border);background:#fff;border-radius:999px;padding:5px 9px;color:var(--muted)}
.pill.strong{border-color:rgba(23,105,170,.35);color:#154c7b;background:rgba(23,105,170,.08)}
.card{background:var(--panel);border:1px solid var(--border);border-radius:8px;padding:16px;margin:15px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:11px;margin:12px 0}
.metric{background:var(--soft);border-left:4px solid var(--blue);border-radius:6px;padding:12px;min-height:86px}
.metric:nth-child(2){border-left-color:var(--teal)}.metric:nth-child(3){border-left-color:var(--orange)}.metric:nth-child(4){border-left-color:var(--red)}
.metric:nth-child(5){border-left-color:var(--green)}.metric:nth-child(6){border-left-color:var(--violet)}
.metric .v{font:700 18px var(--mono);color:var(--ink)}.metric .s{font-size:12px;color:var(--muted);margin-top:4px}
.canvas-wrap{position:relative;overflow-x:auto;border:1px solid var(--border);background:#fff;border-radius:6px;padding:10px;margin:10px 0}
canvas{display:block;width:100%;min-width:760px;height:auto;background:#fff}
.figure-card{padding-bottom:18px}
.figure-caption{font-size:13px;color:var(--muted);margin:4px 0 10px}
.figure-stack{border:1px solid var(--border);background:#fff;border-radius:6px;overflow-x:auto;margin-top:10px}
.figure-stack .canvas-wrap{border:0;border-radius:0;margin:0;padding:8px 10px;overflow:visible;min-width:900px}
.figure-stack .manhattan-wrap{border-bottom:1px solid var(--border);padding-bottom:0}
.figure-stack .heatmap-wrap{padding-top:0}
.figure-stack canvas{min-width:900px}
#matrixCaption{margin:8px 10px 6px}
.figure-legend{min-width:900px;background:var(--soft);border:1px solid var(--border);border-top:0;border-radius:0 0 6px 6px;padding:8px 12px 10px;font:11px var(--mono);display:flex;flex-direction:column;gap:6px}
.legend-row{display:flex;align-items:center;gap:12px;flex-wrap:wrap}
.legend-row .sep{color:#aab7c8;font-size:14px;line-height:1}
.leg-title{color:var(--muted);font-size:10px;text-transform:uppercase;letter-spacing:.04em;min-width:122px}
.leg-gradient{width:180px;height:13px;border-radius:2px;border:1px solid var(--border);background:linear-gradient(90deg,#f8fafc,#deebf7,#a0cebe,#f5c65c,#dd643e,#ae2442,#57142d)}
.leg-gradient.delta{background:linear-gradient(90deg,#2c5493,#6fabcd,#ebf1f7,#f9f9f7,#fae2b3,#de803c,#972d12)}
.leg-labels{display:flex;justify-content:space-between;width:180px;color:var(--muted);font-size:9px}
.leg-item{display:flex;align-items:center;gap:4px;color:var(--ink)}
.leg-sw{width:10px;height:10px;border-radius:2px;border:1px solid var(--border);display:inline-block;flex:0 0 auto}
.leg-dot{width:14px;height:14px;vertical-align:middle;flex:0 0 auto}
.leg-line{display:inline-block;width:20px;height:0;border-top:1.5px dashed #b8323f;vertical-align:middle}
.leg-line.gws{border-top-color:#7f1d1d}.leg-line.suggestive{border-top-color:#c85f14}
.controls{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin:8px 0 10px}
button{border:1px solid #aab7c8;background:#fff;border-radius:6px;padding:7px 10px;cursor:pointer;color:var(--ink);font-weight:650}
button.active{background:var(--blue);border-color:var(--blue);color:#fff}
.small{font-size:12px;color:var(--muted)}
.note{border-left:4px solid var(--blue);background:rgba(23,105,170,.08);border-radius:5px;padding:11px;margin:12px 0;font-size:13px}
.warn{border-left:4px solid var(--orange);background:rgba(200,95,20,.09);border-radius:5px;padding:11px;margin:12px 0;font-size:13px}
table{width:100%;border-collapse:collapse;margin-top:10px;font-size:13px}
th,td{border:1px solid var(--border);padding:7px 8px;text-align:left;vertical-align:top}
th{background:var(--soft)}
.two-col{display:grid;grid-template-columns:1fr 1fr;gap:14px}
.tooltip{position:fixed;display:none;pointer-events:none;background:#172033;color:#fff;border-radius:6px;padding:8px 9px;font-size:12px;z-index:20;max-width:330px;box-shadow:0 10px 25px rgba(0,0,0,.18)}
@media (max-width:820px){body{padding:13px}h1{font-size:21px}.card{padding:13px}.two-col{grid-template-columns:1fr}canvas{min-width:700px}}
</style>
</head>
<body>
<div class="wrap">
  <a href="step16.html">&lt; Back to GWAS Results</a>
  <h1 id="pageTitle"></h1>
  <div class="sub" id="pageSub"></div>
  <div class="status-row" id="badges"></div>

  <div class="card">
    <h2>Dense SNP-Level LD Summary</h2>
    <p id="summaryText"></p>
    <div class="grid" id="metricsGrid"></div>
    <div class="note">UZB cohort LD is from the 1090-sample post-imputation dataset using <code>plink2 --r2-phased square</code>. EUR and SAS reference LD are from Ensembl REST 1000 Genomes phase 3 and are cached with this page build.</div>
  </div>

  <div class="card">
    <h2>Regional GWAS Context</h2>
    <p id="regionalContextText"></p>
    <p class="small" id="regionalContextDetail"></p>
  </div>

  <div class="card figure-card">
    <h2>Zoomed Regional Signal and LD Block</h2>
    <p class="figure-caption" id="figureCaption"></p>
    <div class="controls" aria-label="LD heatmap mode">
      <button type="button" data-mode="eur">EUR 1000G</button>
      <button type="button" data-mode="sas">SAS 1000G</button>
      <button type="button" data-mode="uzb" class="active">Uzbek Cohort</button>
      <button type="button" data-mode="eur_delta">Delta EUR-UZB</button>
      <button type="button" data-mode="sas_delta">Delta SAS-UZB</button>
    </div>
    <p class="small" id="matrixCaption"></p>
    <div class="figure-stack">
      <div class="canvas-wrap manhattan-wrap"><canvas id="manhattanCanvas" width="1120" height="330"></canvas></div>
      <div class="canvas-wrap heatmap-wrap"><canvas id="heatmapCanvas" width="1120" height="760"></canvas></div>
      <div class="figure-legend" id="figureLegend">
        <div class="legend-row">
          <span class="leg-title" id="legendHeatTitle">heatmap r2</span>
          <span class="leg-gradient" id="legendGradient"></span>
          <span class="leg-labels" id="legendHeatLabels"><span>0</span><span>0.25</span><span>0.5</span><span>0.75</span><span>1</span></span>
          <span class="sep">|</span>
          <span class="leg-item"><span class="leg-sw" style="background:#f0f3f8"></span> n/a or missing pair</span>
        </div>
        <div class="legend-row">
          <span class="leg-title">Manhattan dot color</span>
          <span class="leg-item"><span class="leg-sw" style="background:#57142d"></span> r2 to lead &ge;0.75</span>
          <span class="leg-item"><span class="leg-sw" style="background:#ae2442"></span> 0.50-0.75</span>
          <span class="leg-item"><span class="leg-sw" style="background:#dd643e"></span> 0.25-0.50</span>
          <span class="leg-item"><span class="leg-sw" style="background:#a0cebe"></span> 0.03-0.25</span>
          <span class="leg-item"><span class="leg-sw" style="background:#f8fafc"></span> &lt;0.03 / no r2</span>
        </div>
        <div class="legend-row">
          <span class="leg-title">Manhattan dots</span>
          <span class="leg-item"><svg class="leg-dot" viewBox="0 0 14 14"><circle cx="7" cy="7" r="5.2" fill="#b8323f" stroke="#172033" stroke-width="1.4"/></svg> lead <b id="legendLeadId"></b></span>
          <span class="leg-item"><svg class="leg-dot" viewBox="0 0 14 14"><circle cx="7" cy="7" r="4" fill="#dd643e" stroke="#172033" stroke-width="1.6"/></svg> highlighted LD-set SNPs (<span id="legendSelectedCount"></span>)</span>
          <span class="leg-item"><svg class="leg-dot" viewBox="0 0 14 14"><circle cx="7" cy="7" r="2.5" fill="#a0cebe" stroke="#ffffff" stroke-width=".8"/></svg> other regional GWAS rows</span>
        </div>
        <div class="legend-row">
          <span class="leg-title">Thresholds</span>
          <span class="leg-item"><span class="leg-line gws"></span> genome-wide 5e-8 (-log10 p=7.30)</span>
          <span class="leg-item"><span class="leg-line suggestive"></span> suggestive 1e-5 (-log10 p=5.00)</span>
          <span class="leg-item">selected set: <b id="legendPeakCount"></b> peak + <b id="legendExtensionCount"></b> extension variants</span>
        </div>
        <div class="legend-row">
          <span class="leg-title">Genes</span>
          <span class="leg-item"><span class="leg-sw" style="background:#2f7d32"></span> protein-coding</span>
          <span class="leg-item"><span class="leg-sw" style="background:#6d4aa3"></span> lncRNA</span>
          <span class="leg-item"><span class="leg-sw" style="background:#1769aa"></span> other annotation</span>
          <span class="leg-item" style="color:var(--muted)">triangle at gene end = transcription strand</span>
        </div>
        <div class="legend-row">
          <span class="leg-title">Connectors</span>
          <span class="leg-item" style="color:var(--muted);font-size:10.5px">Top heatmap ticks mark each SNP's genomic position; curved connectors map those ticks to uniformly spaced LD matrix columns. These ticks align to the Manhattan x-axis above.</span>
        </div>
        <div class="legend-row" style="font-size:10px;color:var(--muted)">
          <span class="leg-item" id="legendDataCaveat"></span>
        </div>
      </div>
    </div>
  </div>

  <div class="two-col">
    <div class="card">
      <h2>EUR-UZB Discordant Pairs</h2>
      <table id="eurDiscord"><thead><tr><th>Pair</th><th>EUR</th><th>UZB</th><th>|Delta|</th></tr></thead><tbody></tbody></table>
    </div>
    <div class="card">
      <h2>SAS-UZB Discordant Pairs</h2>
      <table id="sasDiscord"><thead><tr><th>Pair</th><th>SAS</th><th>UZB</th><th>|Delta|</th></tr></thead><tbody></tbody></table>
    </div>
  </div>

  <div class="card">
    <h2>Highlighted SNPs</h2>
    <p class="small">Rows include all p&lt;1e-5 regional peak variants and extension variants with UZB r<sup>2</sup> to the lead at or above the extension threshold.</p>
    <table id="variantTable"><thead><tr><th>SNP</th><th>Position</th><th>Role</th><th>Firth p</th><th>OR</th><th>UZB r2 to lead</th><th>EUR r2</th><th>SAS r2</th><th>SAIGE p</th></tr></thead><tbody></tbody></table>
  </div>

  <div class="card">
    <h2>Interpretation Notes</h2>
    <div class="warn">These loci are suggestive, not genome-wide significant, in the primary Firth GWAS. The page is meant to expose LD structure, population differences, and local annotation for follow-up prioritization.</div>
    <p>The chr9 and chr15 pages now use the same dense SNP-level workflow as the chr18 reference: primary GWAS context, UZB cohort LD, EUR/SAS 1000G comparison, delta matrices, discordant-pair summaries, and a locus-level gene track.</p>
    <p><a id="functionalLink"></a> | <a href="../chr18/chr18_ld_eur_sas_uzb.html">chr18 dense LD reference page</a></p>
  </div>
</div>

<div class="tooltip" id="tooltip"></div>
<script id="payload" type="application/json">__PAYLOAD__</script>
<script>
const DATA = JSON.parse(document.getElementById('payload').textContent);
let heatMode = 'uzb';
const tooltip = document.getElementById('tooltip');

function byId(id) { return document.getElementById(id); }
function fmtBp(value) { return Number(value).toLocaleString('en-US'); }
function fmtMb(value) { return (value / 1e6).toFixed(3) + ' Mb'; }
function fmtP(value) {
  if (value === null || value === undefined || !Number.isFinite(value)) return 'NA';
  if (value === 0) return '0';
  if (value < 0.001) return value.toExponential(2);
  return value.toPrecision(3);
}
function fmtNum(value, digits = 3) {
  if (value === null || value === undefined || !Number.isFinite(value)) return 'NA';
  return Number(value).toFixed(digits);
}
function minusLog10(value) { return value > 0 ? -Math.log10(value) : 0; }
function clamp(value, lo, hi) { return Math.max(lo, Math.min(hi, value)); }
function figureWindow() { return DATA.figure_window || DATA.window || DATA.display_window; }
function formatWindowMb(windowRange) { return fmtMb(windowRange[0]) + '-' + fmtMb(windowRange[1]); }
function modeLabel() { return {eur:'EUR 1000G r2', sas:'SAS 1000G r2', uzb:'Uzbek cohort r2', eur_delta:'EUR - UZB r2', sas_delta:'SAS - UZB r2'}[heatMode]; }

function interpolate(stops, value) {
  if (value === null || value === undefined || !Number.isFinite(value)) return '#f0f3f8';
  const v = clamp(value, stops[0][0], stops[stops.length - 1][0]);
  let left = stops[0], right = stops[stops.length - 1];
  for (let i = 0; i < stops.length - 1; i++) {
    if (v >= stops[i][0] && v <= stops[i + 1][0]) { left = stops[i]; right = stops[i + 1]; break; }
  }
  const span = right[0] - left[0] || 1;
  const t = (v - left[0]) / span;
  const rgb = [1, 2, 3].map(i => Math.round(left[i] + (right[i] - left[i]) * t));
  return 'rgb(' + rgb.join(',') + ')';
}
function r2Color(value) {
  return interpolate([
    [0.00, 248, 250, 252], [0.03, 222, 235, 247], [0.10, 160, 206, 190],
    [0.25, 245, 198, 92], [0.50, 221, 100, 62], [0.75, 174, 36, 66], [1.00, 87, 20, 45]
  ], value);
}
function deltaColor(value) {
  return interpolate([
    [-0.70, 44, 84, 147], [-0.25, 111, 171, 205], [-0.04, 235, 241, 247],
    [0.00, 249, 249, 247], [0.04, 250, 226, 179], [0.25, 222, 128, 60], [0.70, 151, 45, 18]
  ], value);
}

function metric(label, value, sub) {
  return '<div class="metric"><h3>' + label + '</h3><div class="v">' + value + '</div><div class="s">' + sub + '</div></div>';
}
function setHeader() {
  byId('pageTitle').textContent = DATA.header;
  byId('pageSub').textContent = DATA.subtitle;
  byId('summaryText').textContent = DATA.context;
  byId('regionalContextText').textContent = DATA.context + '. The full regional GWAS input spans ' + formatWindowMb(DATA.display_window) + ' and contains ' + DATA.counts.regional_gwas + ' Firth rows used for context and filtering.';
  byId('regionalContextDetail').textContent = 'The figure below is deliberately zoomed to ' + formatWindowMb(figureWindow()) + ', the dense LD matrix window. Its Manhattan x-axis, heatmap genomic ruler, and SNP connector ticks use the same genomic scale.';
  byId('figureCaption').textContent = 'Zoomed locus figure: regional Firth GWAS and gene track on top, then the matched multi-population LD block immediately underneath on the same genomic scale.';
  byId('badges').innerHTML = DATA.badges.map((text, i) => '<span class="pill ' + (i === 0 ? 'strong' : '') + '">' + text + '</span>').join('');
  byId('legendLeadId').textContent = DATA.lead.id;
  byId('legendSelectedCount').textContent = DATA.selected.length;
  byId('legendPeakCount').textContent = DATA.counts.peak;
  byId('legendExtensionCount').textContent = DATA.counts.extension;
  byId('legendDataCaveat').textContent = 'UZB LD: 1090 post-imputation samples with plink2 --r2-phased square. EUR/SAS LD: Ensembl REST 1000 Genomes phase 3; blank heatmap cells mean unavailable reference or cohort pairs, and delta modes are blank when either source is missing.';
  byId('metricsGrid').innerHTML = [
    metric('Highlighted SNPs', DATA.selected.length, DATA.counts.peak + ' peak + ' + DATA.counts.extension + ' extension'),
    metric('Lead p / OR', fmtP(DATA.lead.p) + ' / ' + fmtNum(DATA.lead.or, 2), DATA.lead.id + ' at ' + fmtBp(DATA.lead.pos)),
    metric('UZB mean r2', fmtNum(DATA.stats.uzb.mean_pair_r2, 3), DATA.stats.uzb.coverage + ' observed pair cells'),
    metric('EUR vs UZB', fmtNum(DATA.comparisons.eur_uzb.pearson, 3), 'mean |delta| ' + fmtNum(DATA.comparisons.eur_uzb.mean_abs_delta, 3)),
    metric('SAS vs UZB', fmtNum(DATA.comparisons.sas_uzb.pearson, 3), 'mean |delta| ' + fmtNum(DATA.comparisons.sas_uzb.mean_abs_delta, 3)),
    metric('SAIGE support', DATA.saige.best ? fmtP(DATA.saige.best.p) : 'NA', DATA.saige.best ? DATA.saige.best.id : 'no regional row')
  ].join('');
  const link = byId('functionalLink');
  link.href = DATA.functional_page;
  link.textContent = DATA.functional_label;
}

function drawAxes(ctx, x, y, plot, xTicks, yTicks) {
  ctx.strokeStyle = '#b8c3d1';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(plot.left, plot.top);
  ctx.lineTo(plot.left, plot.bottom);
  ctx.lineTo(plot.right, plot.bottom);
  ctx.stroke();
  ctx.fillStyle = '#5f6f84';
  ctx.font = '12px Segoe UI, Arial';
  ctx.textAlign = 'center';
  xTicks.forEach(t => {
    const px = x(t);
    ctx.strokeStyle = '#edf2f7';
    ctx.beginPath(); ctx.moveTo(px, plot.top); ctx.lineTo(px, plot.bottom); ctx.stroke();
    ctx.fillStyle = '#5f6f84'; ctx.fillText(fmtMb(t), px, plot.bottom + 18);
  });
  ctx.textAlign = 'right';
  yTicks.forEach(t => {
    const py = y(t);
    ctx.strokeStyle = '#edf2f7';
    ctx.beginPath(); ctx.moveTo(plot.left, py); ctx.lineTo(plot.right, py); ctx.stroke();
    ctx.fillStyle = '#5f6f84'; ctx.fillText(t.toFixed(0), plot.left - 8, py + 4);
  });
}

function drawManhattan() {
  const canvas = byId('manhattanCanvas');
  const ctx = canvas.getContext('2d');
  const w = canvas.width, h = canvas.height;
  ctx.clearRect(0, 0, w, h); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h);
  const figure = figureWindow();
  const x0 = figure[0], x1 = figure[1];
  const plot = {left:72, right:w - 30, top:118, bottom:260};
  const visibleRows = DATA.gwas.filter(row => row.pos >= x0 && row.pos <= x1);
  const plotRows = visibleRows.length ? visibleRows : DATA.gwas;
  const maxVisibleLogP = Math.max(...plotRows.map(d => minusLog10(d.p)), minusLog10(DATA.lead.p));
  const maxY = Math.max(7.5, Math.ceil(maxVisibleLogP + 0.5));
  const x = pos => plot.left + ((pos - x0) / (x1 - x0)) * (plot.right - plot.left);
  const y = val => plot.bottom - (val / maxY) * (plot.bottom - plot.top);
  const ticks = [];
  const step = Math.max(25_000, Math.round((x1 - x0) / 5 / 5_000) * 5_000);
  for (let t = Math.ceil(x0 / step) * step; t <= x1; t += step) ticks.push(t);
  ctx.fillStyle = '#172033'; ctx.font = '700 13px Segoe UI, Arial'; ctx.textAlign = 'left';
  ctx.fillText('Zoomed regional GWAS and gene track', plot.left, 20);
  ctx.fillStyle = '#5f6f84'; ctx.font = '12px Segoe UI, Arial';
  ctx.fillText('Window ' + formatWindowMb(figure) + ' (GRCh38); this x-scale is reused by the LD ruler below.', plot.left, 39);
  drawGeneTrack(ctx, x, 70, plot.left, plot.right);
  const yTicks = Array.from(new Set([0, 2, 4, 6, Math.floor(maxY)].filter(t => t <= maxY)));
  drawAxes(ctx, x, y, plot, ticks, yTicks);
  ctx.setLineDash([5, 4]);
  ctx.strokeStyle = '#7f1d1d'; ctx.beginPath(); ctx.moveTo(plot.left, y(-Math.log10(5e-8))); ctx.lineTo(plot.right, y(-Math.log10(5e-8))); ctx.stroke();
  ctx.strokeStyle = '#c85f14'; ctx.beginPath(); ctx.moveTo(plot.left, y(5)); ctx.lineTo(plot.right, y(5)); ctx.stroke(); ctx.setLineDash([]);
  plotRows.forEach(row => {
    if (row.pos < x0 || row.pos > x1) return;
    const px = x(row.pos), py = y(minusLog10(row.p));
    const isSelected = row.selected;
    ctx.fillStyle = r2Color(row.r2_to_lead || 0);
    ctx.strokeStyle = isSelected ? '#172033' : '#ffffff';
    ctx.lineWidth = isSelected ? 1.1 : 0.8;
    ctx.beginPath(); ctx.arc(px, py, row.id === DATA.lead.id ? 5.2 : (isSelected ? 4.0 : 2.5), 0, Math.PI * 2); ctx.fill(); ctx.stroke();
  });
  ctx.strokeStyle = '#172033'; ctx.lineWidth = 1.3; ctx.setLineDash([3, 3]);
  const leadX = x(DATA.lead.pos); ctx.beginPath(); ctx.moveTo(leadX, plot.top); ctx.lineTo(leadX, plot.bottom + 84); ctx.stroke(); ctx.setLineDash([]);
  const labels = DATA.selected.filter(v => v.pos >= x0 && v.pos <= x1 && (v.role === 'lead' || v.p < 1e-6 || v.id === 'rs3037402')).slice(0, 12);
  ctx.font = '11px Segoe UI, Arial'; ctx.fillStyle = '#172033'; ctx.textAlign = 'left';
  labels.forEach((v, i) => {
    const px = x(v.pos) + 5;
    const py = y(minusLog10(v.p)) - 7 - ((i % 3) * 11);
    ctx.fillText(v.id, Math.min(px, plot.right - 78), Math.max(plot.top + 11, py));
  });
  ctx.save(); ctx.translate(18, (plot.top + plot.bottom) / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = '#5f6f84'; ctx.textAlign = 'center'; ctx.fillText('-log10(Firth p)', 0, 0); ctx.restore();
  ctx.textAlign = 'center'; ctx.fillStyle = '#5f6f84'; ctx.fillText('Genomic position (Mb; shared with LD block below)', (plot.left + plot.right) / 2, plot.bottom + 42);
}

function drawGeneTrack(ctx, x, yBase, left, right) {
  const figure = figureWindow();
  const genes = DATA.genes.filter(gene => gene.end >= figure[0] && gene.start <= figure[1]);
  ctx.font = '11px Segoe UI, Arial';
  ctx.fillStyle = '#5f6f84'; ctx.textAlign = 'left'; ctx.fillText('Genes / lncRNAs', left, yBase - 18);
  if (!genes.length) {
    ctx.fillText('No annotated genes overlap this zoom window', left, yBase + 2);
    return;
  }
  const colors = {protein_coding:'#2f7d32', lncRNA:'#6d4aa3', processed_pseudogene:'#9a6b26'};
  const laneEnds = [];
  genes.forEach((gene) => {
    const xA = clamp(x(gene.start), left, right);
    const xB = clamp(x(gene.end), left, right);
    const geneLeft = Math.min(xA, xB), geneRight = Math.max(xA, xB);
    let lane = 0;
    while (laneEnds[lane] && laneEnds[lane] > geneLeft - 28) lane++;
    lane = Math.min(lane, 3);
    laneEnds[lane] = Math.max(laneEnds[lane] || 0, geneRight);
    const y = yBase + lane * 18;
    const color = colors[gene.biotype] || '#1769aa';
    ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = 3;
    ctx.beginPath(); ctx.moveTo(geneLeft, y); ctx.lineTo(Math.max(geneLeft + 5, geneRight), y); ctx.stroke();
    ctx.beginPath();
    if (gene.strand < 0) { ctx.moveTo(geneLeft, y); ctx.lineTo(geneLeft + 8, y - 4); ctx.lineTo(geneLeft + 8, y + 4); }
    else { ctx.moveTo(geneRight, y); ctx.lineTo(geneRight - 8, y - 4); ctx.lineTo(geneRight - 8, y + 4); }
    ctx.closePath(); ctx.fill();
    ctx.fillStyle = '#172033'; ctx.textAlign = 'left';
    ctx.fillText(gene.name, Math.min(geneLeft + 3, geneRight - 60), y - 6);
  });
}

function getMatrixForMode() {
  if (heatMode === 'eur') return DATA.matrices.eur;
  if (heatMode === 'sas') return DATA.matrices.sas;
  if (heatMode === 'uzb') return DATA.matrices.uzb;
  const ref = heatMode === 'eur_delta' ? DATA.matrices.eur : DATA.matrices.sas;
  return ref.map((row, i) => row.map((value, j) => {
    const base = DATA.matrices.uzb[i][j];
    return value === null || base === null ? null : value - base;
  }));
}
function drawHeatmap() {
  const canvas = byId('heatmapCanvas');
  const ctx = canvas.getContext('2d');
  const w = canvas.width, h = canvas.height;
  ctx.clearRect(0, 0, w, h); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h);
  const vars = DATA.selected;
  const n = vars.length;
  const rulerLeft = 72, rulerRight = w - 30;
  const figure = figureWindow();
  const cell = Math.min(17, Math.max(12, (w - 260) / n));
  const axisLen = n * cell;
  const axisLeft = Math.max(120, (w - axisLen) / 2);
  const axisTop = 148;
  const diamondTop = axisTop + 7;
  const heatHeight = Math.ceil(axisLen / 2) + cell;
  const labelTop = 118;
  const rulerY = 24;
  const stemTopY = axisTop - 12;
  const mat = getMatrixForMode();
  const label = modeLabel();
  const xGenome = pos => rulerLeft + ((pos - figure[0]) / (figure[1] - figure[0])) * (rulerRight - rulerLeft);
  const colX = index => axisLeft + (index + 0.5) * cell;
  const colorForVariant = v => v.role === 'lead' ? '#b8323f' : r2Color(v.r2_uzb_lead);

  ctx.strokeStyle = '#3a4868'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(rulerLeft, rulerY); ctx.lineTo(rulerRight, rulerY); ctx.stroke();
  ctx.fillStyle = '#5f6f84'; ctx.font = '11px "Courier New", monospace'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
  ctx.fillText('genomic scale shared with Manhattan', rulerLeft, rulerY - 12);
  const tickStep = Math.max(25_000, Math.round((figure[1] - figure[0]) / 5 / 5_000) * 5_000);
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (let tick = Math.ceil(figure[0] / tickStep) * tickStep; tick <= figure[1]; tick += tickStep) {
    const x = xGenome(tick);
    ctx.strokeStyle = '#7c8ba3';
    ctx.beginPath(); ctx.moveTo(x, rulerY - 4); ctx.lineTo(x, rulerY + 5); ctx.stroke();
    ctx.fillStyle = '#5f6f84'; ctx.fillText((tick / 1e6).toFixed(2), x, rulerY + 9);
  }

  vars.forEach((v, index) => {
    const genomicX = xGenome(v.pos);
    const columnX = colX(index);
    const color = colorForVariant(v);
    const isLead = v.id === DATA.lead.id;
    ctx.strokeStyle = color; ctx.lineWidth = isLead ? 1.9 : 1;
    ctx.beginPath(); ctx.moveTo(genomicX, rulerY - 6); ctx.lineTo(genomicX, rulerY + 6); ctx.stroke();
    ctx.globalAlpha = isLead ? 1 : 0.55;
    ctx.lineWidth = isLead ? 1.7 : 0.75;
    ctx.beginPath();
    ctx.moveTo(genomicX, rulerY + 7);
    ctx.bezierCurveTo(genomicX, 70, columnX, 78, columnX, stemTopY);
    ctx.stroke();
    ctx.globalAlpha = 1;
    ctx.lineWidth = isLead ? 1.7 : 0.9;
    ctx.beginPath(); ctx.moveTo(columnX, stemTopY); ctx.lineTo(columnX, axisTop - 2); ctx.stroke();
    ctx.save();
    ctx.translate(columnX, labelTop);
    ctx.rotate(-Math.PI / 2);
    ctx.font = isLead ? '700 9px "Courier New", monospace' : '8px "Courier New", monospace';
    ctx.fillStyle = isLead ? '#b8323f' : (v.role === 'peak' ? '#1769aa' : '#4b5f78');
    ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    ctx.fillText(v.id, 0, 0);
    ctx.restore();
    if (isLead) {
      ctx.fillStyle = '#b8323f';
      ctx.beginPath();
      ctx.moveTo(columnX, stemTopY - 6); ctx.lineTo(columnX - 5, stemTopY); ctx.lineTo(columnX, stemTopY + 5); ctx.lineTo(columnX + 5, stemTopY);
      ctx.closePath(); ctx.fill();
    }
  });

  ctx.strokeStyle = '#3a4868'; ctx.lineWidth = 1.2;
  ctx.beginPath(); ctx.moveTo(axisLeft, axisTop); ctx.lineTo(axisLeft + axisLen, axisTop); ctx.stroke();
  ctx.fillStyle = '#5f6f84'; ctx.font = '11px "Courier New", monospace'; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  ctx.fillText('LD block columns', axisLeft - 8, axisTop);
  for (let index = 0; index < n; index++) {
    const x = colX(index);
    ctx.strokeStyle = '#3a4868';
    ctx.beginPath(); ctx.moveTo(x, axisTop - 4); ctx.lineTo(x, axisTop + 4); ctx.stroke();
  }

  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      const centerX = (colX(i) + colX(j)) / 2;
      const centerY = diamondTop + (colX(j) - colX(i)) / 2;
      const half = cell / 2;
      const value = mat[i][j];
      ctx.beginPath();
      ctx.moveTo(centerX, centerY - half); ctx.lineTo(centerX + half, centerY); ctx.lineTo(centerX, centerY + half); ctx.lineTo(centerX - half, centerY);
      ctx.closePath();
      ctx.fillStyle = heatMode.indexOf('delta') >= 0 ? deltaColor(value) : r2Color(value);
      ctx.fill();
      ctx.strokeStyle = (vars[i].id === DATA.lead.id || vars[j].id === DATA.lead.id) ? 'rgba(184,50,63,.45)' : 'rgba(23,32,51,.08)';
      ctx.lineWidth = (vars[i].id === DATA.lead.id || vars[j].id === DATA.lead.id) ? 0.7 : 0.3;
      ctx.stroke();
    }
  }

  const legendY = Math.min(h - 44, diamondTop + heatHeight + 38);
  ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic'; ctx.font = '12px Segoe UI, Arial'; ctx.fillStyle = '#172033';
  ctx.fillText(label + ' triangular LD block', axisLeft, legendY - 12);
  drawLegend(ctx, axisLeft + 205, legendY - 25, heatMode.indexOf('delta') >= 0);
  canvas._geo = {axisLeft, axisTop, diamondTop, heatHeight, cell, vars};
  byId('matrixCaption').textContent = label + ' triangular LD block across ' + n + ' highlighted variants; top scale uses the same ' + formatWindowMb(figure) + ' genomic window as the Manhattan plot above, and curved connectors map SNP positions to evenly spaced matrix columns.';
  updateFigureLegend();
}
function updateFigureLegend() {
  const isDelta = heatMode.indexOf('delta') >= 0;
  byId('legendHeatTitle').textContent = isDelta ? 'heatmap delta r2' : 'heatmap r2';
  byId('legendGradient').classList.toggle('delta', isDelta);
  byId('legendHeatLabels').innerHTML = isDelta ? '<span>-0.7</span><span>-0.25</span><span>0</span><span>+0.25</span><span>+0.7</span>' : '<span>0</span><span>0.25</span><span>0.5</span><span>0.75</span><span>1</span>';
}
function drawLegend(ctx, x, y, isDelta) {
  const w = 250, h = 14;
  for (let i = 0; i < w; i++) {
    const t = i / (w - 1);
    const value = isDelta ? -0.7 + t * 1.4 : t;
    ctx.fillStyle = isDelta ? deltaColor(value) : r2Color(value);
    ctx.fillRect(x + i, y, 1, h);
  }
  ctx.strokeStyle = '#b8c3d1'; ctx.strokeRect(x, y, w, h);
  ctx.fillStyle = '#5f6f84'; ctx.font = '11px Segoe UI, Arial'; ctx.textAlign = 'center';
  ctx.fillText(isDelta ? '-0.7' : '0', x, y + 30);
  ctx.fillText(isDelta ? '0' : '0.5', x + w / 2, y + 30);
  ctx.fillText(isDelta ? '+0.7' : '1', x + w, y + 30);
}

function heatmapCellFromEvent(event) {
  const canvas = byId('heatmapCanvas');
  const geo = canvas._geo;
  if (!geo) return null;
  const rect = canvas.getBoundingClientRect();
  const scaleX = canvas.width / rect.width;
  const scaleY = canvas.height / rect.height;
  const x = (event.clientX - rect.left) * scaleX;
  const y = (event.clientY - rect.top) * scaleY;
  const n = geo.vars.length;
  if (y < geo.axisTop || y > geo.diamondTop + geo.heatHeight) return null;
  const colX = index => geo.axisLeft + (index + 0.5) * geo.cell;
  let best = null;
  let bestDistance = Infinity;
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      const centerX = (colX(i) + colX(j)) / 2;
      const centerY = geo.diamondTop + (colX(j) - colX(i)) / 2;
      const distance = Math.abs(x - centerX) + Math.abs(y - centerY);
      if (distance <= geo.cell / 2 + 1 && distance < bestDistance) {
        bestDistance = distance;
        best = {i, j};
      }
    }
  }
  return best;
}
function setupTooltip() {
  const canvas = byId('heatmapCanvas');
  canvas.addEventListener('mousemove', event => {
    const cell = heatmapCellFromEvent(event);
    if (!cell) { tooltip.style.display = 'none'; return; }
    const mat = getMatrixForMode();
    const a = DATA.selected[cell.i];
    const b = DATA.selected[cell.j];
    const value = mat[cell.i][cell.j];
    tooltip.innerHTML = '<strong>' + a.id + '</strong> x <strong>' + b.id + '</strong><br>' +
      fmtMb(a.pos) + ' / ' + fmtMb(b.pos) + '<br>' +
      (heatMode.indexOf('delta') >= 0 ? 'delta r2: ' : 'r2: ') + fmtNum(value, 3);
    tooltip.style.display = 'block';
    tooltip.style.left = (event.clientX + 14) + 'px';
    tooltip.style.top = (event.clientY + 14) + 'px';
  });
  canvas.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });
}

function fillDiscord(tableId, rows) {
  const body = byId(tableId).querySelector('tbody');
  body.innerHTML = rows.map(row => '<tr><td>' + row.a + '<br>' + row.b + '</td><td>' + fmtNum(row.ref, 3) + '</td><td>' + fmtNum(row.uzb, 3) + '</td><td>' + fmtNum(row.abs_delta, 3) + '</td></tr>').join('');
}
function fillVariants() {
  const body = byId('variantTable').querySelector('tbody');
  body.innerHTML = DATA.selected.map(v => '<tr><td>' + v.id + '</td><td>' + fmtBp(v.pos) + '</td><td>' + v.role + '</td><td>' + fmtP(v.p) + '</td><td>' + fmtNum(v.or, 3) + '</td><td>' + fmtNum(v.r2_uzb_lead, 3) + '</td><td>' + fmtNum(v.r2_eur_lead, 3) + '</td><td>' + fmtNum(v.r2_sas_lead, 3) + '</td><td>' + fmtP(v.saige_p) + '</td></tr>').join('');
}

document.querySelectorAll('button[data-mode]').forEach(button => {
  button.addEventListener('click', () => {
    heatMode = button.dataset.mode;
    document.querySelectorAll('button[data-mode]').forEach(b => b.classList.toggle('active', b === button));
    drawHeatmap();
  });
});

setHeader();
drawManhattan();
drawHeatmap();
setupTooltip();
fillDiscord('eurDiscord', DATA.discordance.eur_uzb);
fillDiscord('sasDiscord', DATA.discordance.sas_uzb);
fillVariants();
</script>
</body>
</html>
"""


ALIAS_TEMPLATE = """<!DOCTYPE html>
<html lang=\"en\">
<head>
<meta charset=\"utf-8\">
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
<meta http-equiv=\"refresh\" content=\"0; url={target}\">
<title>{title}</title>
<script>location.replace('{target}');</script>
</head>
<body><p>Redirecting to <a href=\"{target}\">{target}</a>.</p></body>
</html>
"""


def to_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if not value or value.upper() == "NA" or value.lower() == "nan":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().lstrip("#").split()
        rows: list[dict[str, str]] = []
        for line in handle:
            parts = line.split()
            if len(parts) != len(header):
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def read_gwas(config: LocusConfig) -> list[dict[str, Any]]:
    rows = []
    for row in read_table(DATA_DIR / f"chr{config.chrom}_region_gwas.tsv"):
        p_value = to_float(row.get("P"))
        position = int(row["POS"])
        if p_value is None:
            continue
        rows.append(
            {
                "chrom": config.chrom,
                "pos": position,
                "id": row["ID"],
                "ref": row.get("REF"),
                "alt": row.get("ALT"),
                "a1": row.get("A1"),
                "a1_freq": to_float(row.get("A1_FREQ")),
                "or": to_float(row.get("OR")),
                "p": p_value,
            }
        )
    return sorted(rows, key=lambda item: item["pos"])


def read_saige(config: LocusConfig) -> dict[str, dict[str, Any]]:
    path = DATA_DIR / f"chr{config.chrom}_region_saige.tsv"
    if not path.exists():
        return {}
    rows: dict[str, dict[str, Any]] = {}
    for row in read_table(path):
        marker = row.get("MarkerID")
        if not marker:
            continue
        p_value = to_float(row.get("p.value"))
        beta = to_float(row.get("BETA"))
        rows[marker] = {
            "id": marker,
            "pos": int(row["POS"]),
            "p": p_value,
            "beta": beta,
            "af_case": to_float(row.get("AF_case")),
            "af_ctrl": to_float(row.get("AF_ctrl")),
        }
    return rows


def read_vars(config: LocusConfig) -> list[str]:
    return (DATA_DIR / f"chr{config.chrom}_uzb_region.phased.vcor2.vars").read_text(encoding="utf-8").splitlines()


def parse_matrix_row(line: str) -> list[float | None]:
    return [to_float(value) for value in line.split()]


def read_matrix_subset(config: LocusConfig, row_indices: list[int], col_indices: list[int]) -> list[list[float | None]]:
    selected_rows = set(row_indices)
    output: list[list[float | None]] = []
    with (DATA_DIR / f"chr{config.chrom}_uzb_region.phased.vcor2").open("r", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if index not in selected_rows:
                continue
            values = parse_matrix_row(line)
            output.append([values[col] for col in col_indices])
    return output


def read_lead_r2(config: LocusConfig, lead_index: int) -> list[float | None]:
    with (DATA_DIR / f"chr{config.chrom}_uzb_region.phased.vcor2").open("r", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if index == lead_index:
                return parse_matrix_row(line)
    raise ValueError(f"Lead row {config.lead_id} not found in chr{config.chrom} matrix")


def select_variants(config: LocusConfig, gwas_rows: list[dict[str, Any]], vars_order: list[str], lead_r2: list[float | None]) -> list[dict[str, Any]]:
    index_by_id = {variant_id: index for index, variant_id in enumerate(vars_order)}
    for row in gwas_rows:
        row["matrix_index"] = index_by_id.get(row["id"])
        row["r2_to_lead"] = lead_r2[row["matrix_index"]] if row["matrix_index"] is not None else None
    selected: list[dict[str, Any]] = []
    for row in gwas_rows:
        if row["matrix_index"] is None:
            continue
        is_peak = row["p"] < PEAK_P
        is_extension = row["p"] >= PEAK_P and row["r2_to_lead"] is not None and row["r2_to_lead"] >= EXTENSION_R2
        is_lead = row["id"] == config.lead_id
        if is_peak or is_extension or is_lead:
            role = "lead" if is_lead else ("peak" if is_peak else "extension")
            item = dict(row)
            item["role"] = role
            selected.append(item)
    return sorted(selected, key=lambda item: item["pos"])


def ensembl_json(url: str) -> Any:
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "ALSU-analysis/1.0"},
    )
    with urllib.request.urlopen(request, timeout=35) as response:
        return json.load(response)


def fetch_reference_matrix(config: LocusConfig, ids: list[str], population: str) -> list[list[float | None]]:
    cache_path = DATA_DIR / f"chr{config.chrom}_1000g_{population}_ld.json"
    if cache_path.exists():
        cached = json.loads(cache_path.read_text(encoding="utf-8"))
    cached_matrix = cached.get("matrix")
    cached_pairs = matrix_stats(cached_matrix)["coverage"] if cached_matrix else 0
    if cached.get("ids") == ids and cached.get("population") == population and cached_pairs > 0:
            return cached["matrix"]

    index_by_id = {variant_id: index for index, variant_id in enumerate(ids)}
    matrix: list[list[float | None]] = [[None for _ in ids] for _ in ids]
    for i in range(len(ids)):
        matrix[i][i] = 1.0

    for query_index, variant_id in enumerate(ids, start=1):
        pop_name = f"1000GENOMES:phase_3:{population}"
        url = (
            f"https://rest.ensembl.org/ld/human/{variant_id}/{pop_name}"
            f"?window_size={REFERENCE_WINDOW_KB}&r2=0.0&d_prime=0.0&content-type=application/json"
        )
        try:
            response = ensembl_json(url)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as error:
            print(f"WARN chr{config.chrom} {population} {variant_id}: {error}")
            time.sleep(0.5)
            continue
        found = 0
        for entry in response:
            first = entry.get("variation1")
            second = entry.get("variation2")
            r2_value = to_float(str(entry.get("r2")))
            if first in index_by_id and second in index_by_id and r2_value is not None:
                i = index_by_id[first]
                j = index_by_id[second]
                value = round(r2_value, 6)
                matrix[i][j] = value
                matrix[j][i] = value
                found += 1
            elif second in index_by_id and r2_value is not None:
                i = index_by_id[variant_id]
                j = index_by_id[second]
                value = round(r2_value, 6)
                matrix[i][j] = value
                matrix[j][i] = value
                found += 1
        print(f"chr{config.chrom} {population} [{query_index}/{len(ids)}] {variant_id}: {found} selected pairs")
        time.sleep(0.15)

    cache_path.write_text(
        json.dumps({"population": population, "ids": ids, "matrix": matrix}, separators=(",", ":")),
        encoding="utf-8",
    )
    return matrix


def fetch_genes(config: LocusConfig) -> list[dict[str, Any]]:
    cache_path = DATA_DIR / f"chr{config.chrom}_ensembl_genes.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text(encoding="utf-8"))
    region = f"{config.chrom}:{config.display_start}-{config.display_end}"
    url = f"https://rest.ensembl.org/overlap/region/human/{region}?feature=gene;content-type=application/json"
    try:
        data = ensembl_json(url)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as error:
        print(f"WARN chr{config.chrom} gene fetch failed: {error}")
        data = []
    genes = []
    for item in data:
        biotype = item.get("biotype") or "gene"
        if biotype not in {"protein_coding", "lncRNA", "processed_pseudogene"}:
            continue
        start = int(item["start"])
        end = int(item["end"])
        if end < config.display_start or start > config.display_end:
            continue
        name = item.get("external_name")
        if not name:
            midpoint = (start + end) / 2e6
            name = f"{biotype.replace('_', ' ')} {midpoint:.2f} Mb"
        genes.append(
            {
                "name": name,
                "biotype": biotype,
                "start": start,
                "end": end,
                "strand": int(item.get("strand") or 1),
            }
        )
    genes = sorted(genes, key=lambda item: (item["start"], item["end"]))[:14]
    cache_path.write_text(json.dumps(genes, indent=2), encoding="utf-8")
    return genes


def round_matrix(matrix: list[list[float | None]]) -> list[list[float | None]]:
    return [[None if value is None else round(float(value), 4) for value in row] for row in matrix]


def matrix_stats(matrix: list[list[float | None]]) -> dict[str, Any]:
    values = []
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            value = matrix[i][j]
            if value is not None:
                values.append(value)
    return {
        "coverage": len(values),
        "mean_pair_r2": statistics.fmean(values) if values else None,
        "median_pair_r2": statistics.median(values) if values else None,
    }


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3:
        return None
    mean_x = statistics.fmean(xs)
    mean_y = statistics.fmean(ys)
    numerator = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    denom_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    denom_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if denom_x == 0 or denom_y == 0:
        return None
    return numerator / (denom_x * denom_y)


def compare_matrices(ref: list[list[float | None]], uzb: list[list[float | None]]) -> dict[str, Any]:
    ref_values = []
    uzb_values = []
    deltas = []
    for i in range(len(uzb)):
        for j in range(i + 1, len(uzb)):
            if ref[i][j] is None or uzb[i][j] is None:
                continue
            ref_values.append(float(ref[i][j]))
            uzb_values.append(float(uzb[i][j]))
            deltas.append(abs(float(ref[i][j]) - float(uzb[i][j])))
    return {
        "pairs": len(deltas),
        "pearson": pearson(ref_values, uzb_values),
        "mean_abs_delta": statistics.fmean(deltas) if deltas else None,
        "max_abs_delta": max(deltas) if deltas else None,
    }


def discordant_pairs(
    selected: list[dict[str, Any]],
    ref: list[list[float | None]],
    uzb: list[list[float | None]],
    limit: int = 12,
) -> list[dict[str, Any]]:
    rows = []
    for i in range(len(selected)):
        for j in range(i + 1, len(selected)):
            if ref[i][j] is None or uzb[i][j] is None:
                continue
            delta = float(ref[i][j]) - float(uzb[i][j])
            rows.append(
                {
                    "a": selected[i]["id"],
                    "b": selected[j]["id"],
                    "ref": round(float(ref[i][j]), 4),
                    "uzb": round(float(uzb[i][j]), 4),
                    "delta": round(delta, 4),
                    "abs_delta": round(abs(delta), 4),
                }
            )
    return sorted(rows, key=lambda item: item["abs_delta"], reverse=True)[:limit]


def build_payload(config: LocusConfig) -> dict[str, Any]:
    gwas_rows = read_gwas(config)
    vars_order = read_vars(config)
    index_by_id = {variant_id: index for index, variant_id in enumerate(vars_order)}
    if config.lead_id not in index_by_id:
        raise ValueError(f"{config.lead_id} is missing from chr{config.chrom} matrix vars")
    lead_r2 = read_lead_r2(config, index_by_id[config.lead_id])
    selected = select_variants(config, gwas_rows, vars_order, lead_r2)
    selected_ids = [row["id"] for row in selected]
    selected_indices = [index_by_id[variant_id] for variant_id in selected_ids]
    uzb_matrix = read_matrix_subset(config, selected_indices, selected_indices)
    eur_matrix = fetch_reference_matrix(config, selected_ids, "EUR")
    sas_matrix = fetch_reference_matrix(config, selected_ids, "SAS")
    genes = fetch_genes(config)
    saige_by_id = read_saige(config)
    peak_count = sum(1 for row in selected if row["role"] in {"lead", "peak"})
    extension_count = sum(1 for row in selected if row["role"] == "extension")
    lead_row = next(row for row in selected if row["id"] == config.lead_id)
    peak_positions = [row["pos"] for row in selected if row["role"] in {"lead", "peak"}]
    lead_index = selected_ids.index(config.lead_id)

    for index, row in enumerate(selected):
        row["r2_uzb_lead"] = uzb_matrix[index][lead_index]
        row["r2_eur_lead"] = eur_matrix[index][lead_index]
        row["r2_sas_lead"] = sas_matrix[index][lead_index]
        saige = saige_by_id.get(row["id"])
        row["saige_p"] = saige["p"] if saige else None

    gwas_payload = []
    selected_id_set = set(selected_ids)
    for row in gwas_rows:
        gwas_payload.append(
            {
                "id": row["id"],
                "pos": row["pos"],
                "p": row["p"],
                "or": row["or"],
                "r2_to_lead": row.get("r2_to_lead"),
                "selected": row["id"] in selected_id_set,
            }
        )

    best_saige = None
    valid_saige = [row for row in saige_by_id.values() if row.get("p") is not None]
    if valid_saige:
        best_saige = min(valid_saige, key=lambda item: item["p"])

    payload = {
        "chrom": config.chrom,
        "header": f"{config.label} - RPL - LD Comparison",
        "subtitle": (
            f"{len(selected)} SNPs - {peak_count} peak + {extension_count} extension - "
            f"{config.subtitle} - UZB LD from 1090 imputed samples"
        ),
        "context": config.context,
        "badges": [
            f"lead: {config.lead_id} - p={lead_row['p']:.3g} - OR={lead_row['or']:.3g}",
            f"No GWS - best p={lead_row['p']:.3g}",
            f"SAIGE lead: {best_saige['id']} - p={best_saige['p']:.3g}" if best_saige else "SAIGE: not available",
            "EUR/SAS: 1000G phase 3",
        ],
        "lead": {"id": config.lead_id, "pos": config.lead_pos, "p": lead_row["p"], "or": lead_row["or"]},
        "window": [config.start, config.end],
        "figure_window": [config.start, config.end],
        "display_window": [config.display_start, config.display_end],
        "peak_span": [min(peak_positions), max(peak_positions)] if peak_positions else [config.lead_pos, config.lead_pos],
        "counts": {"peak": peak_count, "extension": extension_count, "regional_gwas": len(gwas_rows)},
        "selected": [
            {
                "id": row["id"],
                "pos": row["pos"],
                "role": row["role"],
                "p": row["p"],
                "or": row["or"],
                "r2_uzb_lead": row["r2_uzb_lead"],
                "r2_eur_lead": row["r2_eur_lead"],
                "r2_sas_lead": row["r2_sas_lead"],
                "saige_p": row["saige_p"],
            }
            for row in selected
        ],
        "gwas": gwas_payload,
        "genes": genes,
        "matrices": {
            "uzb": round_matrix(uzb_matrix),
            "eur": round_matrix(eur_matrix),
            "sas": round_matrix(sas_matrix),
        },
        "stats": {
            "uzb": matrix_stats(uzb_matrix),
            "eur": matrix_stats(eur_matrix),
            "sas": matrix_stats(sas_matrix),
        },
        "comparisons": {
            "eur_uzb": compare_matrices(eur_matrix, uzb_matrix),
            "sas_uzb": compare_matrices(sas_matrix, uzb_matrix),
        },
        "discordance": {
            "eur_uzb": discordant_pairs(selected, eur_matrix, uzb_matrix),
            "sas_uzb": discordant_pairs(selected, sas_matrix, uzb_matrix),
        },
        "saige": {"best": best_saige},
        "functional_page": config.functional_page,
        "functional_label": f"chr{config.chrom} functional analysis",
        "method": {
            "peak_p": PEAK_P,
            "extension_r2": EXTENSION_R2,
            "reference_window_kb": REFERENCE_WINDOW_KB,
        },
    }
    return payload


def write_page(config: LocusConfig, payload: dict[str, Any]) -> None:
    json_payload = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    html = HTML_TEMPLATE.replace("__TITLE__", payload["header"]).replace("__PAYLOAD__", json_payload)
    canonical = OUT_DIR / config.canonical_page
    canonical.write_text(html, encoding="utf-8")
    alias = ALIAS_TEMPLATE.format(target=config.canonical_page, title=payload["header"])
    (OUT_DIR / config.legacy_page).write_text(alias, encoding="utf-8")
    (DATA_DIR / f"chr{config.chrom}_dense_payload.json").write_text(json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")
    print(f"Wrote {canonical.relative_to(ROOT)} and {config.legacy_page}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for config in LOCI:
        print(f"\n=== chr{config.chrom} {config.lead_id} ===")
        payload = build_payload(config)
        write_page(config, payload)


if __name__ == "__main__":
    main()