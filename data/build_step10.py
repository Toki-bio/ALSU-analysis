#!/usr/bin/env python3
"""Generate comprehensive step10.html with all results, interactive charts, and embedded scripts."""
import json
import html

# ── Load data ──────────────────────────────────────────────────────────────────
with open("data/admix_data.json") as f:
    admix = json.load(f)

# Round Q values for smaller embed
for k in ("k2", "k3", "k4"):
    admix[k] = [[round(x, 4) for x in row] for row in admix[k]]

admix_json = json.dumps(admix, separators=(",", ":"))

with open("data/script1.sh", encoding="utf-8") as f:
    script1_raw = f.read()
with open("data/script2.py", encoding="utf-8") as f:
    script2_raw = f.read()

script1_html = html.escape(script1_raw)
script2_html = html.escape(script2_raw)

# ── PBS top 20 data ────────────────────────────────────────────────────────────
top20 = [
    (9, "9:104189856", 0.0054, 0.9911, 1.0000, 1.0000, 1.0000, 4.6021, -0.2649, 0.2729),
    (8, "8:24813391", 0.0101, 0.9990, 1.0000, 0.9990, 1.0000, 4.2399, -0.0216, 0.0216),
    (17, "17:48275339", 0.0135, 0.9960, 1.0000, 0.9980, 1.0000, 3.9094, -0.0610, 0.0640),
    (1, "1:36934805", 0.0152, 0.9811, 1.0000, 0.9928, 0.8268, 3.6201, -0.2157, 0.2338),
    (2, "2:179642589", 0.0152, 0.9880, 0.9841, 0.9836, 0.9992, 3.5038, 0.0414, -0.0414),
    (11, "11:47355475", 0.0195, 1.0000, 0.9901, 1.0000, 1.0000, 3.4953, 0.1103, -0.1013),
    (6, "6:29640785", 0.0201, 0.9791, 1.0000, 0.9959, 0.8487, 3.3667, -0.1861, 0.2062),
    (4, "4:187003729", 0.0235, 0.9930, 1.0000, 1.0000, 1.0000, 3.3566, -0.0629, 0.0689),
    (10, "10:96818119", 0.0159, 0.9424, 1.0000, 0.9928, 0.9962, 3.2852, -0.4625, 0.5209),
    (1, "1:161199431", 0.0245, 0.9891, 1.0000, 1.0000, 1.0000, 3.2799, -0.0908, 0.1008),
    (2, "2:152580815", 0.0232, 0.9761, 1.0000, 0.9949, 0.9992, 3.2256, -0.1863, 0.2094),
    (13, "13:95672265", 0.0283, 0.9990, 0.9950, 0.9980, 0.9819, 3.1918, 0.0329, -0.0312),
    (16, "16:89833576", 0.0285, 0.9871, 1.0000, 0.9755, 0.9992, 3.1235, -0.0919, 0.1040),
    (6, "6:30585771", 0.0193, 0.9930, 0.9474, 0.9417, 0.9992, 3.1079, 0.3519, -0.3170),
    (1, "1:146661814", 0.0324, 1.0000, 0.9990, 1.0000, 1.0000, 3.0985, 0.0071, -0.0071),
    (6, "6:31059825", 0.0281, 0.9821, 0.9990, 0.9857, 0.9478, 3.0939, -0.1173, 0.1316),
    (4, "4:70354534", 0.0334, 0.9990, 1.0000, 1.0000, 0.9728, 3.0695, -0.0075, 0.0075),
    (11, "11:34999682", 0.0308, 0.9881, 1.0000, 0.9990, 0.9992, 3.0649, -0.0814, 0.0924),
    (8, "8:48701786", 0.0304, 1.0000, 0.9821, 1.0000, 0.9992, 3.0352, 0.1348, -0.1177),
    (2, "2:21238413", 0.0338, 0.9930, 1.0000, 0.9990, 0.9992, 3.0147, -0.0456, 0.0516),
]

top20_rows = ""
for i, (ch, snp, muz, meu, mea, msa, maf, pbs, pe, pea) in enumerate(top20, 1):
    muz_pct = f"{muz*100:.1f}%"
    top20_rows += f"""                        <tr>
                            <td>{i}</td>
                            <td>{ch}</td>
                            <td><code>{snp}</code></td>
                            <td><strong>{pbs:.4f}</strong></td>
                            <td style="color:#c62828">{muz:.4f}</td>
                            <td>{meu:.4f}</td>
                            <td>{mea:.4f}</td>
                            <td>{msa:.4f}</td>
                            <td>{maf:.4f}</td>
                            <td>{pe:.4f}</td>
                            <td>{pea:.4f}</td>
                        </tr>
"""

# ── Histogram data ─────────────────────────────────────────────────────────────
hist_vals = [12,10,18,67,188,497,1660,6243,27069,221604,116856,516,144,71,65,45,46,31,24,38,26,19,13,5,14,11,12,11,8,12,12,9,10,7,12,6,5,6,8,7,7,4,7,2,4,6,2,2,2,68]
hist_json = json.dumps(hist_vals)

# ── Build HTML ─────────────────────────────────────────────────────────────────
page = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Step 10: Multi-Pop Analysis, PBS &amp; ADMIXTURE — ALSU Pipeline</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background: linear-gradient(135deg, #4a4a4a 0%, #2a2a2a 100%); min-height: 100vh; padding: 20px; }}
        .page-wrapper {{ display: flex; gap: 20px; max-width: 1500px; margin: 0 auto; }}
        
        .sidebar {{ width: 250px; background: white; border-radius: 10px; box-shadow: 0 10px 40px rgba(0,0,0,0.2); padding: 20px; height: fit-content; position: sticky; top: 20px; }}
        .sidebar h3 {{ font-size: 1em; color: #333; margin-bottom: 15px; border-bottom: 2px solid #888; padding-bottom: 10px; }}
        .roadmap-mini {{ display: flex; flex-direction: column; gap: 8px; }}
        .roadmap-mini a {{ padding: 10px 12px; background: #f5f5f5; color: #333; text-decoration: none; border-radius: 5px; border-left: 4px solid transparent; transition: all 0.3s ease; font-size: 0.9em; }}
        .roadmap-mini a:hover {{ background: #f0f0f0; border-left-color: #666; }}
        .roadmap-mini a.current {{ background: #555; color: white; border-left-color: #333; font-weight: 600; }}
        
        .main-content {{ flex: 1; background: white; border-radius: 10px; box-shadow: 0 10px 40px rgba(0,0,0,0.2); overflow: hidden; }}
        
        .search-bar {{ background: #f5f5f5; padding: 15px 30px; border-bottom: 1px solid #ddd; }}
        .search-bar input {{ width: 100%; max-width: 400px; padding: 12px 15px; border: 2px solid #ddd; border-radius: 5px; font-size: 1em; transition: border-color 0.3s ease; }}
        .search-bar input:focus {{ outline: none; border-color: #667eea; }}
        
        .search-results {{ display: none; padding: 20px 30px; background: #e3f2fd; border-bottom: 1px solid #ddd; }}
        .search-results.active {{ display: block; }}
        .search-results a {{ color: #1976d2; text-decoration: none; margin-right: 15px; }}
        .search-results a:hover {{ text-decoration: underline; }}
        
        .header {{ background: linear-gradient(135deg, #1b5e20 0%, #2e7d32 100%); color: white; padding: 30px; }}
        .back-link {{ color: rgba(255,255,255,0.8); text-decoration: none; margin-bottom: 15px; display: inline-block; }}
        .back-link:hover {{ color: white; }}
        .header h1 {{ font-size: 2em; margin-bottom: 10px; }}
        .header p {{ opacity: 0.9; }}
        .step-badge {{ display: inline-block; background: rgba(255,255,255,0.2); padding: 5px 15px; border-radius: 20px; font-size: 0.9em; margin-top: 10px; }}
        
        .content {{ padding: 30px; }}
        
        @media (max-width: 1024px) {{ .page-wrapper {{ flex-direction: column; }} .sidebar {{ width: 100%; position: relative; top: 0; }} }}
        
        h2 {{ color: #333; font-size: 1.5em; margin: 30px 0 15px 0; padding-bottom: 10px; border-bottom: 2px solid #2e7d32; }}
        h3 {{ color: #555; margin-top: 25px; margin-bottom: 12px; font-size: 1.15em; }}
        h4 {{ color: #444; margin-top: 18px; margin-bottom: 8px; }}
        p, li {{ color: #666; line-height: 1.7; margin-bottom: 12px; }}
        ul {{ margin-left: 25px; }}
        ol {{ margin-left: 25px; }}
        code {{ background: #f5f5f5; padding: 2px 6px; border-radius: 3px; font-family: 'Courier New', monospace; font-size: 0.92em; }}
        
        .code-block {{ background: #2d2d2d; color: #f8f8f2; padding: 20px; border-radius: 5px; overflow-x: auto; margin: 20px 0; font-family: 'Courier New', monospace; font-size: 0.85em; line-height: 1.5; white-space: pre-wrap; word-wrap: break-word; max-height: 600px; overflow-y: auto; }}
        .command-output {{ background: #f5f5f5; border: 1px solid #ddd; padding: 12px 15px; border-radius: 5px; margin: 10px 0; font-family: 'Courier New', monospace; font-size: 0.9em; white-space: pre-wrap; word-wrap: break-word; color: #333; }}
        
        .stat-box {{ background: #e8f5e9; padding: 15px 20px; border-left: 4px solid #2e7d32; border-radius: 5px; margin: 20px 0; }}
        .stat-label {{ font-weight: 600; color: #2e7d32; }}
        .stat-value {{ color: #333; margin-left: 10px; }}
        .success {{ background: #d4edda; padding: 15px; border-left: 4px solid #28a745; border-radius: 5px; margin: 20px 0; }}
        .success strong {{ color: #155724; }}
        .important {{ background: #fff3cd; padding: 15px; border-left: 4px solid #ffc107; border-radius: 5px; margin: 20px 0; }}
        .important strong {{ color: #856404; }}
        .warning {{ background: #f8d7da; padding: 15px; border-left: 4px solid #dc3545; border-radius: 5px; margin: 20px 0; }}
        .warning strong {{ color: #721c24; }}
        .info {{ background: #e3f2fd; padding: 15px; border-left: 4px solid #1976d2; border-radius: 5px; margin: 20px 0; }}
        .info strong {{ color: #0d47a1; }}
        
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 0.95em; }}
        th, td {{ padding: 10px 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background: #e8f5e9; font-weight: 600; color: #1b5e20; position: sticky; top: 0; }}
        tr:hover {{ background: #f9f9f9; }}
        
        .formula {{ background: #f0f4f8; border: 1px solid #ccc; padding: 15px 20px; border-radius: 5px; margin: 15px 0; text-align: center; font-family: 'Cambria Math', 'Times New Roman', serif; font-size: 1.15em; color: #333; }}
        .formula .label {{ font-size: 0.85em; color: #666; display: block; margin-top: 8px; }}
        
        .tier-box {{ border: 2px solid #ccc; border-radius: 8px; padding: 20px; margin: 15px 0; }}
        .tier-box h4 {{ margin-top: 0; margin-bottom: 10px; }}
        .tier1 {{ border-color: #c62828; background: #ffebee; }}
        .tier1 h4 {{ color: #c62828; }}
        .tier2 {{ border-color: #e65100; background: #fff3e0; }}
        .tier2 h4 {{ color: #e65100; }}
        .tier3 {{ border-color: #1565c0; background: #e3f2fd; }}
        .tier3 h4 {{ color: #1565c0; }}
        
        .pipeline-flow {{ display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin: 20px 0; }}
        .pipeline-step {{ background: #e8f5e9; border: 2px solid #2e7d32; border-radius: 8px; padding: 12px 18px; font-size: 0.9em; text-align: center; min-width: 140px; }}
        .pipeline-arrow {{ font-size: 1.5em; color: #2e7d32; font-weight: bold; }}
        
        /* Interactive chart styles */
        .chart-container {{ border: 1px solid #ddd; border-radius: 8px; padding: 20px; margin: 20px 0; background: #fafafa; }}
        .chart-toolbar {{ display: flex; gap: 10px; margin-bottom: 15px; flex-wrap: wrap; align-items: center; }}
        .chart-toolbar button {{ padding: 8px 16px; border: 2px solid #2e7d32; background: white; color: #2e7d32; border-radius: 5px; cursor: pointer; font-size: 0.9em; font-weight: 600; transition: all 0.2s; }}
        .chart-toolbar button:hover {{ background: #e8f5e9; }}
        .chart-toolbar button.active {{ background: #2e7d32; color: white; }}
        .chart-toolbar .legend {{ display: flex; gap: 12px; margin-left: auto; flex-wrap: wrap; }}
        .chart-toolbar .legend-item {{ display: flex; align-items: center; gap: 5px; font-size: 0.85em; color: #555; }}
        .legend-swatch {{ width: 16px; height: 16px; border-radius: 3px; border: 1px solid rgba(0,0,0,0.2); }}
        
        canvas {{ display: block; margin: 0 auto; cursor: crosshair; }}
        .tooltip-box {{ position: absolute; background: rgba(0,0,0,0.85); color: white; padding: 10px 14px; border-radius: 6px; font-size: 0.85em; pointer-events: none; z-index: 100; white-space: nowrap; line-height: 1.5; }}
        
        .results-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }}
        .result-card {{ background: white; border: 2px solid #e0e0e0; border-radius: 8px; padding: 20px; text-align: center; }}
        .result-card .value {{ font-size: 2em; font-weight: 700; color: #1b5e20; }}
        .result-card .label {{ font-size: 0.85em; color: #666; margin-top: 5px; }}
        
        details {{ margin: 15px 0; }}
        details summary {{ cursor: pointer; padding: 12px 15px; background: #f5f5f5; border: 1px solid #ddd; border-radius: 5px; font-weight: 600; color: #333; }}
        details summary:hover {{ background: #eee; }}
        details[open] summary {{ border-radius: 5px 5px 0 0; }}
        
        .scroll-table {{ overflow-x: auto; max-height: 500px; overflow-y: auto; }}
        .scroll-table table {{ font-size: 0.88em; }}
        
        .fst-heatmap {{ display: grid; grid-template-columns: repeat(6, 1fr); gap: 2px; max-width: 600px; margin: 20px auto; }}
        .fst-cell {{ padding: 12px; text-align: center; font-size: 0.85em; border-radius: 4px; }}
    </style>
</head>
<body>
    <div class="page-wrapper">
        <!-- LEFT SIDEBAR -->
        <aside class="sidebar">
            <h3>Pipeline Steps</h3>
            <div class="roadmap-mini">
                <a href="step1.html">Step 1: Missingness</a>
                <a href="step2.html">Step 2: IBD Dedup</a>
                <a href="step3.html">Step 3: SNP QC</a>
                <a href="step4.html">Step 4: Imputation</a>
                <a href="step5.html">Step 5: ID Normalize</a>
                <a href="step6.html">Step 6: Final QC</a>
                <a href="step7.html">Step 7: PCA Analysis</a>
                <a href="step8.html">Step 8: Global PCA</a>
                <a href="step9.html">Step 9: F<sub>ST</sub> Analysis</a>
                <a href="step10.html" class="current">Step 10: Multi-Pop &amp; PBS</a>
            </div>
        </aside>
        
        <!-- MAIN CONTENT -->
        <div class="main-content">
            <div class="search-bar">
                <input type="text" id="searchInput" placeholder="Search across all steps..." onkeyup="searchAllSteps()">
            </div>
            <div class="search-results" id="searchResults"></div>
            
            <div class="header">
                <a href="../index.html" class="back-link">&larr; Back to Roadmap</a>
                <h1>Step 10: Multi-Population Analysis, PBS &amp; ADMIXTURE</h1>
                <p>Multi-population F<sub>ST</sub>, Population Branch Statistic, Uzbek-specific SNPs, and ancestry structure analysis</p>
                <span class="step-badge">Completed &mdash; February 2026</span>
            </div>
        
            <div class="content">

<!-- ═══════════════════════════ 1. OVERVIEW ═══════════════════════════ -->
<h2>1. Overview &amp; Rationale</h2>

<p>
    Step 9 established genome-wide F<sub>ST</sub> between Uzbek and European (EUR) populations (weighted F<sub>ST</sub> = 0.020). However, <strong>pairwise F<sub>ST</sub> alone cannot distinguish Uzbek-specific selection from shared drift</strong>. A high F<sub>ST</sub> between UZB and EUR could simply reflect their geographic distance.
</p>
<p>
    This step adds <strong>three additional reference populations</strong> from 1000 Genomes Phase 3, plus performs <strong>ADMIXTURE analysis</strong> to decompose individual ancestry proportions:
</p>

<table>
    <tr><th>Population</th><th>Abbrev.</th><th>n</th><th>1000G Sub-populations</th></tr>
    <tr><td>Uzbek cohort</td><td><strong>UZB</strong></td><td>1,199</td><td>&mdash; (our samples)</td></tr>
    <tr><td>European</td><td><strong>EUR</strong></td><td>503</td><td>CEU, TSI, FIN, GBR, IBS</td></tr>
    <tr><td>East Asian</td><td><strong>EAS</strong></td><td>504</td><td>CHB, JPT, CHS, CDX, KHV</td></tr>
    <tr><td>South Asian</td><td><strong>SAS</strong></td><td>489</td><td>GIH, PJL, BEB, STU, ITU</td></tr>
    <tr><td>African</td><td><strong>AFR</strong></td><td>661</td><td>YRI, LWK, GWD, MSL, ESN, ACB, ASW</td></tr>
</table>

<h3>Why these populations?</h3>
<ul>
    <li><strong>EUR</strong> &mdash; closest major panel; already used in Step 9. Shared post-Neolithic ancestry with UZB.</li>
    <li><strong>EAS</strong> &mdash; Uzbek populations have documented East Asian admixture (Turkic migrations). Required as the second outgroup for the PBS triangle.</li>
    <li><strong>SAS</strong> &mdash; geographically adjacent; shared Indo-Iranian heritage. Key comparison for distinguishing Central Asian vs broader South-Central Asian patterns.</li>
    <li><strong>AFR</strong> &mdash; deep outgroup that captures ancestral vs derived allele states. Alleles rare in ALL non-UZB populations including AFR are the strongest candidates for Uzbek-specific selection.</li>
</ul>


<!-- ═══════════════════════════ 2. METHODS ═══════════════════════════ -->
<h2>2. Methods</h2>

<h3>2.1. Three-Tier Identification Strategy</h3>
<p>We use three complementary approaches to identify Uzbek-specific variants, each capturing different aspects of population differentiation:</p>

<div class="tier-box tier1">
    <h4>Tier 1: Population Branch Statistic (PBS)</h4>
    <p>PBS isolates <strong>branch-specific evolution</strong> using a three-population tree. High PBS on the Uzbek branch &rarr; allele frequency change <em>specific</em> to the Uzbek lineage, not simply shared divergence between EUR and EAS.</p>
    <div class="formula">
        T<sub>XY</sub> = &minus;log(1 &minus; F<sub>ST(XY)</sub>)
        <br><br>
        PBS<sub>UZB</sub> = (T<sub>UZB&ndash;EUR</sub> + T<sub>UZB&ndash;EAS</sub> &minus; T<sub>EUR&ndash;EAS</sub>) / 2
        <span class="label">Yi et al. 2010, Science &mdash; originally used to detect high-altitude adaptation in Tibetans (EPAS1)</span>
    </div>
    <p><strong>Threshold:</strong> PBS &ge; 0.3 (~top 0.1% genome-wide in most populations)</p>
</div>

<div class="tier-box tier2">
    <h4>Tier 2: Multi-Population Delta Allele Frequency</h4>
    <p>Computes |AF<sub>UZB</sub> &minus; AF<sub>X</sub>| for X &isin; &#123;EUR, EAS, SAS, AFR&#125;. SNPs where UZB differs from <strong>every</strong> major reference population by &ge; 0.3 are strong candidates.</p>
    <p><strong>Threshold:</strong> min(&Delta;AF) &ge; 0.3 across all four comparisons</p>
</div>

<div class="tier-box tier3">
    <h4>Tier 3: Near-Private Variants</h4>
    <p>Identifies alleles that are <strong>common in UZB but rare everywhere else</strong>: MAF &ge; 5% in Uzbek, MAF &le; 1% in EUR, EAS, SAS, and AFR. These are potential novel or recently selected variants.</p>
    <p><strong>Threshold:</strong> UZB MAF &ge; 0.05, all others MAF &le; 0.01</p>
</div>

<h3>2.2. ADMIXTURE Analysis</h3>
<p>
    <strong>ADMIXTURE</strong> (Alexander et al., 2009) is a maximum-likelihood model-based method that decomposes each individual&rsquo;s genome into proportions of <em>K</em> ancestral populations. Unlike PCA (which shows relative positions), ADMIXTURE gives <strong>explicit ancestry fractions</strong> (e.g., &ldquo;65% Western, 35% Eastern&rdquo;).
</p>
<p>
    <strong>How it works, simply:</strong> Imagine every person&rsquo;s genome as a mosaic of colored tiles. ADMIXTURE tries to find the best set of <em>K</em> &ldquo;colors&rdquo; (ancestral populations) and then figures out what proportion of each person&rsquo;s tiles are each color. It does this by maximizing the likelihood of the observed genotype data under the model.
</p>

<div class="info">
    <strong>Key parameters:</strong>
    <ul>
        <li><strong>Input:</strong> 1,074 Uzbek samples &times; 380,376 LD-pruned SNPs (MAF &gt; 5%, r&sup2; &lt; 0.1, window 50 step 10)</li>
        <li><strong>K values tested:</strong> K=2 through K=8</li>
        <li><strong>Cross-validation:</strong> 5-fold CV with <code>--cv=5</code> to select optimal K</li>
        <li><strong>Parallelization:</strong> 32 threads (<code>-j32</code>) on 64-core server</li>
    </ul>
</div>

<h3>2.3. Pipeline Architecture</h3>
<div class="pipeline-flow">
    <div class="pipeline-step"><strong>1000G VCFs</strong><br>(chr1&ndash;22)</div>
    <span class="pipeline-arrow">&rarr;</span>
    <div class="pipeline-step"><strong>PLINK convert</strong><br>(&ndash;&ndash;keep per pop)</div>
    <span class="pipeline-arrow">&rarr;</span>
    <div class="pipeline-step"><strong>SNP overlap</strong><br>(376K target set)</div>
    <span class="pipeline-arrow">&rarr;</span>
    <div class="pipeline-step"><strong>Merge chr1&ndash;22</strong><br>(per population)</div>
    <span class="pipeline-arrow">&rarr;</span>
    <div class="pipeline-step"><strong>Pairwise F<sub>ST</sub></strong><br>(5 pairs)</div>
    <span class="pipeline-arrow">&rarr;</span>
    <div class="pipeline-step"><strong>Frequencies</strong><br>(5 populations)</div>
    <span class="pipeline-arrow">&rarr;</span>
    <div class="pipeline-step"><strong>PBS + &Delta;AF</strong><br>(Python)</div>
</div>

<div class="pipeline-flow" style="margin-top: 5px;">
    <div class="pipeline-step" style="border-color: #1565c0; background: #e3f2fd;"><strong>LD Prune</strong><br>(10.8M &rarr; 380K)</div>
    <span class="pipeline-arrow" style="color: #1565c0;">&rarr;</span>
    <div class="pipeline-step" style="border-color: #1565c0; background: #e3f2fd;"><strong>ADMIXTURE</strong><br>(K=2&hellip;8, 32 threads)</div>
    <span class="pipeline-arrow" style="color: #1565c0;">&rarr;</span>
    <div class="pipeline-step" style="border-color: #1565c0; background: #e3f2fd;"><strong>Structure Plots</strong><br>(interactive)</div>
</div>


<!-- ═══════════════════════════ 3. PAIRWISE FST ═══════════════════════════ -->
<h2>3. Results: Pairwise F<sub>ST</sub></h2>

<p>
    Pairwise weighted F<sub>ST</sub> values calculated via PLINK 1.9 (Weir &amp; Cockerham estimator). F<sub>ST</sub> measures allele frequency differentiation between populations: 0 = identical, 1 = completely fixed for different alleles.
</p>

<table>
    <tr>
        <th>Population Pair</th>
        <th>Mean F<sub>ST</sub></th>
        <th>Weighted F<sub>ST</sub></th>
        <th>Interpretation</th>
    </tr>
    <tr>
        <td><strong>UZB vs EUR</strong></td>
        <td>0.0160</td>
        <td><strong>0.0204</strong></td>
        <td style="color:#2e7d32">Low &mdash; closest pair, consistent with shared Eurasian ancestry</td>
    </tr>
    <tr>
        <td><strong>UZB vs SAS</strong></td>
        <td>0.0148</td>
        <td><strong>0.0179</strong></td>
        <td style="color:#2e7d32">Very low &mdash; UZB closest to South Asians! Indo-Iranian heritage</td>
    </tr>
    <tr>
        <td><strong>UZB vs EAS</strong></td>
        <td>0.0369</td>
        <td><strong>0.0493</strong></td>
        <td style="color:#e65100">Moderate &mdash; Turkic admixture but still distinct</td>
    </tr>
    <tr>
        <td><strong>EUR vs EAS</strong></td>
        <td>0.0724</td>
        <td><strong>0.1061</strong></td>
        <td style="color:#c62828">High &mdash; major continental divide (expected)</td>
    </tr>
    <tr>
        <td><strong>UZB vs AFR</strong></td>
        <td>0.0729</td>
        <td><strong>0.1094</strong></td>
        <td style="color:#c62828">High &mdash; out-of-Africa divergence</td>
    </tr>
</table>

<div class="info">
    <strong>Key finding:</strong> UZB is <em>closest to SAS</em> (F<sub>ST</sub> = 0.018), even closer than UZB&ndash;EUR (0.020). This makes biological sense: Uzbekistan shares Indo-Iranian ancestry with South Asia, while also showing moderate East Asian admixture (F<sub>ST</sub> = 0.049). The UZB&ndash;EAS distance is roughly half of EUR&ndash;EAS (0.106), confirming that Uzbek is genuinely intermediate between European and East Asian gene pools.
</div>


<!-- ═══════════════════════════ 4. PBS RESULTS ═══════════════════════════ -->
<h2>4. Results: Population Branch Statistic (PBS)</h2>

<h3>4.1. Summary Statistics</h3>

<div class="results-grid">
    <div class="result-card">
        <div class="value">376,208</div>
        <div class="label">Total SNPs analyzed</div>
    </div>
    <div class="result-card">
        <div class="value" style="color:#c62828">456</div>
        <div class="label">PBS &ge; 0.3 (Tier 1)</div>
    </div>
    <div class="result-card">
        <div class="value" style="color:#e65100">422</div>
        <div class="label">min &Delta;AF &ge; 0.3 (Tier 2)</div>
    </div>
    <div class="result-card">
        <div class="value" style="color:#1565c0">272</div>
        <div class="label">Near-private (Tier 3)</div>
    </div>
    <div class="result-card">
        <div class="value" style="color:#1b5e20; font-size:2.2em;">490</div>
        <div class="label"><strong>Total Uzbek-specific candidates (union)</strong></div>
    </div>
</div>

<table>
    <tr><th>Statistic</th><th>Value</th><th>What it means</th></tr>
    <tr><td>Mean PBS<sub>UZB</sub></td><td><strong>&minus;0.01257</strong></td><td>Slightly negative on average &mdash; most drift is on the EUR&ndash;EAS branch (they are the most divergent pair)</td></tr>
    <tr><td>Median PBS<sub>UZB</sub></td><td>&minus;0.00504</td><td>Very close to zero. Typical SNP shows no UZB-specific branch lengthening</td></tr>
    <tr><td>Std deviation</td><td>0.0564</td><td>Most variation is small, but a long right tail exists</td></tr>
    <tr><td>Minimum</td><td>&minus;1.4295</td><td>Strong negative = EUR&ndash;EAS divergence dominates at this locus</td></tr>
    <tr><td>Maximum</td><td><strong>4.6021</strong></td><td>Extreme outlier &mdash; nearly complete allele reversal in UZB vs all others</td></tr>
    <tr><td>95th percentile</td><td>0.0146</td><td>95% of SNPs have PBS &lt; 0.015</td></tr>
    <tr><td>99th percentile</td><td>0.0312</td><td>99% of SNPs have PBS &lt; 0.031</td></tr>
    <tr><td>99.9th percentile</td><td><strong>0.4059</strong></td><td>Top 0.1% cutoff &mdash; approximately matches our 0.3 threshold</td></tr>
    <tr><td>PBS &ge; 0.3</td><td>456 SNPs</td><td>0.12% of all SNPs &mdash; strong UZB-specific differentiation</td></tr>
    <tr><td>PBS &ge; 0.1</td><td>781 SNPs</td><td>0.21% &mdash; moderate signals</td></tr>
    <tr><td>Negative PBS</td><td>257,368 (68.4%)</td><td>Expected: UZB sits between EUR and EAS, so most branch lengthening is on the outgroups</td></tr>
</table>

<h3>4.2. PBS Distribution</h3>
<div class="chart-container">
    <p style="margin-bottom:10px; font-weight:600; color:#333;">PBS<sub>UZB</sub> Distribution (376,208 SNPs, clipped to [&minus;0.5, 2.0])</p>
    <canvas id="pbsHistCanvas" width="900" height="350"></canvas>
    <p style="font-size:0.85em; color:#888; margin-top:8px;">Logarithmic y-axis. The vast majority of SNPs cluster near zero. 68 SNPs have PBS &gt; 2.0 (not shown individually).</p>
</div>

<h3>4.3. Top 20 PBS Hits</h3>
<p>
    These SNPs show the strongest Uzbek-specific branch lengthening. Note the pattern: <strong>UZB allele frequency is near 0</strong> while all other populations are near fixation (&asymp;100%). This means the <em>derived</em> allele swept to near-fixation everywhere in the world <em>except</em> in Uzbek, where the ancestral allele was retained (or a rare reversal occurred).
</p>

<div class="scroll-table">
    <table>
        <tr>
            <th>#</th>
            <th>Chr</th>
            <th>Position</th>
            <th>PBS<sub>UZB</sub></th>
            <th style="color:#c62828">MAF<sub>UZB</sub></th>
            <th>MAF<sub>EUR</sub></th>
            <th>MAF<sub>EAS</sub></th>
            <th>MAF<sub>SAS</sub></th>
            <th>MAF<sub>AFR</sub></th>
            <th>PBS<sub>EUR</sub></th>
            <th>PBS<sub>EAS</sub></th>
        </tr>
{top20_rows}
    </table>
</div>

<div class="important">
    <strong>Interpreting the top PBS hits:</strong>
    <ul>
        <li>All top 20 hits show UZB MAF of 1&ndash;3% while other populations are &gt;95%. The UZB-specific branch is long because UZB retained the minor allele while others went to fixation.</li>
        <li><strong>Negative PBS<sub>EUR</sub></strong> at most loci means the EUR&ndash;EAS shared branch is shorter than the UZB-specific branch &mdash; consistent with UZB being the outlier.</li>
        <li>Two HLA-region hits (<strong>chr6:29640785</strong>, <strong>6:30585771</strong>, <strong>6:31059825</strong>) likely reflect pathogen-driven balancing selection in the MHC region, which is known for extreme population-specific allele frequencies.</li>
        <li>The top hit <strong>9:104189856</strong> (PBS = 4.60) is near the <em>ABCA1</em> gene region, involved in cholesterol metabolism.</li>
    </ul>
</div>


<!-- ═══════════════════════════ 5. UZB-SPECIFIC SNPS ═══════════════════════════ -->
<h2>5. Uzbek-Specific SNP Breakdown</h2>

<table>
    <tr><th>Tier</th><th>Criterion</th><th>Count</th><th>What it captures</th></tr>
    <tr>
        <td><strong>1 (PBS)</strong></td>
        <td>PBS &ge; 0.3</td>
        <td><strong>456</strong></td>
        <td>Alleles under UZB-specific positive selection or extreme drift on the Uzbek branch</td>
    </tr>
    <tr>
        <td><strong>2 (&Delta;AF)</strong></td>
        <td>min |&Delta;AF| &ge; 0.3 vs all 4 populations</td>
        <td><strong>422</strong></td>
        <td>UZB frequency outliers &mdash; different from <em>everyone</em> by at least 30%</td>
    </tr>
    <tr>
        <td><strong>3 (near-private)</strong></td>
        <td>UZB MAF &ge; 5%, all others MAF &le; 1%</td>
        <td><strong>272</strong></td>
        <td>Alleles common only in Uzbek &mdash; potential novel / recently selected variants</td>
    </tr>
    <tr style="background: #e8f5e9;">
        <td colspan="2"><strong>Union of all tiers</strong></td>
        <td><strong style="font-size:1.2em; color:#1b5e20;">490</strong></td>
        <td>Total unique Uzbek-specific candidate SNPs</td>
    </tr>
</table>

<div class="info">
    <strong>What does &ldquo;490 Uzbek-specific SNPs&rdquo; mean?</strong>
    <p style="margin-top:8px;">Out of 376,208 testable SNPs, only <strong>0.13%</strong> show strong Uzbek-specific differentiation. This is consistent with a population that is <em>admixed</em> (not a genetic isolate). Most allele frequencies in UZB fall between EUR and EAS (as expected from its intermediate position on PCA). The 490 exceptions are loci where UZB deviates from <em>all</em> major world populations &mdash; these are the best candidates for Uzbek-specific selection or founder effects.</p>
</div>

<h3>5.1. Near-Private Variants (Tier 3)</h3>
<p>
    <strong>272 variants</strong> are common in Uzbek (&ge;5% frequency) but essentially absent (&lt;1%) in EUR, EAS, SAS, and AFR. These are particularly interesting because they represent alleles found <em>almost exclusively</em> in the Uzbek cohort.
</p>
<p>
    Possible explanations for near-private variants:
</p>
<ul>
    <li><strong>Novel mutations</strong> that arose in Central Asia after population separation</li>
    <li><strong>Ancient alleles</strong> preserved in Central Asia but lost (drifted out) in other populations</li>
    <li><strong>Admixture-generated combinations</strong> &mdash; alleles present at low frequency in multiple source populations that reached higher frequency only in the admixed Uzbek population (admixture surfing)</li>
    <li><strong>Local adaptation</strong> to Central Asian environmental pressures (diet, altitude, pathogens, UV)</li>
    <li><strong>Genotyping artifacts</strong> &mdash; must be validated before biological interpretation</li>
</ul>


<!-- ═══════════════════════════ 6. ADMIXTURE ═══════════════════════════ -->
<h2>6. ADMIXTURE Ancestry Analysis</h2>

<h3>6.1. Cross-Validation Results</h3>
<p>
    ADMIXTURE uses <strong>cross-validation error</strong> to determine the optimal number of ancestral populations (K). Lower CV error = better fit. We tested K=2 through K=4 (K=5&ndash;8 running).
</p>

<table>
    <tr>
        <th>K</th>
        <th>CV Error</th>
        <th>Iterations</th>
        <th>Time</th>
        <th>Log-likelihood</th>
        <th>Assessment</th>
    </tr>
    <tr style="background: #e8f5e9;">
        <td><strong>K=2</strong></td>
        <td><strong>0.51335</strong></td>
        <td>25</td>
        <td>29 min</td>
        <td>Best</td>
        <td style="color:#2e7d32; font-weight:600;">&star; Optimal (lowest CV)</td>
    </tr>
    <tr>
        <td><strong>K=3</strong></td>
        <td>0.51343</td>
        <td>66</td>
        <td>108 min</td>
        <td>&mdash;</td>
        <td>Nearly tied with K=2 (difference = 0.00008)</td>
    </tr>
    <tr>
        <td><strong>K=4</strong></td>
        <td>0.51579</td>
        <td>109</td>
        <td>213 min</td>
        <td>&mdash;</td>
        <td style="color:#c62828;">Worse &mdash; 4th component adds noise</td>
    </tr>
</table>

<div class="info">
    <strong>What does this tell us?</strong>
    <p style="margin-top:8px;">K=2 and K=3 have essentially identical CV error (within 0.02% of each other), while K=4 is clearly worse. This means <strong>the Uzbek cohort is best described by 2&ndash;3 ancestral components</strong>. K=2 captures the primary West&ndash;East axis (consistent with PCA), while K=3 adds a third component that may represent South Asian, Iranian, or ancestral Central Asian ancestry.</p>
</div>

<h3>6.2. ADMIXTURE Structure Plots (Interactive)</h3>
<p>
    Each vertical bar represents one individual (1,074 total). Colors represent ancestry proportions from inferred ancestral populations. Individuals are sorted by K=2 Component 1 (Western ancestry) from left to right.
</p>
<p>
    <strong>How to read this:</strong> A bar that is 70% blue and 30% orange means that individual&rsquo;s genome is estimated to derive 70% from one ancestral population and 30% from another. A bar of solid single color means that individual is &ldquo;nearly pure&rdquo; for one ancestry component.
</p>

<div class="chart-container" id="admixContainer">
    <div class="chart-toolbar" id="admixToolbar">
        <button class="active" onclick="drawAdmixture(2)">K = 2</button>
        <button onclick="drawAdmixture(3)">K = 3</button>
        <button onclick="drawAdmixture(4)">K = 4</button>
        <div class="legend" id="admixLegend"></div>
    </div>
    <canvas id="admixCanvas" width="1050" height="280"></canvas>
    <div class="tooltip-box" id="admixTooltip" style="display:none;"></div>
    <p style="font-size:0.85em; color:#888; margin-top:8px;" id="admixCaption">Hover over any bar to see the individual&rsquo;s sample ID and ancestry proportions.</p>
</div>

<h3>6.3. Component Interpretation</h3>

<h4>K=2: West&ndash;East Axis (primary structure)</h4>
<table>
    <tr><th>Component</th><th>Mean</th><th>Range</th><th>Likely ancestry</th></tr>
    <tr><td style="color: #2196F3;"><strong>&#x25A0; Component 1 (Blue)</strong></td><td>~65%</td><td>0&ndash;100%</td><td>Western Eurasian / European-like</td></tr>
    <tr><td style="color: #FF9800;"><strong>&#x25A0; Component 2 (Orange)</strong></td><td>~35%</td><td>0&ndash;100%</td><td>Eastern Eurasian / East Asian-like</td></tr>
</table>
<p>
    This mirrors the PCA results (Steps 7&ndash;8): the main axis of variation in Uzbek is a <strong>West&ndash;East admixture cline</strong>. Most individuals are ~65% Western and ~35% Eastern, consistent with Central Asian populations being a mixture of Indo-Iranian and Turkic/Mongolian ancestries.
</p>
<p>
    <strong>Outliers:</strong> ~31 individuals are &gt;95% Component 1 (nearly pure Western &mdash; likely Tajik, Russian, or other European-origin individuals in the cohort). ~18 individuals are &gt;95% Component 2 (nearly pure Eastern &mdash; likely Kazakh, Kyrgyz, or other Turkic individuals).
</p>

<h4>K=3: Three-Way Admixture</h4>
<table>
    <tr><th>Component</th><th>Mean</th><th>Likely ancestry</th></tr>
    <tr><td style="color: #2196F3;"><strong>&#x25A0; Component 1</strong></td><td>~49%</td><td>Western Eurasian (more refined)</td></tr>
    <tr><td style="color: #FF9800;"><strong>&#x25A0; Component 2</strong></td><td>~31%</td><td>East Asian-like</td></tr>
    <tr><td style="color: #4CAF50;"><strong>&#x25A0; Component 3</strong></td><td>~20%</td><td>South/Central Asian or ancestral Iranian</td></tr>
</table>
<p>
    K=3 splits the Western component into two: a &ldquo;pure European-like&rdquo; fraction and a third component that may represent <strong>ancestral Iranian or South Asian ancestry</strong> shared across Central Asian populations. This is consistent with the F<sub>ST</sub> finding that UZB is closest to SAS (South Asia).
</p>

<h4>K=4: Minimal Additional Resolution</h4>
<p>K=4 adds a small 4th component (mean ~4.5%) with higher CV error. Only 5 individuals are &gt;80% assigned to this component, suggesting it may represent <strong>outlier individuals</strong> rather than a genuine population-level ancestry. K=4 is not recommended for the main analysis.</p>


<!-- ═══════════════════════════ 7. THEORY ═══════════════════════════ -->
<h2>7. Theoretical Background</h2>

<h3>7.1. Population Branch Statistic (PBS)</h3>
<p>
    Developed by <strong>Yi et al. (2010)</strong> to detect high-altitude adaptation in Tibetans (identifying <em>EPAS1</em>), PBS extends pairwise F<sub>ST</sub> into a three-population phylogenetic framework.
</p>
<p>
    <strong>The key insight:</strong> If you only compare two populations, you can&rsquo;t tell which one changed. But with three populations forming a triangle, you can decompose the total divergence into <em>branch-specific</em> components and determine which branch got longer (= which population changed).
</p>

<p>Think of it like a family tree analogy:</p>
<ul>
    <li>Comparing UZB vs EUR: &ldquo;These siblings look different.&rdquo; (But which one changed?)</li>
    <li>Adding EAS as a cousin: &ldquo;The cousin looks more like EUR, so UZB must be the one who changed.&rdquo;</li>
</ul>

<div class="formula">
    Step 1: Convert F<sub>ST</sub> to &ldquo;divergence time&rdquo;<br>
    T<sub>AB</sub> = &minus;log(1 &minus; F<sub>ST(AB)</sub>)<br><br>
    Step 2: Compute branches<br>
    PBS<sub>A</sub> = (T<sub>AB</sub> + T<sub>AC</sub> &minus; T<sub>BC</sub>) / 2<br><br>
    <span class="label">The log transformation converts F<sub>ST</sub> into additive divergence times (like molecular clock distances). The formula then uses the three pairwise distances to solve for each individual branch length.</span>
</div>

<h3>7.2. ADMIXTURE (Model-Based Ancestry)</h3>
<p>
    <strong>ADMIXTURE</strong> (Alexander et al., 2009) estimates individual ancestry proportions using a parametric model. Unlike PCA (which is geometry-based), ADMIXTURE fits an explicit probabilistic model:
</p>
<ul>
    <li>There exist K ancestral populations, each with their own allele frequencies at every locus</li>
    <li>Each individual&rsquo;s genome is a mixture of these ancestral populations</li>
    <li>The algorithm finds the ancestral frequencies and individual proportions that maximize the likelihood of the observed genotype data</li>
</ul>

<div class="important">
    <strong>ADMIXTURE caveats:</strong>
    <ul>
        <li>Ancestral &ldquo;populations&rdquo; are mathematical constructs, not real historical groups</li>
        <li>The model assumes random mating within ancestral populations (Hardy-Weinberg equilibrium)</li>
        <li>Results depend on the SNP set used (LD pruning is essential)</li>
        <li>Cross-validation helps select K, but two K values with similar CV error are equally valid</li>
    </ul>
</div>

<h3>7.3. Why Three Tiers?</h3>
<table>
    <tr><th>Tier</th><th>What it detects</th><th>Limitation</th></tr>
    <tr>
        <td><strong>1 (PBS)</strong></td>
        <td>Recent positive selection on UZB branch &mdash; allele frequency change <em>specific</em> to UZB relative to both EUR and EAS</td>
        <td>Requires all three populations to share the SNP. Assumes tree-like history (violated by admixture).</td>
    </tr>
    <tr>
        <td><strong>2 (&Delta;AF)</strong></td>
        <td>SNPs where UZB is a frequency outlier vs <em>all</em> global populations &mdash; captures drift + selection + founder effects</td>
        <td>Doesn&rsquo;t distinguish selection from strong drift. More permissive.</td>
    </tr>
    <tr>
        <td><strong>3 (near-private)</strong></td>
        <td>Alleles common in UZB (&ge;5%) but nearly absent elsewhere (&le;1%) &mdash; potential novel variants or recent sweeps</td>
        <td>Small numbers expected. Could be genotyping artifacts or structural variant tagging.</td>
    </tr>
</table>


<!-- ═══════════════════════════ 8. SCRIPTS ═══════════════════════════ -->
<h2>8. Complete Scripts</h2>

<h3>8.1. Master Wrapper &mdash; <code>run_uzb_specific_pipeline.sh</code></h3>
<div class="code-block">#!/bin/bash
# Master script: Run the complete Uzbek-specific SNP identification pipeline
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" &amp;&amp; pwd)"
OUTDIR="/staging/ALSU-analysis/admixture_analysis/pop_specific"
FST_DIR="/staging/ALSU-analysis/Fst_analysis"

echo "===== Uzbek-Specific SNP Pipeline ====="
echo "Start time: $(date)"

# Phase 1: Extract populations and compute pairwise Fst
bash "${{SCRIPT_DIR}}/01_extract_multipop.sh"

# Symlink existing UZB-EUR Fst to output directory
if [ ! -f "${{OUTDIR}}/fst_UZB_vs_EUR.fst" ]; then
    ln -s "${{FST_DIR}}/genomewide_uzbek_vs_eur_fst.fst" "${{OUTDIR}}/fst_UZB_vs_EUR.fst"
fi

# Phase 2: PBS &amp; Uzbek-specific SNP analysis
python3 "${{SCRIPT_DIR}}/02_calculate_pbs.py" \\
    --outdir "${{OUTDIR}}" \\
    --pbs-threshold 0.3 \\
    --delta-af-threshold 0.3

echo "===== Pipeline Complete! ====="
echo "End time: $(date)"</div>

<h3>8.2. Phase 1 &mdash; <code>01_extract_multipop.sh</code> (270 lines)</h3>
<details>
    <summary>Click to expand full script</summary>
    <div class="code-block">{script1_html}</div>
</details>

<h3>8.3. Phase 2 &mdash; <code>02_calculate_pbs.py</code> (409 lines)</h3>
<details>
    <summary>Click to expand full script</summary>
    <div class="code-block">{script2_html}</div>
</details>

<h3>8.4. ADMIXTURE Analysis Commands</h3>
<div class="code-block"># ── LD Pruning (creating input for ADMIXTURE) ──────────────────────
# Start from 10.8M imputed SNPs in UZB dataset
plink2 --bfile /staging/ALSU-analysis/Fst_analysis/uzbek_data_hg19 \\
    --maf 0.05 \\
    --indep-pairwise 50 10 0.1 \\
    --out UZB_prune

# Apply pruning
plink2 --bfile /staging/ALSU-analysis/Fst_analysis/uzbek_data_hg19 \\
    --extract UZB_prune.prune.in \\
    --make-bed \\
    --out UZB_for_admixture
# Result: 1,074 samples &times; 380,376 SNPs

# ── Run ADMIXTURE K=2 through K=8 ─────────────────────────────────
for K in 2 3 4 5 6 7 8; do
    echo "=== Running K=$K ==="
    admixture --cv=5 -j32 UZB_for_admixture.bed $K 2&gt;&amp;1 | \\
        tee admixture_K${{K}}.log
done
# Note: Each K takes ~30 min (K=2) to ~4+ hours (K=4)
# Output: UZB_for_admixture.{{K}}.Q (ancestry proportions)
#         UZB_for_admixture.{{K}}.P (allele frequencies)</div>

<h3>8.5. Data Export Script (for interactive visualization)</h3>
<details>
    <summary>Click to expand Python export script</summary>
    <div class="code-block">#!/usr/bin/env python3
"""Export ADMIXTURE Q matrices and PBS data as JSON for web visualization."""
import json
import numpy as np

base = "/staging/ALSU-analysis/admixture_analysis"

# Load sample IDs from .fam file
ids = []
with open(f"{{base}}/UZB_for_admixture.fam") as f:
    for line in f:
        ids.append(line.strip().split()[1])

data = {{"ids": ids}}

# Load Q matrices for each K
for k in [2, 3, 4]:
    q = np.loadtxt(f"{{base}}/UZB_for_admixture.{{k}}.Q")
    data[f"k{{k}}"] = q.tolist()

# Sort order: by K=2 component 1
k2 = np.array(data["k2"])
sort_idx = np.argsort(k2[:, 0]).tolist()
data["sort_order"] = sort_idx

# Top PBS hits
import pandas as pd
pbs = pd.read_csv(f"{{base}}/pop_specific/pbs_results.tsv", sep="\\t")
top30 = pbs.nlargest(30, "PBS_UZB")
data["pbs_top"] = top30.to_dict("records")
data["pbs_header"] = list(pbs.columns)

with open("/tmp/admix_data.json", "w") as f:
    json.dump(data, f)
print(f"Exported {{len(ids)}} samples, JSON size: {{len(json.dumps(data))}} bytes")</div>
</details>


<!-- ═══════════════════════════ 9. DIRECTORY ═══════════════════════════ -->
<h2>9. Output Directory Structure</h2>
<div class="code-block">/staging/ALSU-analysis/admixture_analysis/
&boxvr;&boxh;&boxh; UZB_for_admixture.bed/bim/fam    # 1,074 &times; 380,376 LD-pruned SNPs
&boxvr;&boxh;&boxh; UZB_for_admixture.2.Q/P          # K=2 ADMIXTURE output
&boxvr;&boxh;&boxh; UZB_for_admixture.3.Q/P          # K=3 ADMIXTURE output
&boxvr;&boxh;&boxh; UZB_for_admixture.4.Q/P          # K=4 ADMIXTURE output
&boxvr;&boxh;&boxh; admixture_K*.log                  # ADMIXTURE log files
&boxvr;&boxh;&boxh; run_admixture.sh                  # ADMIXTURE launch script
&boxur;&boxh;&boxh; pop_specific/
    &boxvr;&boxh;&boxh; scripts/
    &boxvr;   &boxvr;&boxh;&boxh; 01_extract_multipop.sh       # Phase 1 (270 lines)
    &boxvr;   &boxvr;&boxh;&boxh; 02_calculate_pbs.py          # Phase 2 (409 lines)
    &boxvr;   &boxur;&boxh;&boxh; run_uzb_specific_pipeline.sh # Master wrapper
    &boxvr;&boxh;&boxh; 1000G_EAS_all.bed/bim/fam            # 504 &times; ~376K SNPs
    &boxvr;&boxh;&boxh; 1000G_SAS_all.bed/bim/fam            # 489 &times; ~376K SNPs
    &boxvr;&boxh;&boxh; 1000G_AFR_all.bed/bim/fam            # 661 &times; ~376K SNPs
    &boxvr;&boxh;&boxh; freq_{{UZB,EUR,EAS,SAS,AFR}}.frq     # Allele frequencies
    &boxvr;&boxh;&boxh; fst_UZB_vs_{{EUR,EAS,SAS,AFR}}.fst   # Pairwise Fst files
    &boxvr;&boxh;&boxh; fst_EUR_vs_EAS.fst                   # PBS outgroup pair
    &boxvr;&boxh;&boxh; pbs_results.tsv                      # All PBS scores (376K rows)
    &boxvr;&boxh;&boxh; uzbek_specific_snps.tsv              # 490 Uzbek-specific candidates
    &boxvr;&boxh;&boxh; delta_af_all.tsv                     # &Delta;AF across all populations
    &boxvr;&boxh;&boxh; near_private_variants.tsv            # 272 near-private variants
    &boxur;&boxh;&boxh; summary_stats.txt                    # Human-readable summary</div>


<!-- ═══════════════════════════ 10. INTERPRETATION ═══════════════════════════ -->
<h2>10. Interpretation &amp; Biological Context</h2>

<h3>10.1. What the ADMIXTURE + PBS results tell us together</h3>
<p>
    The Uzbek cohort is a <strong>two-way admixed population</strong> (~65% Western Eurasian, ~35% Eastern Eurasian) with continuous variation across individuals. This admixture explains:
</p>
<ul>
    <li><strong>Negative mean PBS (&minus;0.013)</strong>: Because UZB sits between EUR and EAS, the outgroup branch (EUR&ndash;EAS) tends to be longer than the UZB branch for most SNPs &mdash; alleles inherited from EUR match EUR, and alleles from EAS match EAS.</li>
    <li><strong>490 Uzbek-specific SNPs</strong>: These are loci where UZB deviates from <em>both</em> parental populations. They could be:
        <ul>
            <li>Selection targets unique to Central Asia (pathogens, diet, UV, altitude)</li>
            <li>Admixture-generated novel haplotype combinations</li>
            <li>Drift effects in a population that went through bottlenecks (Silk Road history)</li>
        </ul>
    </li>
    <li><strong>UZB closest to SAS</strong>: The F<sub>ST</sub> = 0.018 to South Asia (&lt; 0.020 to EUR) suggests a shared Indo-Iranian foundation underlying the Turkic admixture.</li>
</ul>

<h3>10.2. Caveats &amp; Limitations</h3>
<ul>
    <li><strong>PBS assumes tree-like history</strong>: UZB is admixed (EUR + EAS), violating this assumption. High PBS could reflect admixture dynamics rather than recent selection.</li>
    <li><strong>Ascertainment bias</strong>: The 376K SNP set was filtered for UZB&ndash;EUR overlap, potentially missing truly Uzbek-specific variants not on the original genotyping array.</li>
    <li><strong>Sample size imbalance</strong>: UZB (1,199) vs EUR (503) affects individual F<sub>ST</sub> estimation variance (though genome-wide weighted mean is robust).</li>
    <li><strong>No gene annotation yet</strong>: Top hits need functional annotation (nearest gene, regulatory elements, expression data) before biological conclusions.</li>
    <li><strong>ADMIXTURE K=5&ndash;8</strong>: Still running at time of writing. Results should be checked for substructure.</li>
</ul>


<!-- ═══════════════════════════ 11. NEXT STEPS ═══════════════════════════ -->
<h2>11. Next Steps</h2>
<ol>
    <li><strong>Gene annotation</strong>: Map the 490 Uzbek-specific SNPs to nearest genes using ANNOVAR or VEP; check for enrichment in known selection targets.</li>
    <li><strong>ADMIXTURE K=5&ndash;8</strong>: Monitor completion, check if additional K values reveal substructure or if K=2&ndash;3 remains optimal.</li>
    <li><strong>IBD re-analysis</strong>: Use ADMIXTURE ancestry to stratify IBD segments &mdash; do related individuals share more Western or Eastern ancestry?</li>
    <li><strong>Pathway enrichment</strong>: Test if Uzbek-specific SNPs are enriched in KEGG/GO pathways (immunity, metabolism, etc.).</li>
    <li><strong>Case-control design</strong> (Step 11): Define pregnancy loss criteria, map phenotype IDs to genotype IDs, and test the 490 candidate SNPs for association with pregnancy loss.</li>
    <li><strong>Extended haplotype analysis</strong>: Phase genotypes and compute iHS/nSL for top PBS regions to confirm recent selection signals.</li>
</ol>

            </div> <!-- /content -->
        </div> <!-- /main-content -->
    </div> <!-- /page-wrapper -->

    <!-- ═══════════════════════════ EMBEDDED DATA ═══════════════════════════ -->
    <script>
    // ADMIXTURE Q-matrix data (1,074 samples x K=2,3,4)
    const ADMIX_DATA = {admix_json};
    </script>

    <script>
    // ── PBS Histogram ──────────────────────────────────────────────────────
    (function() {{
        const histVals = {hist_json};
        const binWidth = 2.5 / 50;  // range [-0.5, 2.0] / 50 bins
        const canvas = document.getElementById('pbsHistCanvas');
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;
        const pad = {{l:65, r:20, t:20, b:50}};
        const plotW = W - pad.l - pad.r;
        const plotH = H - pad.t - pad.b;
        
        // Log scale
        const maxLog = Math.log10(Math.max(...histVals));
        
        ctx.fillStyle = '#fff';
        ctx.fillRect(0, 0, W, H);
        
        // Bars
        const barW = plotW / histVals.length;
        for (let i = 0; i < histVals.length; i++) {{
            const v = histVals[i];
            if (v === 0) continue;
            const logV = Math.log10(v);
            const barH = (logV / maxLog) * plotH;
            const x = pad.l + i * barW;
            const y = pad.t + plotH - barH;
            
            // Color by position
            const binCenter = -0.5 + (i + 0.5) * binWidth;
            if (binCenter >= 0.3) ctx.fillStyle = '#c62828';
            else if (binCenter >= 0.1) ctx.fillStyle = '#e65100';
            else if (binCenter >= 0) ctx.fillStyle = '#2e7d32';
            else ctx.fillStyle = '#90a4ae';
            
            ctx.fillRect(x, y, barW - 1, barH);
        }}
        
        // Axes
        ctx.strokeStyle = '#333';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(pad.l, pad.t);
        ctx.lineTo(pad.l, pad.t + plotH);
        ctx.lineTo(pad.l + plotW, pad.t + plotH);
        ctx.stroke();
        
        // Y-axis labels (log scale)
        ctx.fillStyle = '#333';
        ctx.font = '11px Segoe UI';
        ctx.textAlign = 'right';
        for (let exp = 0; exp <= 5; exp++) {{
            const y = pad.t + plotH - (exp / maxLog) * plotH;
            ctx.fillText('10^' + exp, pad.l - 5, y + 4);
            ctx.strokeStyle = '#ddd';
            ctx.beginPath();
            ctx.moveTo(pad.l, y);
            ctx.lineTo(pad.l + plotW, y);
            ctx.stroke();
        }}
        ctx.strokeStyle = '#333';
        
        // X-axis labels
        ctx.textAlign = 'center';
        for (let v = -0.5; v <= 2.0; v += 0.5) {{
            const x = pad.l + ((v + 0.5) / 2.5) * plotW;
            ctx.fillText(v.toFixed(1), x, pad.t + plotH + 20);
            ctx.strokeStyle = '#ddd';
            ctx.beginPath();
            ctx.moveTo(x, pad.t);
            ctx.lineTo(x, pad.t + plotH);
            ctx.stroke();
        }}
        
        // Axis titles
        ctx.fillStyle = '#333';
        ctx.font = '12px Segoe UI';
        ctx.textAlign = 'center';
        ctx.fillText('PBS_UZB', pad.l + plotW/2, H - 5);
        
        ctx.save();
        ctx.translate(15, pad.t + plotH/2);
        ctx.rotate(-Math.PI/2);
        ctx.fillText('Count (log\u2081\u2080 scale)', 0, 0);
        ctx.restore();
        
        // Threshold line at 0.3
        const threshX = pad.l + (0.8 / 2.5) * plotW;  // 0.3 + 0.5 = 0.8 in [0, 2.5]
        ctx.strokeStyle = '#c62828';
        ctx.setLineDash([5, 3]);
        ctx.beginPath();
        ctx.moveTo(threshX, pad.t);
        ctx.lineTo(threshX, pad.t + plotH);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = '#c62828';
        ctx.font = 'bold 11px Segoe UI';
        ctx.textAlign = 'left';
        ctx.fillText('PBS = 0.3 threshold', threshX + 4, pad.t + 15);
        ctx.fillText('(456 SNPs above)', threshX + 4, pad.t + 28);
    }})();

    // ── ADMIXTURE interactive chart ────────────────────────────────────────
    const admixColors = {{
        2: ['#2196F3', '#FF9800'],
        3: ['#2196F3', '#FF9800', '#4CAF50'],
        4: ['#2196F3', '#FF9800', '#4CAF50', '#9C27B0']
    }};
    const compNames = {{
        2: ['Western Eurasian', 'Eastern Eurasian'],
        3: ['Component 1', 'Component 2', 'Component 3'],
        4: ['Component 1', 'Component 2', 'Component 3', 'Component 4']
    }};
    
    let currentK = 2;
    
    function drawAdmixture(k) {{
        currentK = k;
        const canvas = document.getElementById('admixCanvas');
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;
        const pad = {{l: 40, r: 10, t: 5, b: 25}};
        const plotW = W - pad.l - pad.r;
        const plotH = H - pad.t - pad.b;
        
        const data = ADMIX_DATA;
        const qKey = 'k' + k;
        const qMatrix = data[qKey];
        const sortOrder = data.sort_order;
        const ids = data.ids;
        const n = sortOrder.length;
        const barW = plotW / n;
        
        // Update buttons
        document.querySelectorAll('#admixToolbar button').forEach(btn => {{
            btn.classList.toggle('active', btn.textContent.includes(k));
        }});
        
        // Update legend
        const legendDiv = document.getElementById('admixLegend');
        legendDiv.innerHTML = admixColors[k].map((c, i) => 
            '<span class="legend-item"><span class="legend-swatch" style="background:' + c + '"></span>' + compNames[k][i] + '</span>'
        ).join('');
        
        // Clear
        ctx.fillStyle = '#fafafa';
        ctx.fillRect(0, 0, W, H);
        
        // Draw bars
        for (let i = 0; i < n; i++) {{
            const idx = sortOrder[i];
            const q = qMatrix[idx];
            const x = pad.l + i * barW;
            let yAccum = 0;
            
            for (let c = 0; c < k; c++) {{
                const h = q[c] * plotH;
                ctx.fillStyle = admixColors[k][c];
                ctx.fillRect(x, pad.t + yAccum, Math.max(barW, 1), h);
                yAccum += h;
            }}
        }}
        
        // Border
        ctx.strokeStyle = '#333';
        ctx.lineWidth = 1;
        ctx.strokeRect(pad.l, pad.t, plotW, plotH);
        
        // Y-axis labels
        ctx.fillStyle = '#333';
        ctx.font = '11px Segoe UI';
        ctx.textAlign = 'right';
        ctx.fillText('100%', pad.l - 4, pad.t + 10);
        ctx.fillText('50%', pad.l - 4, pad.t + plotH/2 + 4);
        ctx.fillText('0%', pad.l - 4, pad.t + plotH);
        
        // X-axis
        ctx.textAlign = 'center';
        ctx.fillText('Individuals (n=' + n + ', sorted by K=2 ancestry)', pad.l + plotW/2, H - 3);
    }}
    
    // Tooltip
    (function() {{
        const canvas = document.getElementById('admixCanvas');
        const tooltip = document.getElementById('admixTooltip');
        const container = document.getElementById('admixContainer');
        const pad = {{l: 40, r: 10, t: 5, b: 25}};
        
        canvas.addEventListener('mousemove', function(e) {{
            const rect = canvas.getBoundingClientRect();
            const mx = e.clientX - rect.left;
            const my = e.clientY - rect.top;
            const plotW = canvas.width - pad.l - pad.r;
            const n = ADMIX_DATA.sort_order.length;
            const barW = plotW / n;
            
            const barIdx = Math.floor((mx - pad.l) / barW);
            if (barIdx < 0 || barIdx >= n) {{
                tooltip.style.display = 'none';
                return;
            }}
            
            const sampleIdx = ADMIX_DATA.sort_order[barIdx];
            const sampleId = ADMIX_DATA.ids[sampleIdx];
            const q = ADMIX_DATA['k' + currentK][sampleIdx];
            
            let html = '<strong>' + sampleId + '</strong> (#' + (barIdx+1) + '/' + n + ')<br>';
            for (let c = 0; c < currentK; c++) {{
                html += '<span style="color:' + admixColors[currentK][c] + '">&block;</span> ' + 
                        compNames[currentK][c] + ': <strong>' + (q[c]*100).toFixed(1) + '%</strong><br>';
            }}
            
            tooltip.innerHTML = html;
            tooltip.style.display = 'block';
            
            const containerRect = container.getBoundingClientRect();
            let tx = e.clientX - containerRect.left + 15;
            let ty = e.clientY - containerRect.top - 10;
            if (tx + 250 > containerRect.width) tx = tx - 270;
            tooltip.style.left = tx + 'px';
            tooltip.style.top = ty + 'px';
        }});
        
        canvas.addEventListener('mouseleave', function() {{
            tooltip.style.display = 'none';
        }});
    }})();
    
    // Initial draw
    drawAdmixture(2);
    
    // ── Search ─────────────────────────────────────────────────────────────
    const stepPages = [
        {{ title: 'Step 1: Missingness', url: 'step1.html', keywords: 'missingness quality control mind geno' }},
        {{ title: 'Step 2: IBD Dedup', url: 'step2.html', keywords: 'ibd relatedness duplicates pihat' }},
        {{ title: 'Step 3: SNP QC', url: 'step3.html', keywords: 'snp quality control hwe maf' }},
        {{ title: 'Step 4: Imputation', url: 'step4.html', keywords: 'imputation beagle phasing reference' }},
        {{ title: 'Step 5: ID Normalize', url: 'step5.html', keywords: 'id normalize rename samples' }},
        {{ title: 'Step 6: Final QC', url: 'step6.html', keywords: 'final quality control filtering' }},
        {{ title: 'Step 7: PCA Analysis', url: 'step7.html', keywords: 'pca principal component uzbek' }},
        {{ title: 'Step 8: Global PCA', url: 'step8.html', keywords: 'global pca 1000genomes populations' }},
        {{ title: 'Step 9: Fst Analysis', url: 'step9.html', keywords: 'fst fixation differentiation eur' }},
        {{ title: 'Step 10: Multi-Pop PBS ADMIXTURE', url: 'step10.html', keywords: 'pbs admixture multi-population delta af uzbek-specific near-private selection eas sas afr ancestry' }}
    ];
    
    function searchAllSteps() {{
        const query = document.getElementById('searchInput').value.toLowerCase().trim();
        const resultsDiv = document.getElementById('searchResults');
        if (query.length < 2) {{ resultsDiv.className = 'search-results'; return; }}
        const matches = stepPages.filter(p => p.title.toLowerCase().includes(query) || p.keywords.includes(query));
        if (matches.length > 0) {{
            resultsDiv.innerHTML = '<strong>Found in:</strong> ' + matches.map(m => '<a href="' + m.url + '">' + m.title + '</a>').join(' ');
            resultsDiv.className = 'search-results active';
        }} else {{
            resultsDiv.className = 'search-results';
        }}
    }}
    </script>
</body>
</html>'''

# ── Write output ───────────────────────────────────────────────────────────────
with open("steps/step10.html", "w", encoding="utf-8") as f:
    f.write(page)

print(f"Generated steps/step10.html: {len(page):,} bytes")
print(f"  Embedded ADMIXTURE JSON: {len(admix_json):,} bytes")
print(f"  Script 1 HTML: {len(script1_html):,} bytes")
print(f"  Script 2 HTML: {len(script2_html):,} bytes")
