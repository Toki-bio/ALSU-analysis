import os, re, subprocess, html, json
S = "C:/work/ALSU-analysis/"
STYLE = open(S + "steps/deprecated_old_runs.html", encoding="utf-8").read()
STYLE = STYLE[STYLE.index("<style>"):STYLE.index("</style>") + 8]
orph = ["START_HERE.html", "old_logs/alsu_pipeline_report.html", "old_logs/alsu_pipeline_report_graph.html", "steps/step1_old.html", "steps/step2_old.html"]
arch = sorted("steps/archive/" + f for f in os.listdir(S + "steps/archive") if f.endswith(".html") and not f.startswith("index") and "original_cohort_1047" not in f)
orig = sorted("steps/archive/" + f for f in os.listdir(S + "steps/archive") if f.endswith(".html") and "original_cohort_1047" in f)


def git_date(p):
    r = subprocess.run(["git", "-C", S, "log", "-1", "--format=%ad", "--date=short", "--", p], capture_output=True, text=True)
    return r.stdout.strip() or "?"


def title(p):
    t = open(S + p, encoding="utf-8", errors="ignore").read()
    m = re.search(r"<title>(.*?)</title>", t, flags=re.S)
    return html.unescape(re.sub(r"\s+", " ", m.group(1)).strip()) if m else p


def kind(p):
    n = os.path.basename(p).lower()
    if "original_cohort" in n:
        return "original-cohort (1,047) version of a step page; linked from the live page"
    if n.startswith("start_here"):
        return "obsolete setup notice (workflow documentation system, 2025)"
    if "alsu_pipeline_report" in n:
        return "historical imputation report, November&ndash;December 2025"
    if n.endswith("_old.html") or "_old" in n:
        return "previous version of the page"
    if "backup" in n:
        return "backup copy"
    if "_new" in n or "restructured" in n or "_clean" in n or "_deploy" in n or "pre_bio" in n:
        return "draft or intermediate version of the page (December 2025 &ndash; January 2026)"
    return "archived page"


def superseded_by(p):
    n = os.path.basename(p)
    m = re.match(r"(step\d+)", n.lower())
    if m:
        return f'<a href="../{m.group(1)}.html">{m.group(1)}</a>'
    if "pipeline_report" in n:
        return '<a href="../../steps/step4.html">step 4</a> (imputation) and following steps'
    return '<a href="../../index.html">pipeline index</a>'


rows = []
for grp, lst in (("Not linked from the pipeline (found by the link audit)", orph + arch),):
    for p in lst:
        rel = os.path.relpath(S + p, S + "steps/archive").replace("\\", "/")
        rows.append((p, f'<tr><td><a href="{rel}">{p}</a></td><td>{html.escape(title(p))[:80]}</td><td>{kind(p)}</td><td>{git_date(p)}</td><td>{superseded_by(p)}</td><td>{os.path.getsize(S + p) // 1024} KB</td></tr>'))
orow = []
for p in orig:
    rel = os.path.basename(p)
    orow.append(f'<tr><td><a href="{rel}">{p}</a></td><td>{html.escape(title(p))[:80]}</td><td>{git_date(p)}</td><td>{superseded_by(p)}</td></tr>')

page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Archive index - ALSU Pipeline</title>
{STYLE}
</head>
<body><div class="wrap">
<p class="back"><a href="../../index.html">&larr; Back to pipeline</a> &middot; <a href="../deprecated_old_runs.html">Deprecated old runs and pending updates</a></p>
<h1>Archive index</h1>
<p>Created 2026-10-04 by the link audit of the site (script <code>linkgraph.py</code>). Before this page, {len(orph) + len(arch)} HTML files were not reachable by any link from the pipeline index. They are drafts, backups and older versions of step pages, an obsolete setup notice and two 2025 imputation reports. They are kept for provenance only, were <strong>not</strong> updated, may contain superseded numbers, and should not be cited. Each carries a banner pointing here. Nothing was deleted.</p>
<h2>1. Files that were not linked from anywhere</h2>
<table><thead><tr><th>File</th><th>Title</th><th>What it is</th><th>Last commit</th><th>Current page</th><th>Size</th></tr></thead><tbody>
{"".join(r for _, r in rows)}
</tbody></table>
<h2>2. Original-cohort (1,047 samples, January 2026) versions of the step pages</h2>
<p>These are linked from the live step pages (&ldquo;view original cohort results&rdquo;). They describe the 1,047-sample cohort on older references. Their internal navigation links pointed to non-existent files inside this folder and were repaired on 2026-10-04 to point to the live pages.</p>
<table><thead><tr><th>File</th><th>Title</th><th>Last commit</th><th>Current page</th></tr></thead><tbody>
{"".join(orow)}
</tbody></table>
<h2>3. What the audit also found</h2>
<ul>
<li>379 broken internal links, all inside the files above (navigation links to <code>stepN.html</code> written relative to the wrong folder). Repaired where the target exists. Live pages had no broken internal links to HTML, Markdown, script, JSON or image files.</li>
<li><code>steps/step5.html</code> had a JavaScript syntax error (extra closing brace) that stopped its page script; <code>steps/step8.html</code> had a truncated search function. Both fixed.</li>
<li>Participant-level identifiers: the &ldquo;extreme individuals&rdquo; table on step 15 showed internal sample codes containing initials; replaced by ranks. One participant name remains in the git history of the repository (earlier commits); it is not on any live page.</li>
</ul>
</div></body></html>
"""
os.makedirs(S + "steps/archive", exist_ok=True)
open(S + "steps/archive/index.html", "w", encoding="utf-8", newline="").write(page)

# banners on the unlinked files
BAN = '<div style="background:#fee2e2;border:2px solid #dc2626;padding:10px 14px;margin:10px 0;font:14px/1.4 Segoe UI,Arial,sans-serif;color:#1f2937"><strong>Archived page (2026-10-04 audit).</strong> This file is a draft, backup or older version that is not part of the current pipeline. Numbers may be outdated or withdrawn. See the <a href="{idx}">archive index</a> and the <a href="{dep}">list of deprecated runs</a>.</div>'
for p in orph + arch:
    t = open(S + p, encoding="utf-8", errors="ignore", newline="").read()
    if "Archived page (2026-10-04 audit)" in t:
        continue
    depth = p.count("/")
    idx = ("../" * depth) + "steps/archive/index.html"
    dep = ("../" * depth) + "steps/deprecated_old_runs.html"
    m = re.search(r"<body[^>]*>", t)
    if not m:
        continue
    t = t[:m.end()] + BAN.format(idx=idx, dep=dep) + t[m.end():]
    open(S + p, "w", encoding="utf-8", newline="").write(t)

# repair nav links inside steps/archive/*.html
fixed = 0
for f in os.listdir(S + "steps/archive"):
    if not f.endswith(".html") or f == "index.html":
        continue
    p = S + "steps/archive/" + f
    t = open(p, encoding="utf-8", errors="ignore", newline="").read()
    def fix(m):
        global fixed
        h = m.group(2)
        if re.match(r"(https?:|mailto:|#|/|javascript:|data:)", h):
            return m.group(0)
        base = h.split("#")[0].split("?")[0]
        if not base:
            return m.group(0)
        d = S + "steps/archive/"
        if os.path.exists(os.path.normpath(d + base)):
            return m.group(0)
        if os.path.exists(os.path.normpath(d + "../" + base)):
            fixed += 1
            return f'{m.group(1)}../{h}'
        if os.path.exists(os.path.normpath(d + "../../" + base)):
            fixed += 1
            return f'{m.group(1)}../../{h}'
        return m.group(0)
    t2 = re.sub(r'((?:href|src)\s*=\s*["\'])([^"\']+)', fix, t)
    if t2 != t:
        open(p, "w", encoding="utf-8", newline="").write(t2)

# loose root files not referenced by any page or script
corpus = ""
for dp, dn, fn in os.walk(S):
    if ".git" in dp.replace("\\", "/").split("/"):
        continue
    for f in fn:
        if f.endswith((".html", ".md", ".py", ".sh", ".json", ".txt")) and os.path.getsize(os.path.join(dp, f)) < 3_000_000:
            try:
                corpus += open(os.path.join(dp, f), encoding="utf-8", errors="ignore").read() + " "
            except Exception:
                pass
loose = sorted(f for f in os.listdir(S) if os.path.isfile(S + f) and not f.endswith(".html") and corpus.count(f) <= 1 and f != ".gitignore")
lrows = "".join(f"<tr><td><code>{f}</code></td><td>{os.path.getsize(S + f) // 1024 or '&lt;1'} KB</td><td>{git_date(f)}</td></tr>" for f in loose)
sec = f"""<h2>4. Loose files in the repository root that no page or script refers to</h2>
<p>{len(loose)} non-HTML files in the repository root (analysis scripts of the superseded V2 runs, one-off check scripts, investigation logs, older notes) are not referenced by any page or other script. They were <strong>not moved or deleted</strong>: the V2 PBS/FST scripts document how the withdrawn results were produced, and moving files is a repository-structure decision for the project owner. Suggested action: move them into a <code>scripts/legacy/</code> folder in one commit.</p>
<table><thead><tr><th>File</th><th>Size</th><th>Last commit</th></tr></thead><tbody>{lrows}</tbody></table>
"""
pg = open(S + "steps/archive/index.html", encoding="utf-8").read()
pg = pg.replace("</div></body></html>", sec + "</div></body></html>")
open(S + "steps/archive/index.html", "w", encoding="utf-8", newline="").write(pg)
print(len(loose), "loose files listed")

print("archive index written;", len(orph) + len(arch), "banners;", fixed, "links repaired")
