#!/usr/bin/env python3
"""Render the admixture-aware relatedness markdown report to a static HTML page."""

from __future__ import annotations

import html
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MD_PATH = ROOT / "ADMIXTURE_AWARE_RELATEDNESS_METHODS_DISCUSSION.md"
OUT_PATH = ROOT / "steps" / "admixture_aware_relatedness.html"


def slugify(text: str) -> str:
    text = re.sub(r"<[^>]+>", "", text)
    text = re.sub(r"[^a-zA-Z0-9]+", "-", text.lower()).strip("-")
    return text or "section"


def inline_markup(text: str) -> str:
    escaped = html.escape(text)
    escaped = re.sub(r"`([^`]+)`", r"<code>\1</code>", escaped)
    escaped = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", escaped)
    escaped = escaped.replace("&gt;=", "&ge;").replace("&lt;=", "&le;")
    return escaped


def close_list(output: list[str], list_type: str | None) -> str | None:
    if list_type:
        output.append(f"</{list_type}>")
    return None


def render_markdown(markdown: str) -> tuple[str, list[tuple[int, str, str]]]:
    output: list[str] = []
    headings: list[tuple[int, str, str]] = []
    list_type: str | None = None
    in_code = False
    code_lang = ""
    code_lines: list[str] = []

    for raw_line in markdown.splitlines():
        line = raw_line.rstrip("\n")
        fence = re.match(r"^```(.*)$", line)
        if fence:
            if in_code:
                output.append(
                    f'<pre><code class="language-{html.escape(code_lang)}">'
                    + html.escape("\n".join(code_lines))
                    + "</code></pre>"
                )
                code_lines = []
                code_lang = ""
                in_code = False
            else:
                list_type = close_list(output, list_type)
                in_code = True
                code_lang = fence.group(1).strip()
            continue

        if in_code:
            code_lines.append(line)
            continue

        if not line.strip():
            list_type = close_list(output, list_type)
            continue

        heading = re.match(r"^(#{1,4})\s+(.+)$", line)
        if heading:
            list_type = close_list(output, list_type)
            level = len(heading.group(1))
            text = heading.group(2).strip()
            section_id = slugify(text)
            headings.append((level, text, section_id))
            output.append(f'<h{level} id="{section_id}">{inline_markup(text)}</h{level}>')
            continue

        unordered = re.match(r"^\s*-\s+(.+)$", line)
        ordered = re.match(r"^\s*\d+\.\s+(.+)$", line)
        if unordered or ordered:
            desired = "ul" if unordered else "ol"
            if list_type != desired:
                list_type = close_list(output, list_type)
                output.append(f"<{desired}>")
                list_type = desired
            output.append(f"<li>{inline_markup((unordered or ordered).group(1))}</li>")
            continue

        list_type = close_list(output, list_type)
        output.append(f"<p>{inline_markup(line)}</p>")

    if in_code:
        output.append("<pre><code>" + html.escape("\n".join(code_lines)) + "</code></pre>")
    close_list(output, list_type)
    return "\n".join(output), headings


def build_sidebar(headings: list[tuple[int, str, str]]) -> str:
    links = ['<a href="step15.html">Step 15: ROH &amp; IBD</a>', '<a href="step11.html">Step 11: ADMIXTURE</a>']
    for level, text, section_id in headings:
        if level == 2:
            links.append(f'<a href="#{section_id}">{html.escape(text)}</a>')
    return "\n".join(links)


def main() -> None:
    markdown = MD_PATH.read_text(encoding="utf-8")
    markdown = re.sub(r"^# .+\n+", "", markdown, count=1)
    body, headings = render_markdown(markdown)
    sidebar = build_sidebar(headings)
    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Admixture-Aware Relatedness — ALSU Analysis</title>
<style>
*{{box-sizing:border-box}}body{{margin:0;font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;background:#f3f5f4;color:#263238}}.page{{display:flex;gap:22px;max-width:1480px;margin:0 auto;padding:22px}}.sidebar{{width:260px;flex:0 0 260px;background:#fff;border:1px solid #dfe6e2;border-radius:8px;padding:18px;height:calc(100vh - 44px);position:sticky;top:22px;overflow:auto}}.sidebar h3{{margin:0 0 12px;color:#263238;font-size:1rem;border-bottom:2px solid #00695c;padding-bottom:8px}}.sidebar a{{display:block;text-decoration:none;color:#37474f;background:#f7faf8;border-left:4px solid transparent;border-radius:5px;padding:8px 10px;margin:5px 0;font-size:.88rem}}.sidebar a:hover{{border-left-color:#00695c;background:#eef6f2}}main{{flex:1;background:#fff;border:1px solid #dfe6e2;border-radius:8px;overflow:hidden}}.hero{{background:linear-gradient(135deg,#004d40,#263238);color:#fff;padding:32px 38px}}.hero a{{color:#c8e6c9;text-decoration:none}}.hero h1{{font-size:2rem;margin:8px 0 8px}}.hero p{{margin:0;max-width:980px;line-height:1.6;color:#e0f2f1}}.content{{padding:34px 42px}}h1{{font-size:1.9rem;margin:0 0 16px}}h2{{font-size:1.35rem;color:#263238;margin:32px 0 14px;padding-bottom:8px;border-bottom:2px solid #e0e6e3}}h3{{font-size:1.08rem;color:#37474f;margin:22px 0 10px}}h4{{font-size:1rem;color:#37474f;margin:18px 0 8px}}p,li{{line-height:1.68;color:#4f5b62}}ul,ol{{padding-left:24px}}code{{background:#eef3f1;border:1px solid #dbe5e0;border-radius:4px;padding:1px 5px;font-family:Consolas,'Courier New',monospace;font-size:.9em}}pre{{background:#172124;color:#d9ece5;border-radius:8px;padding:16px 18px;overflow:auto;line-height:1.55}}pre code{{background:transparent;border:0;color:inherit;padding:0}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin:22px 0}}.card{{background:#f7faf8;border:1px solid #dfe6e2;border-radius:8px;padding:16px}}.card strong{{display:block;font-size:1.65rem;color:#00695c;margin-bottom:4px}}.callout{{background:#e0f2f1;border-left:5px solid #00695c;border-radius:8px;padding:16px 18px;margin:20px 0;color:#004d40}}.figure{{margin:28px 0;background:#f7faf8;border:1px solid #dfe6e2;border-radius:8px;padding:14px}}.figure img{{display:block;width:100%;height:auto;border-radius:6px}}.caption{{font-size:.86rem;color:#607d6f;margin:8px 4px 0}}table{{border-collapse:collapse;width:100%;margin:16px 0;font-size:.9rem}}th,td{{border-bottom:1px solid #e3e8e5;padding:9px 10px;text-align:left}}th{{background:#eef3f1;color:#263238}}@media(max-width:900px){{.page{{display:block;padding:12px}}.sidebar{{position:relative;width:auto;height:auto;margin-bottom:14px}}.content{{padding:24px 20px}}}}
</style>
</head>
<body>
<div class="page">
<aside class="sidebar"><h3>Report Links</h3>{sidebar}</aside>
<main>
<div class="hero"><a href="../index.html">Back to Home</a><h1>Admixture-Aware Relatedness</h1><p>Cross-method assessment of PI_HAT, KING, Step 11 ADMIXTURE distances, PC-AiR/PC-Relate, and RelateAdmix for the merged ALSU + GWAS2026 Uzbek GSA datasets.</p></div>
<div class="content">
<div class="cards">
  <div class="card"><strong>6,228</strong>PI_HAT &ge; 0.185 pairs</div>
  <div class="card"><strong>4,494</strong>PI_HAT-only discordant pairs</div>
  <div class="card"><strong>5</strong>PI_HAT-only pairs supported by PC-Relate</div>
    <div class="card"><strong>1</strong>PI_HAT-only pair supported by RelateAdmix</div>
  <div class="card"><strong>992</strong>Step 11 samples mapped by birthplace</div>
</div>
<div class="callout"><strong>Current recommendation:</strong> keep KING as the conservative operational screen, use PC-Relate and RelateAdmix as ancestry-aware validation layers, and do not remove samples solely on PI_HAT-only calls in this admixed merged dataset.</div>
<div class="figure"><img src="../images/alsu_uzbekistan_k5_birthplace_map.png" alt="ALSU K=5 ADMIXTURE by Uzbekistan birthplace region"><div class="caption">Birthplace-region map generated from the verified Step 11 K=5 ADMIXTURE source files and the GWAS survey geography fields.</div></div>
{body}
</div>
</main>
</div>
</body>
</html>
"""
    OUT_PATH.write_text(html_text, encoding="utf-8")
    print(f"Wrote {OUT_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()