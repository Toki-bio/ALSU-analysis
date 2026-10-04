import json, re, html
D = "C:/work/alsu/nygc_site_2026-10-04/"
S = "C:/work/ALSU-analysis/"
R = json.load(open(D + "results.json"))
STYLE = """<style>
.nygc{margin:18px 0}
.nygc h2{color:#1b5e20;font-size:1.35em;margin:26px 0 10px;padding-bottom:6px;border-bottom:3px solid #16a34a}
.nygc h3{color:#333;font-size:1.1em;margin:18px 0 8px}
.nygc p,.nygc li{color:#444;line-height:1.65;margin-bottom:9px;font-size:.96em}
.nygc ul,.nygc ol{margin-left:22px}
.nygc table{width:100%;border-collapse:collapse;font-size:.86em;margin:12px 0}
.nygc th{background:#f0f0f0;padding:8px 10px;text-align:left;border-bottom:2px solid #ddd;color:#333}
.nygc td{padding:7px 10px;border-bottom:1px solid #eee;color:#444}
.nygc tbody tr:nth-child(even){background:#fafafa}
.nygc .cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:16px 0}
.nygc .card{background:#f9f9f9;border:1px solid #e8e8e8;border-radius:8px;padding:12px;text-align:center}
.nygc .card .num{font-size:1.7em;font-weight:700;color:#1b5e20}
.nygc .card .lab{font-size:.72em;color:#777;text-transform:uppercase;letter-spacing:.4px}
.nygc .box{border-radius:8px;padding:13px 17px;margin:14px 0;font-size:.92em;line-height:1.6}
.nygc .key{background:#e8f5e9;border-left:4px solid #2e7d32;color:#1b5e20}
.nygc .info{background:#e3f2fd;border-left:4px solid #1565c0;color:#0d47a1}
.nygc .warn{background:#fff3e0;border-left:4px solid #ef6c00;color:#7a3b00}
.nygc .meth{background:#f3e5f5;border-left:4px solid #7b1fa2;color:#4a148c}
.nygc img{max-width:100%;height:auto;display:block;margin:12px auto;border:1px solid #e5e5e5;border-radius:6px}
.nygc .cap{font-size:.84em;color:#666;text-align:center;margin-top:-4px}
.nygc code{background:#f3f3f3;padding:1px 5px;border-radius:3px;font-size:.88em}
.nygc pre{background:#1e1e1e;color:#d4d4d4;padding:14px 16px;border-radius:8px;overflow-x:auto;font-size:.82em;line-height:1.5;margin:12px 0}
details.nygc-old{margin:26px 0;border:1px solid #d9d9d9;border-radius:8px;padding:6px 14px;background:#fcfcfc}
details.nygc-old>summary{cursor:pointer;font-weight:600;color:#7a3b00;padding:8px 0}
</style>"""
def f(x, d=4): return f"{x:.{d}f}"
def esc(s): return html.escape(str(s))
def cards(items): return '<div class="cards">' + "".join(f'<div class="card"><div class="num">{n}</div><div class="lab">{l}</div></div>' for n, l in items) + "</div>"
def box(kind, inner): return f'<div class="box {kind}">{inner}</div>'
def table(head, rows): return "<table><thead><tr>" + "".join(f"<th>{h}</th>" for h in head) + "</tr></thead><tbody>" + "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows) + "</tbody></table>"
def img(src, cap): return f'<img src="../images/{src}" alt="{esc(cap)}"><p class="cap">{cap}</p>'
def read(p): return open(S + p, encoding="utf-8", newline="").read()
def write(p, t): open(S + p, "w", encoding="utf-8", newline="").write(t)
def wrap(block, bid): return f"<!--NYGC:{bid}:BEGIN-->{STYLE}<div class=\"nygc\">{block}</div><!--NYGC:{bid}:END-->"
def put(t, bid, block, after=None, before=None, replace_range=None):
    """idempotently place wrapped block. after: literal string to insert after; before: literal to insert before;
    replace_range: (start_literal, end_literal) -> everything from start up to (not incl.) end is replaced."""
    w = wrap(block, bid)
    t = re.sub(rf"<!--NYGC:{bid}:BEGIN-->.*?<!--NYGC:{bid}:END-->", "", t, flags=re.S)
    if replace_range:
        a, b = replace_range; i = t.index(a); j = t.index(b, i)
        return t[:i] + w + t[j:]
    if after: i = t.index(after) + len(after); return t[:i] + w + t[i:]
    if before: i = t.index(before); return t[:i] + w + t[i:]
    raise ValueError
def old_wrap(t, bid, start_lit, end_lit, title, end_regex=None):
    if f"NYGC:{bid}:OLDBEGIN" in t: return t
    i = t.index(start_lit)
    if end_regex:
        m = re.compile(end_regex, re.S).search(t, i); j = m.start()
    else:
        j = t.index(end_lit, i)
    return t[:i] + f'<!--NYGC:{bid}:OLDBEGIN--><details class="nygc-old"><summary>{title}</summary>' + t[i:j] + f"</details><!--NYGC:{bid}:OLDEND-->" + t[j:]

def banner(t, html_inner, kind="ok"):
    """replace the first correction-banner div (page-level) with an 'updated' banner; if none, insert after <body...>"""
    col = {"ok": ("#ecfdf5", "#16a34a"), "warn": ("#fff4e5", "#d97706")}[kind]
    new = f'<div class="correction-banner" style="background:{col[0]};border:2px solid {col[1]};padding:12px 16px;margin:12px 0;font-size:14px;line-height:1.45;color:#1f2937">{html_inner}</div>'
    m = re.search(r'<div class="correction-banner"[^>]*>.*?</div>', t, flags=re.S)
    if m:
        return t[:m.start()] + new + t[m.end():]
    m = re.search(r"<body[^>]*>", t)
    return t[:m.end()] + new + t[m.end():]

def wsub(t, old, new, count=1):
    toks = re.split(r"\s+", old.strip())
    pat = r"\s+".join(re.escape(x) for x in toks)
    m = re.search(pat, t)
    assert m, "wsub: not found: " + old[:60]
    return t[:m.start()] + new + t[m.end():]
