#!/usr/bin/env python3
"""Parse dl_all.txt: decode base64 blocks to files, extract text sections."""
import base64, re, os

with open("dl_all.txt", encoding="utf-8") as f:
    raw = f.read()

# Strip SSH banner (lines before first ===)
start = raw.find("=== UZB K6-8 Q BASE64 ===")
if start < 0:
    raise ValueError("Marker not found")
raw = raw[start:]

# Split into sections by === ... ===
sections = re.split(r'^(=== .+? ===)\s*$', raw, flags=re.MULTILINE)
# sections[0] is empty, then alternating: marker, content, marker, content...

data = {}
for i in range(1, len(sections), 2):
    key = sections[i].strip("= \n")
    val = sections[i+1] if i+1 < len(sections) else ""
    data[key] = val.strip()

# --- Decode UZB K=6,7,8 ---
for k in [6, 7, 8]:
    marker = f"--- UZB_K={k} ---"
    block = data.get("UZB K6-8 Q BASE64", "")
    # Split by sub-markers
    parts = re.split(r'^--- UZB_K=\d+ ---\s*$', block, flags=re.MULTILINE)
    sub_markers = re.findall(r'^--- UZB_K=(\d+) ---\s*$', block, flags=re.MULTILINE)
    kmap = {}
    for j, km in enumerate(sub_markers):
        kmap[int(km)] = parts[j+1].strip() if j+1 < len(parts) else ""
    
    if k in kmap:
        decoded = base64.b64decode(kmap[k])
        outpath = f"data/UZB_v2_admix.{k}.Q"
        with open(outpath, "wb") as out:
            out.write(decoded)
        lines = decoded.decode().strip().split('\n')
        print(f"UZB K={k}: {len(lines)} lines, {len(decoded)} bytes -> {outpath}")

# --- Decode Global K=2..8 ---
gblock = data.get("GLOBAL Q BASE64", "")
gsub = re.findall(r'^--- GLOBAL_K=(\d+) ---\s*$', gblock, flags=re.MULTILINE)
gparts = re.split(r'^--- GLOBAL_K=\d+ ---\s*$', gblock, flags=re.MULTILINE)
os.makedirs("data/global", exist_ok=True)
for j, km in enumerate(gsub):
    k = int(km)
    b64 = gparts[j+1].strip() if j+1 < len(gparts) else ""
    decoded = base64.b64decode(b64)
    outpath = f"data/global/global_v2_admix.{k}.Q"
    with open(outpath, "wb") as out:
        out.write(decoded)
    lines = decoded.decode().strip().split('\n')
    print(f"Global K={k}: {len(lines)} lines, {len(decoded)} bytes -> {outpath}")

# --- Decode Global FAM ---
fam_b64 = data.get("GLOBAL FAM", "")
if fam_b64:
    decoded = base64.b64decode(fam_b64)
    outpath = "data/global/global_v2_admix.fam"
    with open(outpath, "wb") as out:
        out.write(decoded)
    lines = decoded.decode().strip().split('\n')
    print(f"Global FAM: {len(lines)} lines -> {outpath}")

# --- sNMF CE ---
ce = data.get("sNMF CE", "")
if ce:
    with open("data/snmf_v2_ce.csv", "w") as out:
        out.write(ce + "\n")
    print(f"sNMF CE saved")
    for line in ce.strip().split('\n'):
        print(f"  {line}")

# --- Global log-likelihoods ---
ll = data.get("GLOBAL LOG-LIKELIHOODS", "")
if ll:
    with open("data/global_v2_loglik.txt", "w") as out:
        out.write(ll + "\n")
    print(f"Global log-likelihoods saved")
    for line in ll.strip().split('\n'):
        print(f"  {line}")

print("\nDone!")
