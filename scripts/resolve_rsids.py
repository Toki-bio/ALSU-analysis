import requests, time, json

SNPS = [
    ("12", 22967890, "G", "T"),
    ("12", 125520190, "C", "T"),
    ("12", 5664803, "T", "G"),
    ("11", 20665570, "T", "G"),
    ("5", 53879140, "C", "A"),
    ("3", 133749168, "A", "G"),
    ("11", 207698, "C", "T"),
    ("10", 8512594, "G", "A"),
]
H = {"Content-Type": "application/json"}
rsid_map = {}
for c, p, ref, alt in SNPS:
    r = requests.get(f"https://rest.ensembl.org/overlap/region/human/{c}:{p}-{p}?feature=variation", headers=H, timeout=15)
    if r.status_code == 200:
        found = False
        for v in r.json():
            vid = v.get("id", "")
            if vid.startswith("rs"):
                rsid_map[f"{c}:{p}"] = vid
                print(f"{c}:{p}:{ref}:{alt} -> {vid}")
                found = True
                break
        if not found:
            print(f"{c}:{p}:{ref}:{alt} -> no rsID in overlap")
    else:
        print(f"{c}:{p}:{ref}:{alt} -> HTTP {r.status_code}")
    time.sleep(0.1)

with open("data/rsid_map.json", "w") as f:
    json.dump(rsid_map, f, indent=2)
print(f"\nResolved {len(rsid_map)}/{len(SNPS)} rsIDs")
