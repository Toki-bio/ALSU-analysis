import urllib.request, json, time, os
D = "C:/work/alsu/nygc_site_2026-10-04/"
R = json.load(open(D + "results.json"))
cp = D + "annot_cache.json"
cache = json.load(open(cp)) if os.path.exists(cp) else {}


def get(u):
    for k in range(4):
        try:
            r = urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "alsu", "Accept": "application/json"}), timeout=40)
            return json.load(r)
        except Exception as e:
            err = e
            time.sleep(2)
    return {"_error": str(err)}


def annot(snp):
    ch, pos = snp.split(":")[:2]
    v = get(f"https://rest.ensembl.org/overlap/region/human/{ch}:{pos}-{pos}?feature=variation;content-type=application/json")
    out = {"rsids": [], "cons": [], "clin": []}
    if isinstance(v, list):
        for x in v:
            if x.get("start") == int(pos):
                out["rsids"].append(x["id"]); out["cons"].append(x.get("consequence_type")); out["clin"] += x.get("clinical_significance") or []
    else:
        out["error"] = str(v)
    gw = []
    for rs in out["rsids"][:2]:
        a = get(f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rs}/associations?projection=associationBySnp")
        for e in (a.get("_embedded", {}) or {}).get("associations", []) if isinstance(a, dict) else []:
            for t in e.get("efoTraits", []) or []:
                gw.append(t.get("trait"))
        time.sleep(.3)
    out["gwas_n"] = len(gw)
    out["gwas_traits"] = sorted(set(x for x in gw if x))[:8]
    return out


snps = [r["snp"] for r in R["ue_top30"]] + [r["snp"] for r in R["pbs_top15"]]
for s in snps:
    if s not in cache:
        cache[s] = annot(s)
        time.sleep(.2)
        json.dump(cache, open(cp, "w"))
print(len(cache))
for s in snps[:6]:
    print(s, cache[s])
