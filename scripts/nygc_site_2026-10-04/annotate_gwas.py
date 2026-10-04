import urllib.request, urllib.error, json, time
D = "C:/work/alsu/nygc_site_2026-10-04/"
c = json.load(open(D + "annot_cache.json"))
def q(rs):
    u = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rs}/associations?projection=associationBySnp"
    for k in range(8):
        try:
            a = json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "alsu", "Accept": "application/json"}), timeout=40))
            ev = (a.get("_embedded") or {}).get("associations", [])
            return "ok", [t.get("trait") for e in ev for t in (e.get("efoTraits") or [])], len(ev)
        except urllib.error.HTTPError as e:
            if e.code == 404: return "ok", [], 0          # SNP known to catalog API but no associations / not in catalog
            time.sleep(5 * (k + 1))
        except Exception:
            time.sleep(5 * (k + 1))
    return "error", [], 0
# positive control first
st, tr, n = q("rs1800414")
print("control rs1800414:", st, n, sorted(set(tr))[:4])
assert st == "ok" and n > 0, "control failed; do not trust zeros"
for s, v in c.items():
    if v.get("gwas_status") == "ok": continue
    tr_all = []; n_all = 0; stat = "ok"
    for rs in v["rsids"][:2]:
        st, tr, n = q(rs); time.sleep(1.2)
        if st == "error": stat = "error"
        tr_all += tr; n_all += n
    v["gwas_status"] = stat; v["gwas_n"] = n_all; v["gwas_traits"] = sorted(set(x for x in tr_all if x))[:8]
    json.dump(c, open(D + "annot_cache.json", "w"))
print("ok:", sum(1 for v in c.values() if v["gwas_status"] == "ok"), "error:", sum(1 for v in c.values() if v["gwas_status"] == "error"), "with hits:", sum(1 for v in c.values() if v["gwas_n"]))
for s, v in c.items():
    if v["gwas_n"]: print(s, v["rsids"], v["gwas_n"], v["gwas_traits"][:4])
