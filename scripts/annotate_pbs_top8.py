#!/usr/bin/env python3
"""Annotate the 8 PBS candidate SNPs via Ensembl VEP, GWAS Catalog, ClinVar, GTEx."""
import json, time, sys, xml.etree.ElementTree as ET
import requests

HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
OUT = "data/pbs_annotations.json"

# 8 PBS candidates from step12
SNPS = [
    {"id": "12:22967890:G:T",  "chr": "12", "pos": 22967890,  "ref": "G", "alt": "T", "pbs": 2.9885},
    {"id": "12:125520190:C:T", "chr": "12", "pos": 125520190, "ref": "C", "alt": "T", "pbs": 2.6867},
    {"id": "12:5664803:T:G",   "chr": "12", "pos": 5664803,   "ref": "T", "alt": "G", "pbs": 2.6689},
    {"id": "11:20665570:T:G",  "chr": "11", "pos": 20665570,  "ref": "T", "alt": "G", "pbs": 1.1429},
    {"id": "5:53879140:C:A",   "chr": "5",  "pos": 53879140,  "ref": "C", "alt": "A", "pbs": 0.5306},
    {"id": "3:133749168:A:G",  "chr": "3",  "pos": 133749168, "ref": "A", "alt": "G", "pbs": 0.4878},
    {"id": "11:207698:C:T",    "chr": "11", "pos": 207698,    "ref": "C", "alt": "T", "pbs": 0.3799},
    {"id": "10:8512594:G:A",   "chr": "10", "pos": 8512594,   "ref": "G", "alt": "A", "pbs": 0.3229},
]

GTEX_TISSUES = ["Whole_Blood", "Liver", "Artery_Aorta", "Brain_Frontal_Cortex_BA9", "Adipose_Subcutaneous"]

def resolve_rsids(snps):
    """Use Ensembl overlap API to get proper rsIDs."""
    print("Phase 0: Resolving rsIDs via Ensembl overlap...")
    for s in snps:
        r = requests.get(
            f"https://rest.ensembl.org/overlap/region/human/{s['chr']}:{s['pos']}-{s['pos']}?feature=variation",
            headers=HEADERS, timeout=15)
        if r.status_code == 200:
            for v in r.json():
                vid = v.get("id", "")
                if vid.startswith("rs"):
                    s["rsid"] = vid
                    print(f"  {s['id']} -> {vid}")
                    break
        time.sleep(0.1)

def get(url, params=None, retries=3):
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=20)
            if r.status_code == 429:
                time.sleep(2 * (attempt + 1)); continue
            if r.status_code == 200:
                return r.json()
            if r.status_code in (400, 404):
                return None
        except Exception:
            time.sleep(1)
    return None

def post(url, payload, retries=3):
    for attempt in range(retries):
        try:
            r = requests.post(url, json=payload, headers=HEADERS, timeout=30)
            if r.status_code == 429:
                time.sleep(2 * (attempt + 1)); continue
            if r.status_code in (200, 201):
                return r.json()
        except Exception:
            time.sleep(1)
    return None

# --- 1. Ensembl VEP batch ---
def run_vep(snps):
    print("Phase 1: Ensembl VEP batch...")
    variants = []
    for s in snps:
        variants.append(f"{s['chr']} {s['pos']} {s['id']} {s['ref']} {s['alt']} . . .")
    data = post("https://rest.ensembl.org/vep/human/region", {"variants": variants})
    results = {}
    if not data:
        print("  VEP batch returned nothing, trying individual lookups")
        for s in snps:
            url = f"https://rest.ensembl.org/vep/human/region/{s['chr']}:{s['pos']}-{s['pos']}:1/{s['alt']}"
            d = get(url)
            time.sleep(0.1)
            if d and isinstance(d, list) and len(d) > 0:
                results[s['id']] = parse_vep_item(d[0])
                print(f"  {s['id']} -> {results[s['id']].get('rsid', '?')}")
        return results

    for item in data:
        inp = item.get("input", "")
        parts = inp.split()
        snp_id = parts[2] if len(parts) >= 3 else ""
        results[snp_id] = parse_vep_item(item)
        print(f"  {snp_id} -> rsid={results[snp_id].get('rsid','?')}, gene={results[snp_id].get('gene','?')}")
    return results

def parse_vep_item(item):
    rsid = item.get("id", ".")
    if rsid.startswith("."): rsid = ""
    colocated = item.get("colocated_variants", [])
    if not rsid and colocated:
        for cv in colocated:
            rid = cv.get("id", "")
            if rid.startswith("rs"):
                rsid = rid; break
    gene, consequence, biotype, sift, polyphen = "", "", "", "", ""
    tc = item.get("transcript_consequences", [])
    if tc:
        # prefer protein_coding
        best = tc[0]
        for t in tc:
            if t.get("biotype") == "protein_coding":
                best = t; break
        consequence = ", ".join(best.get("consequence_terms", []))
        gene = best.get("gene_symbol", "")
        biotype = best.get("biotype", "")
        sift = best.get("sift_prediction", "")
        polyphen = best.get("polyphen_prediction", "")
    else:
        ig = item.get("intergenic_consequences", [])
        if ig:
            consequence = ", ".join(ig[0].get("consequence_terms", []))
    return {"rsid": rsid, "gene": gene, "consequence": consequence,
            "biotype": biotype, "sift": sift, "polyphen": polyphen}

# --- 2. GWAS Catalog ---
def run_gwas(snps, vep):
    print("Phase 2: GWAS Catalog...")
    results = {}
    for s in snps:
        rsid = s.get("rsid", "") or vep.get(s["id"], {}).get("rsid", "")
        if not rsid or not rsid.startswith("rs"):
            print(f"  {s['id']}: no rsID, skipping GWAS")
            results[s["id"]] = {"hits": 0, "associations": []}
            continue
        data = get(f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations")
        time.sleep(0.3)
        assocs = []
        if data and "_embedded" in data:
            for a in data["_embedded"].get("associations", []):
                traits = [t.get("trait", "") for t in a.get("efoTraits", [])]
                pval = a.get("pvalue", "")
                assocs.append({"traits": traits, "pvalue": pval,
                               "betaNum": a.get("betaNum"), "betaUnit": a.get("betaUnit")})
        results[s["id"]] = {"hits": len(assocs), "associations": assocs[:10]}
        print(f"  {rsid}: {len(assocs)} associations")
    return results

# --- 3. ClinVar ---
def run_clinvar(snps, vep):
    print("Phase 3: ClinVar...")
    results = {}
    for s in snps:
        rsid = s.get("rsid", "") or vep.get(s["id"], {}).get("rsid", "")
        if not rsid or not rsid.startswith("rs"):
            results[s["id"]] = {"significance": "", "conditions": []}
            continue
        # esearch
        data = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                    {"db": "clinvar", "term": rsid, "retmode": "json"})
        time.sleep(0.35)
        ids = []
        if data and "esearchresult" in data:
            ids = data["esearchresult"].get("idlist", [])
        if not ids:
            results[s["id"]] = {"significance": "", "conditions": []}
            print(f"  {rsid}: not in ClinVar")
            continue
        # esummary
        sdata = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                     {"db": "clinvar", "id": ",".join(ids[:5]), "retmode": "json"})
        time.sleep(0.35)
        sig, conds = "", []
        if sdata and "result" in sdata:
            for uid in ids[:5]:
                rec = sdata["result"].get(uid, {})
                if rec:
                    sig = rec.get("clinical_significance", {}).get("description", "") if isinstance(rec.get("clinical_significance"), dict) else str(rec.get("clinical_significance", ""))
                    title = rec.get("title", "")
                    if title: conds.append(title)
        results[s["id"]] = {"significance": sig, "conditions": conds}
        print(f"  {rsid}: {sig or 'no significance'}")
    return results

# --- 4. GTEx eQTL ---
def run_gtex(snps, vep):
    print("Phase 4: GTEx eQTL...")
    results = {}
    for s in snps:
        rsid = vep.get(s["id"], {}).get("rsid", "")
        eqtls = []
        variant_id = f"chr{s['chr']}_{s['pos']}_{s['ref']}_{s['alt']}_b38"
        for tissue in GTEX_TISSUES:
            data = get("https://gtexportal.org/api/v2/association/singleTissueEqtl",
                       {"variantId": variant_id, "tissueSiteDetailId": tissue, "datasetId": "gtex_v8"})
            time.sleep(0.2)
            if data and data.get("data"):
                for hit in data["data"][:3]:
                    eqtls.append({"tissue": tissue, "gene": hit.get("geneSymbol", ""),
                                  "pvalue": hit.get("pValue"), "nes": hit.get("nes")})
        results[s["id"]] = {"hits": len(eqtls), "eqtls": eqtls}
        print(f"  {s['id']}: {len(eqtls)} eQTL hits")
    return results

# --- Main ---
def main():
    resolve_rsids(SNPS)
    vep = run_vep(SNPS)
    gwas = run_gwas(SNPS, vep)
    clinvar = run_clinvar(SNPS, vep)
    gtex = run_gtex(SNPS, vep)

    combined = []
    for s in SNPS:
        sid = s["id"]
        v = vep.get(sid, {})
        g = gwas.get(sid, {})
        c = clinvar.get(sid, {})
        e = gtex.get(sid, {})
        combined.append({
            **s,
            "rsid": s.get("rsid", "") or v.get("rsid", ""),
            "gene": v.get("gene", ""),
            "consequence": v.get("consequence", ""),
            "biotype": v.get("biotype", ""),
            "sift": v.get("sift", ""),
            "polyphen": v.get("polyphen", ""),
            "gwas_hits": g.get("hits", 0),
            "gwas_associations": g.get("associations", []),
            "clinvar_significance": c.get("significance", ""),
            "clinvar_conditions": c.get("conditions", []),
            "gtex_hits": e.get("hits", 0),
            "gtex_eqtls": e.get("eqtls", []),
        })

    with open(OUT, "w") as f:
        json.dump(combined, f, indent=2)
    print(f"\nDone. {len(combined)} SNPs annotated -> {OUT}")
    for s in combined:
        print(f"  {s['id']} rsid={s['rsid']} gene={s['gene']} cons={s['consequence']} gwas={s['gwas_hits']} clinvar={s['clinvar_significance'] or '-'} gtex={s['gtex_hits']}")

if __name__ == "__main__":
    main()
