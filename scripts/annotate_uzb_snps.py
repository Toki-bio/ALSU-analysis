#!/usr/bin/env python3
"""
Annotate 490 Uzbek-specific PBS SNPs using:
  1. Ensembl REST API  → rsID, gene, consequence
  2. GWAS Catalog API  → disease associations
  3. NCBI ClinVar      → clinical significance
  4. GTEx API          → eQTL signals (Whole Blood, Artery, Liver)

Output: uzbek_snp_annotations.json + uzbek_snp_annotations.tsv

Usage: python3 annotate_uzb_snps.py
"""

import json
import time
import csv
import sys
import os
import requests
from pathlib import Path

# ─── Config ────────────────────────────────────────────────────────────────────
INPUT_TSV  = "/staging/ALSU-analysis/admixture_analysis/pop_specific/uzbek_specific_snps.tsv"
OUT_DIR    = "/staging/ALSU-analysis/admixture_analysis/pop_specific/annotation"
BUILD      = "GRCh38"   # coordinates in the input file

ENSEMBL_BASE  = "https://rest.ensembl.org"
GWAS_BASE     = "https://www.ebi.ac.uk/gwas/rest/api"
NCBI_BASE     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GTEX_BASE     = "https://gtexportal.org/api/v2"

GTEX_TISSUES = [
    "Whole_Blood",
    "Liver",
    "Artery_Aorta",
    "Brain_Frontal_Cortex_Ba9",
    "Adipose_Subcutaneous",
]

HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

RATE_ENSEMBL = 0.07   # ~15 req/s
RATE_GWAS    = 0.3
RATE_NCBI    = 0.35   # 3 req/s without API key
RATE_GTEX    = 0.2

Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
CHECKPOINT = os.path.join(OUT_DIR, "checkpoint.json")

# ─── Helpers ───────────────────────────────────────────────────────────────────
def get(url, params=None, retries=4, wait=1.5):
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=20)
            if r.status_code == 429:
                time.sleep(wait * (attempt + 2))
                continue
            if r.status_code == 200:
                return r.json()
            if r.status_code in (400, 404, 500):
                return None
        except Exception:
            time.sleep(wait)
    return None


def post_json(url, payload, retries=4, wait=1.5):
    for attempt in range(retries):
        try:
            r = requests.post(url, json=payload, headers=HEADERS, timeout=30)
            if r.status_code == 429:
                time.sleep(wait * (attempt + 2))
                continue
            if r.status_code in (200, 201):
                return r.json()
        except Exception:
            time.sleep(wait)
    return None


# ─── 1. Load SNPs ──────────────────────────────────────────────────────────────
def f(val, default=0.0):
    """Safe float conversion."""
    try:
        return float(str(val).strip())
    except (ValueError, TypeError):
        return default


def load_snps(path):
    snps = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # Normalize headers (strip whitespace)
        reader.fieldnames = [h.strip() for h in reader.fieldnames]
        for row in reader:
            snp_id = row.get("SNP", "").strip()
            if not snp_id or ":" not in snp_id:
                continue
            parts  = snp_id.split(":")
            chrom  = parts[0]
            pos    = int(parts[1])
            maf_eur = f(row.get("MAF_EUR"))
            maf_eas = f(row.get("MAF_EAS"))
            maf_sas = f(row.get("MAF_SAS"))
            maf_afr = f(row.get("MAF_AFR"))
            maf_uzb = f(row.get("MAF_UZB"))
            snps.append({
                "id":       snp_id,
                "chrom":    chrom,
                "pos":      pos,
                "A1":       row.get("A1", "").strip(),
                "A2":       row.get("A2", "").strip(),
                "MAF_UZB":  maf_uzb,
                "MAF_EUR":  maf_eur,
                "MAF_EAS":  maf_eas,
                "MAF_SAS":  maf_sas,
                "MAF_AFR":  maf_afr,
                "PBS_UZB":  f(row.get("PBS_UZB")),
                "PBS_EUR":  f(row.get("PBS_EUR")),
                "PBS_EAS":  f(row.get("PBS_EAS")),
                "near_private": (maf_eur <= 0.01 and maf_eas <= 0.01 and
                                 maf_sas <= 0.01 and maf_afr <= 0.01 and
                                 maf_uzb >= 0.05),
            })
    return snps


# ─── 2. Ensembl: coordinates → rsID + gene + VEP consequence ──────────────────
def ensembl_overlap(chrom, pos):
    """Return list of variants overlapping a single position."""
    url = f"{ENSEMBL_BASE}/overlap/region/human/{chrom}:{pos}-{pos}"
    params = {"feature": "variation", "content-type": "application/json"}
    return get(url, params=params) or []


def ensembl_vep_region(chrom, pos, allele):
    """VEP annotation for a single position."""
    url = f"{ENSEMBL_BASE}/vep/human/region/{chrom}:{pos}-{pos}:{allele}"
    return get(url) or []


def annotate_ensembl_batch(snps, results):
    """POST batch VEP for up to 200 variants at a time."""
    chunk_size = 100
    for i in range(0, len(snps), chunk_size):
        chunk   = snps[i:i + chunk_size]
        need    = [s for s in chunk if "rsid" not in results.get(s["id"], {})]
        if not need:
            continue

        variants = []
        for s in need:
            allele = s["A1"] if s["A1"] else "A"
            variants.append(f"{s['chrom']} {s['pos']} . {s['A2'] or 'N'} {allele} . . .")

        payload = {"variants": variants}
        data    = post_json(f"{ENSEMBL_BASE}/vep/human/region", payload)
        time.sleep(RATE_ENSEMBL)

        if not data:
            continue

        for item in data:
            # input field is like "9 104189856 . G C . . ."
            parts = item.get("input", "").split()
            if len(parts) >= 2:
                key = f"{parts[0]}:{parts[1]}"
                rsid = item.get("id", "")
                gene_sym = ""
                consequence = ""
                tc = item.get("transcript_consequences", [])
                if tc:
                    consequence = tc[0].get("consequence_terms", [""])[0]
                    gene_sym    = tc[0].get("gene_symbol", "")
                    biotype     = tc[0].get("biotype", "")
                    sift        = tc[0].get("sift_prediction", "")
                    polyphen    = tc[0].get("polyphen_prediction", "")
                else:
                    biotype  = ""
                    sift     = ""
                    polyphen = ""

                ig = item.get("intergenic_consequences", [])
                if ig and not gene_sym:
                    consequence = ig[0].get("consequence_terms", [""])[0]

                results.setdefault(key, {}).update({
                    "rsid":        rsid,
                    "gene":        gene_sym,
                    "consequence": consequence,
                    "sift":        sift,
                    "polyphen":    polyphen,
                })

    return results


def fallback_overlap(snp, results):
    """Single-position overlap lookup if VEP batch missed this SNP."""
    key  = snp["id"]
    data = ensembl_overlap(snp["chrom"], snp["pos"])
    time.sleep(RATE_ENSEMBL)
    if data:
        # prefer rs* IDs
        rs_candidates = [v["id"] for v in data if str(v.get("id","")).startswith("rs")]
        rsid = rs_candidates[0] if rs_candidates else (data[0].get("id","") if data else "")
        results.setdefault(key, {})["rsid"] = rsid
    return results


# ─── 3. GWAS Catalog ───────────────────────────────────────────────────────────
def gwas_lookup(rsid):
    """Return list of (trait, p-value, pubmed) for a rsID."""
    if not rsid or not rsid.startswith("rs"):
        return []
    url    = f"{GWAS_BASE}/singleNucleotidePolymorphisms/{rsid}/associations"
    data   = get(url)
    time.sleep(RATE_GWAS)
    if not data:
        return []
    assocs = data.get("_embedded", {}).get("associations", [])
    out    = []
    seen   = set()
    for a in assocs:
        traits = [t.get("trait", "") for t in a.get("efoTraits", [])]
        trait  = "; ".join(traits) if traits else a.get("description", "")
        pval   = a.get("pvalue", "")
        pmid   = ""
        study  = a.get("study", {})
        if study:
            pmids = study.get("publicationInfo", {}).get("pubmedId", "")
            pmid  = str(pmids)
        key = (trait, pval)
        if key not in seen:
            seen.add(key)
            out.append({"trait": trait, "pvalue": pval, "pmid": pmid})
    return out[:5]  # top 5


# ─── 4. ClinVar ────────────────────────────────────────────────────────────────
def clinvar_lookup(rsid):
    """Return clinical significance string and condition."""
    if not rsid or not rsid.startswith("rs"):
        return {}
    search_url = f"{NCBI_BASE}/esearch.fcgi"
    params     = {"db": "clinvar", "term": f"{rsid}[rs]", "retmode": "json", "retmax": 5}
    data       = get(search_url, params=params)
    time.sleep(RATE_NCBI)
    if not data:
        return {}
    ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return {}

    # Fetch summary for first ID
    sum_url  = f"{NCBI_BASE}/esummary.fcgi"
    s_params = {"db": "clinvar", "id": ids[0], "retmode": "json"}
    sdata    = get(sum_url, params=s_params)
    time.sleep(RATE_NCBI)
    if not sdata:
        return {}

    result = sdata.get("result", {})
    uid    = ids[0]
    entry  = result.get(uid, {})
    clinsig = entry.get("clinical_significance", {})
    sig     = clinsig.get("description", "") if isinstance(clinsig, dict) else str(clinsig)
    conds   = entry.get("trait_set", [])
    cond_names = []
    for c in (conds if isinstance(conds, list) else []):
        for n in c.get("trait_name", []):
            if isinstance(n, dict):
                cond_names.append(n.get("value", ""))
            elif isinstance(n, str):
                cond_names.append(n)
    return {
        "clinvar_sig":       sig,
        "clinvar_condition": "; ".join(cond_names[:3]),
        "clinvar_id":        uid,
    }


# ─── 5. GTEx ───────────────────────────────────────────────────────────────────
def gtex_lookup(snp):
    """Check GTEx eQTL signal. Returns list of {tissue, gene, pval, effect}."""
    # Build variant ID in GTEx format: chr9_104189856_C_G_b38
    chrom = f"chr{snp['chrom']}"
    pos   = snp["pos"]
    a2    = snp["A2"] or "N"   # ref allele
    a1    = snp["A1"] or "N"   # alt allele
    varid = f"{chrom}_{pos}_{a2}_{a1}_b38"

    out = []
    for tissue in GTEX_TISSUES:
        url    = f"{GTEX_BASE}/association/singleTissueEqtl"
        params = {"variantId": varid, "tissueSiteDetailId": tissue,
                  "datasetId": "gtex_v8"}
        data   = get(url, params=params)
        time.sleep(RATE_GTEX)
        if not data:
            continue
        hits = data.get("data", [])
        for h in hits[:2]:
            pval = h.get("pValue", None)
            if pval is not None and float(pval) < 0.05:
                out.append({
                    "tissue":    tissue,
                    "gene":      h.get("geneSymbol", ""),
                    "pvalue":    pval,
                    "effect":    h.get("nes", ""),
                })
    return out


# ─── Main ──────────────────────────────────────────────────────────────────────
def main():
    print(f"Loading SNPs from {INPUT_TSV}...")
    snps = load_snps(INPUT_TSV)
    print(f"  {len(snps)} SNPs loaded")

    # Load checkpoint
    if os.path.exists(CHECKPOINT):
        with open(CHECKPOINT) as f:
            results = json.load(f)
        print(f"  Resuming from checkpoint ({len(results)} done)")
    else:
        results = {}

    total = len(snps)

    # ── Phase 1: Ensembl VEP batch ──────────────────────────────────────────
    print("\n[1/4] Ensembl VEP batch annotation...")
    need_vep = [s for s in snps if "rsid" not in results.get(s["id"], {})]
    print(f"  {len(need_vep)} SNPs need Ensembl annotation")
    annotate_ensembl_batch(need_vep, results)

    # Fallback for any still missing rsID
    still_missing = [s for s in snps if not results.get(s["id"], {}).get("rsid")]
    print(f"  Fallback overlap for {len(still_missing)} SNPs...")
    for i, s in enumerate(still_missing):
        fallback_overlap(s, results)
        if i % 50 == 0:
            print(f"    {i}/{len(still_missing)}")

    with open(CHECKPOINT, "w") as f:
        json.dump(results, f)
    print("  Ensembl done. Checkpoint saved.")

    # ── Phase 2: GWAS Catalog ───────────────────────────────────────────────
    print("\n[2/4] GWAS Catalog lookup...")
    need_gwas = [s for s in snps if "gwas" not in results.get(s["id"], {})]
    print(f"  {len(need_gwas)} SNPs to query")
    for i, s in enumerate(need_gwas):
        rsid = results.get(s["id"], {}).get("rsid", "")
        hits = gwas_lookup(rsid)
        results.setdefault(s["id"], {})["gwas"] = hits
        if i % 50 == 0:
            print(f"    {i}/{len(need_gwas)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(results, f)
    with open(CHECKPOINT, "w") as f:
        json.dump(results, f)
    print("  GWAS done. Checkpoint saved.")

    # ── Phase 3: ClinVar ────────────────────────────────────────────────────
    print("\n[3/4] ClinVar lookup...")
    need_cv = [s for s in snps if "clinvar_sig" not in results.get(s["id"], {})]
    print(f"  {len(need_cv)} SNPs to query")
    for i, s in enumerate(need_cv):
        rsid = results.get(s["id"], {}).get("rsid", "")
        cv   = clinvar_lookup(rsid)
        results.setdefault(s["id"], {}).update(cv)
        if i % 50 == 0:
            print(f"    {i}/{len(need_cv)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(results, f)
    with open(CHECKPOINT, "w") as f:
        json.dump(results, f)
    print("  ClinVar done. Checkpoint saved.")

    # ── Phase 4: GTEx ───────────────────────────────────────────────────────
    print("\n[4/4] GTEx eQTL lookup...")
    need_gtex = [s for s in snps if "gtex" not in results.get(s["id"], {})]
    print(f"  {len(need_gtex)} SNPs to query")
    for i, s in enumerate(need_gtex):
        hits = gtex_lookup(s)
        results.setdefault(s["id"], {})["gtex"] = hits
        if i % 25 == 0:
            print(f"    {i}/{len(need_gtex)}")
            with open(CHECKPOINT, "w") as f:
                json.dump(results, f)
    with open(CHECKPOINT, "w") as f:
        json.dump(results, f)
    print("  GTEx done. Checkpoint saved.")

    # ── Merge & output ─────────────────────────────────────────────────────
    print("\nMerging results...")
    output = []
    for s in snps:
        ann = results.get(s["id"], {})
        row = {
            "snp_id":       s["id"],
            "chrom":        s["chrom"],
            "pos":          s["pos"],
            "A1":           s["A1"],
            "A2":           s["A2"],
            "PBS_UZB":      round(s["PBS_UZB"], 4),
            "MAF_UZB":      round(s["MAF_UZB"], 4),
            "MAF_EUR":      round(s["MAF_EUR"], 4),
            "MAF_EAS":      round(s["MAF_EAS"], 4),
            "near_private": s["near_private"],
            "rsid":         ann.get("rsid", ""),
            "gene":         ann.get("gene", ""),
            "consequence":  ann.get("consequence", ""),
            "sift":         ann.get("sift", ""),
            "polyphen":     ann.get("polyphen", ""),
            "gwas_hits":    len(ann.get("gwas", [])),
            "gwas_top":     ann.get("gwas", [{}])[0].get("trait", "") if ann.get("gwas") else "",
            "clinvar_sig":  ann.get("clinvar_sig", ""),
            "clinvar_cond": ann.get("clinvar_condition", ""),
            "gtex_hits":    len(ann.get("gtex", [])),
            "gtex_top":     ann.get("gtex", [{}])[0].get("tissue", "") if ann.get("gtex") else "",
            "gwas_detail":  ann.get("gwas", []),
            "gtex_detail":  ann.get("gtex", []),
        }
        output.append(row)

    # Sort by PBS_UZB descending
    output.sort(key=lambda x: x["PBS_UZB"], reverse=True)

    # Save JSON
    json_path = os.path.join(OUT_DIR, "uzbek_snp_annotations.json")
    with open(json_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"  JSON saved: {json_path}")

    # Save TSV
    tsv_path = os.path.join(OUT_DIR, "uzbek_snp_annotations.tsv")
    tsv_cols  = ["snp_id","chrom","pos","A1","A2","PBS_UZB","MAF_UZB","MAF_EUR",
                 "MAF_EAS","near_private","rsid","gene","consequence","sift",
                 "polyphen","gwas_hits","gwas_top","clinvar_sig","clinvar_cond",
                 "gtex_hits","gtex_top"]
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=tsv_cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output)
    print(f"  TSV saved: {tsv_path}")

    # Summary stats
    gwas_count    = sum(1 for r in output if r["gwas_hits"] > 0)
    clinvar_count = sum(1 for r in output if r["clinvar_sig"])
    gtex_count    = sum(1 for r in output if r["gtex_hits"] > 0)
    rsid_count    = sum(1 for r in output if r["rsid"].startswith("rs"))

    print(f"""
Summary:
  Total SNPs:           {total}
  With rsID:            {rsid_count}
  GWAS catalog hits:    {gwas_count}
  ClinVar entries:      {clinvar_count}
  GTEx eQTL hits:       {gtex_count}
""")

    # Quick summary for HTML embedding
    summary = {
        "total": total, "rsid_count": rsid_count,
        "gwas_count": gwas_count, "clinvar_count": clinvar_count,
        "gtex_count": gtex_count,
        "top20": output[:20],
    }
    with open(os.path.join(OUT_DIR, "annotation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("Done!")


if __name__ == "__main__":
    main()
