#!/bin/bash
set -e

# =======================================================================
# V2 PBS + LD Recomputation (SELF-CONTAINED: builds pop_mapping on server)
# =======================================================================

WORKDIR=~/v2/pbs_v2
MERGED=~/v2/global_pca/UZB_1kG_v2_merged
UZB_QC=~/v2/plink/UZB_v2_qc

mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "===== V2 PBS RECOMPUTATION ====="
echo "Start: $(date)"

# -----------------------------------------------------------------------
# 1. Build pop_mapping from FAM + 1000G panel info
# -----------------------------------------------------------------------
echo "=== Step 1: Build population assignments ==="

# The merged FAM has 1000G samples (named like HG00096, NA18525 etc)
# and UZB samples (named like 2_01-02, 13_02-82 etc)
# We need: 1000G -> their superpop, all others -> UZB

# Use the pop_mapping that's already on the server at ~/v2/pca/
# But that one has fine-grained labels. So we'll build superpop from it.

# Strategy: match IID from FAM against 1000G panel on server
# If exists in 1000G -> use that superpop
# Otherwise -> UZB

# The global_admixture already had correct assignments, let's check
if [ -f ~/v2/global_admixture/pop_mapping.txt ]; then
    echo "  Found global_admixture pop_mapping"
    head -3 ~/v2/global_admixture/pop_mapping.txt
    echo "  ..."
    awk '{print $NF}' ~/v2/global_admixture/pop_mapping.txt | sort | uniq -c | sort -rn
fi

# Check what pop_mapping files exist
echo ""
echo "  Available pop_mapping files:"
find ~/v2 -name "pop_mapping*" -ls 2>/dev/null

# Check the one used for global PCA (step 8)
echo ""
echo "  global_pca pop_mapping:"
if [ -f ~/v2/global_pca/pop_mapping.txt ]; then
    head -5 ~/v2/global_pca/pop_mapping.txt
    echo "..."
    awk '{print $NF}' ~/v2/global_pca/pop_mapping.txt | sort | uniq -c | sort -rn
fi

# Build correct superpop mapping from the FAM file
# 1000G sample IDs start with letters (HG, NA); UZB start with digits
awk '{
    iid = $2
    if (iid ~ /^(HG|NA)/) {
        # This is a 1000G sample - we need their superpop
        # We will fill this from the 1000G panel
        print iid, "1KG"
    } else {
        print iid, "UZB"
    }
}' "${MERGED}.fam" > /tmp/fam_prelabel.txt

N_1KG=$(grep -c "1KG" /tmp/fam_prelabel.txt || echo 0)
N_UZB=$(grep -c "UZB" /tmp/fam_prelabel.txt || echo 0)
echo ""
echo "  FAM prelabeling: $N_1KG 1000G samples, $N_UZB UZB samples"

# Now find the 1000G panel file to get superpops
echo "  Looking for 1000G panel..."
PANEL=""
for p in ~/v2/global_pca/pop_mapping.txt ~/v2/pca/pop_mapping.txt ~/v2/global_admixture/pop_mapping.txt; do
    if [ -f "$p" ]; then
        # Check if it has EUR/EAS/SAS/AFR/AMR labels
        if grep -q "EUR" "$p"; then
            PANEL="$p"
            break
        fi
    fi
done

if [ -n "$PANEL" ]; then
    echo "  Using panel: $PANEL"
    # Build final pop_mapping: for each sample in FAM, find in panel; default UZB
    awk 'NR==FNR {
        # Read panel: could be "sample pop" or "sample\tpop" format
        gsub(/\r/,"")
        pop[$1] = $NF
        next
    }
    {
        iid = $2
        if (iid in pop && pop[iid] ~ /^(EUR|EAS|SAS|AFR|AMR)$/) {
            print $1, iid, pop[iid]
        } else {
            print $1, iid, "UZB"
        }
    }' "$PANEL" "${MERGED}.fam" > clusters.txt
else
    echo "  ERROR: No suitable pop_mapping found!"
    echo "  Will try to derive from sample naming conventions + integrated panel"
    
    # Last resort: use integrated_call_samples file if it exists
    PANEL1K=$(find /staging ~/v2 -name "integrated_call_samples*panel*" 2>/dev/null | head -1)
    if [ -n "$PANEL1K" ]; then
        echo "  Found 1000G panel: $PANEL1K"
        awk 'NR==FNR && NR>1 {
            pop[$1] = $3  # sample, pop, super_pop
            next
        }
        {
            iid = $2
            if (iid in pop) print $1, iid, pop[iid]
            else print $1, iid, "UZB"
        }' "$PANEL1K" "${MERGED}.fam" > clusters.txt
    else
        echo "  FATAL: Cannot determine population assignments"
        exit 1
    fi
fi

echo ""
echo "  Population counts:"
awk '{print $3}' clusters.txt | sort | uniq -c | sort -rn
echo ""

N_UZB_FINAL=$(awk '$3=="UZB"' clusters.txt | wc -l)
N_EUR_FINAL=$(awk '$3=="EUR"' clusters.txt | wc -l)
echo "  UZB=$N_UZB_FINAL EUR=$N_EUR_FINAL"

# Sanity check
if [ "$N_UZB_FINAL" -lt 500 ]; then
    echo "  ERROR: Too few UZB samples ($N_UZB_FINAL). Something wrong with pop assignment."
    echo "  First 10 cluster entries:"
    head -10 clusters.txt
    exit 1
fi

# -----------------------------------------------------------------------
# 2-4. Population splits, frequencies, FST
# -----------------------------------------------------------------------
echo "=== Step 2: Per-population sample lists ==="
for POP in UZB EUR EAS SAS AFR; do
    awk -v p="$POP" '$3==p {print $1, $2}' clusters.txt > keep_${POP}.txt
    echo "  $POP: $(wc -l < keep_${POP}.txt) samples"
done
echo ""

echo "=== Step 3: Allele frequencies ==="
for POP in UZB EUR EAS SAS AFR; do
    plink --bfile "$MERGED" \
          --keep keep_${POP}.txt \
          --freq \
          --out freq_${POP} \
          --allow-no-sex --silent
    echo "  $POP done"
done
echo ""

echo "=== Step 4: Pairwise per-SNP FST ==="
compute_fst() {
    POP1=$1; POP2=$2
    cat keep_${POP1}.txt keep_${POP2}.txt > keep_${POP1}_${POP2}.txt
    awk -v p="$POP1" '{print $1, $2, p}' keep_${POP1}.txt > within_${POP1}_${POP2}.txt
    awk -v p="$POP2" '{print $1, $2, p}' keep_${POP2}.txt >> within_${POP1}_${POP2}.txt
    plink --bfile "$MERGED" \
          --keep keep_${POP1}_${POP2}.txt \
          --fst \
          --within within_${POP1}_${POP2}.txt \
          --out fst_${POP1}_${POP2} \
          --allow-no-sex --silent
    echo "  $POP1 vs $POP2: $(grep 'Weighted' fst_${POP1}_${POP2}.log | awk '{print $NF}')"
}

compute_fst UZB EUR
compute_fst UZB EAS
compute_fst EUR EAS
compute_fst UZB SAS
compute_fst UZB AFR
echo ""

# -----------------------------------------------------------------------
# 5. PBS calculation (Python)
# -----------------------------------------------------------------------
echo "=== Step 5: PBS calculation ==="

python3 << 'PYEOF'
import math, json, os, statistics
import numpy as np
from collections import Counter

def read_frq(path):
    freqs = {}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p)>=5: freqs[p[1]] = float(p[4]) if p[4]!='NA' else None
    return freqs

def read_fst(path):
    fst = {}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p)>=5: fst[p[1]] = float(p[4]) if p[4] not in ('nan','NA') else 0.0
    return fst

def fst_to_T(f):
    return -math.log(1.0 - max(0.0, min(f, 0.999)))

freq = {p: read_frq(f"freq_{p}.frq") for p in ['UZB','EUR','EAS','SAS','AFR']}
fst_data = {}
for pair in ['UZB_EUR','UZB_EAS','EUR_EAS','UZB_SAS','UZB_AFR']:
    fst_data[pair] = read_fst(f"fst_{pair}.fst")

bim = {}
with open(os.path.expanduser("~/v2/global_pca/UZB_1kG_v2_merged.bim")) as f:
    for line in f:
        p = line.split(); bim[p[1]] = (p[0], int(p[3]))

common = set(fst_data['UZB_EUR']) & set(fst_data['UZB_EAS']) & set(fst_data['EUR_EAS'])
for p in ['UZB','EUR','EAS','SAS','AFR']: common &= set(freq[p])
print(f"  Common SNPs: {len(common)}")

results = []
for snp in common:
    T_ue = fst_to_T(fst_data['UZB_EUR'][snp])
    T_ua = fst_to_T(fst_data['UZB_EAS'][snp])
    T_ea = fst_to_T(fst_data['EUR_EAS'][snp])
    pbs_uzb = (T_ue + T_ua - T_ea) / 2.0
    pbs_eur = (T_ue + T_ea - T_ua) / 2.0
    pbs_eas = (T_ua + T_ea - T_ue) / 2.0
    mafs = {p: freq[p].get(snp) or 0 for p in ['UZB','EUR','EAS','SAS','AFR']}
    delta_af = min(abs(mafs['UZB']-mafs[p]) for p in ['EUR','EAS','SAS','AFR'])
    near_priv = mafs['UZB']>=0.05 and all(mafs[p]<=0.01 for p in ['EUR','EAS','SAS','AFR'])
    chrom, bp = bim.get(snp, ('?',0))
    results.append({
        'snp':snp,'chr':chrom,'bp':bp,
        'pbs_uzb':pbs_uzb,'pbs_eur':pbs_eur,'pbs_eas':pbs_eas,
        'maf_uzb':mafs['UZB'],'maf_eur':mafs['EUR'],'maf_eas':mafs['EAS'],
        'maf_sas':mafs['SAS'],'maf_afr':mafs['AFR'],
        'delta_af':delta_af,'near_private':near_priv,
        'tier1':pbs_uzb>=0.3,'tier2':delta_af>=0.3,'tier3':near_priv
    })

pbs_vals = sorted(r['pbs_uzb'] for r in results)
n = len(pbs_vals)
tier1=[r for r in results if r['tier1']]
tier2=[r for r in results if r['tier2']]
tier3=[r for r in results if r['tier3']]
cands=[r for r in results if r['tier1'] or r['tier2'] or r['tier3']]

stats = {
    'n_snps':n,'mean':statistics.mean(pbs_vals),'median':statistics.median(pbs_vals),
    'stdev':statistics.stdev(pbs_vals),'min':min(pbs_vals),'max':max(pbs_vals),
    'p95':pbs_vals[int(0.95*n)],'p99':pbs_vals[int(0.99*n)],'p999':pbs_vals[int(0.999*n)],
    'n_pbs_ge_03':len(tier1),'n_pbs_ge_01':len([v for v in pbs_vals if v>=0.1]),
    'n_negative':len([v for v in pbs_vals if v<0]),
    'pct_negative':len([v for v in pbs_vals if v<0])/n*100,
    'n_tier1':len(tier1),'n_tier2':len(tier2),'n_tier3':len(tier3),'n_candidates':len(cands)
}
print("\n=== PBS SUMMARY ===")
for k,v in stats.items():
    print(f"  {k}: {v:.6f}" if isinstance(v,float) else f"  {k}: {v}")

json.dump(stats, open("pbs_stats.json","w"), indent=2)
cands.sort(key=lambda r: -r['pbs_uzb'])
json.dump(cands, open("pbs_candidates.json","w"), indent=2)

with open("pbs_all.tsv","w") as f:
    f.write("SNP\tCHR\tBP\tPBS_UZB\tPBS_EUR\tPBS_EAS\tMAF_UZB\tMAF_EUR\tMAF_EAS\tMAF_SAS\tMAF_AFR\tDELTA_AF\tNEAR_PRIVATE\tTIER1\tTIER2\tTIER3\n")
    for r in sorted(results, key=lambda x: -x['pbs_uzb']):
        f.write(f"{r['snp']}\t{r['chr']}\t{r['bp']}\t{r['pbs_uzb']:.6f}\t{r['pbs_eur']:.6f}\t{r['pbs_eas']:.6f}\t{r['maf_uzb']:.4f}\t{r['maf_eur']:.4f}\t{r['maf_eas']:.4f}\t{r['maf_sas']:.4f}\t{r['maf_afr']:.4f}\t{r['delta_af']:.4f}\t{1 if r['near_private'] else 0}\t{1 if r['tier1'] else 0}\t{1 if r['tier2'] else 0}\t{1 if r['tier3'] else 0}\n")

print(f"\n=== TOP 20 PBS HITS ===")
print(f"{'#':>3} {'SNP':>20} {'Chr':>3} {'Pos':>12} {'PBS_UZB':>10} {'MAF_UZB':>8} {'MAF_EUR':>8} {'MAF_EAS':>8} {'MAF_SAS':>8} {'MAF_AFR':>8}")
for i,r in enumerate(cands[:20]):
    print(f"{i+1:3d} {r['snp']:>20} {r['chr']:>3} {r['bp']:>12d} {r['pbs_uzb']:10.4f} {r['maf_uzb']:8.4f} {r['maf_eur']:8.4f} {r['maf_eas']:8.4f} {r['maf_sas']:8.4f} {r['maf_afr']:8.4f}")

bins=np.arange(-1.5,5.05,0.05)
hist,edges=np.histogram(pbs_vals,bins=bins)
json.dump([{"x":round(float(edges[i]),3),"y":int(hist[i])} for i in range(len(hist)) if hist[i]>0],open("pbs_histogram.json","w"))

max_pbs=max(r['pbs_uzb'] for r in cands) if cands else 1
with open("pbs_clump_pval.txt","w") as f:
    f.write("SNP\tP\n")
    for r in cands:
        p=max(1.0-r['pbs_uzb']/(max_pbs+0.001),1e-300)
        f.write(f"{r['snp']}\t{p:.15e}\n")
with open("pbs_candidate_snps.txt","w") as f:
    for r in cands: f.write(f"{r['snp']}\n")

fst_summ={}
for pair in ["UZB_EUR","UZB_EAS","UZB_SAS","UZB_AFR","EUR_EAS"]:
    lf=f"fst_{pair}.log"
    if os.path.exists(lf):
        for line in open(lf):
            if "Weighted" in line: fst_summ[f"{pair}_weighted"]=float(line.split()[-1])
            elif "Mean" in line and "Fst" in line: fst_summ[f"{pair}_mean"]=float(line.split()[-1])
json.dump(fst_summ,open("fst_summary.json","w"),indent=2)
print(f"\nCandidates: {len(cands)}")
PYEOF
echo ""

# -----------------------------------------------------------------------
# 6. LD clumping
# -----------------------------------------------------------------------
echo "=== Step 6: LD Clumping ==="
N_CAND=$(wc -l < pbs_candidate_snps.txt)
echo "  Candidates: $N_CAND"

if [ "$N_CAND" -gt 0 ]; then
    plink --bfile "$MERGED" \
          --keep keep_UZB.txt \
          --clump pbs_clump_pval.txt \
          --clump-snp-field SNP --clump-field P \
          --clump-p1 1 --clump-p2 1 \
          --clump-r2 0.5 --clump-kb 1000 \
          --allow-no-sex --silent \
          --out pbs_clumped 2>&1 || true

    if [ -f pbs_clumped.clumped ]; then
        N_LOCI=$(awk 'NR>1 && NF>0 {c++} END {print c+0}' pbs_clumped.clumped)
        echo "  Independent loci: $N_LOCI"
    fi

    plink --bfile "$MERGED" \
          --keep keep_UZB.txt \
          --extract pbs_candidate_snps.txt \
          --r2 --ld-window-kb 5000 --ld-window 99999 --ld-window-r2 0.1 \
          --allow-no-sex --silent \
          --out pbs_ld 2>&1 || true

    if [ -f pbs_ld.ld ]; then
        echo "  LD pairs: $(tail -n +2 pbs_ld.ld | wc -l)"
    fi
fi
echo ""

# -----------------------------------------------------------------------
# 7. LD Decay
# -----------------------------------------------------------------------
echo "=== Step 7: LD Decay ==="
awk 'BEGIN{srand(42)} {if(rand()<3000/NR) a[NR]=$2} END{for(i in a) print a[i]}' ${UZB_QC}.bim > decay_snps.txt
echo "  Decay SNPs: $(wc -l < decay_snps.txt)"

plink --bfile "$UZB_QC" \
      --extract decay_snps.txt \
      --r2 --ld-window-kb 2000 --ld-window 99999 --ld-window-r2 0.0 \
      --allow-no-sex --silent \
      --out ld_decay 2>&1 || true

if [ -f ld_decay.ld ]; then
    awk 'NR>1 && $1==$4 {
        d=($5>$2)?($5-$2):($2-$5); b=int(d/50000)*50;
        s[b]+=$7; c[b]++
    } END {for(b in s) printf "%d\t%.6f\t%d\n",b,s[b]/c[b],c[b]}' ld_decay.ld | sort -n > ld_decay_binned.tsv
    echo "  Decay bins: $(wc -l < ld_decay_binned.tsv)"
else
    echo "  No LD decay output"
fi
echo ""

# -----------------------------------------------------------------------
# 8. Compile JSON
# -----------------------------------------------------------------------
echo "=== Step 8: Compile JSON ==="

python3 << 'PYEOF2'
import json, os
from collections import Counter

clumps=[]
if os.path.exists("pbs_clumped.clumped"):
    cm={c['snp']:c for c in json.load(open("pbs_candidates.json"))}
    with open("pbs_clumped.clumped") as f:
        f.readline()
        for line in f:
            p=line.split()
            if len(p)<3: continue
            snp=p[2]; proxies=[]
            if len(p)>11 and p[11]!="NONE":
                for t in p[11].split(","):
                    ps=t.split("(")[0].strip()
                    if ps and ps!="NONE": proxies.append(ps)
            c=cm.get(snp,{})
            clumps.append({"chr":p[0],"bp":int(p[3]) if len(p)>3 else 0,"snp":snp,"n_total":len(proxies),"pbs":c.get('pbs_uzb',0),"maf":c.get('maf_uzb',0),"proxies":proxies})
    clumps.sort(key=lambda x:-x['pbs'])

ld_pairs=[]
if os.path.exists("pbs_ld.ld"):
    with open("pbs_ld.ld") as f:
        f.readline()
        for line in f:
            p=line.split()
            if len(p)>=7:
                ld_pairs.append({"snp_a":p[2],"snp_b":p[5],"kb":round(abs(int(p[4])-int(p[1]))/1000,1),"r2":round(float(p[6]),3)})
    ld_pairs.sort(key=lambda x:-x['r2'])

decay=[]
if os.path.exists("ld_decay_binned.tsv"):
    with open("ld_decay_binned.tsv") as f:
        for line in f:
            p=line.strip().split('\t')
            if len(p)>=3:
                kb=int(p[0])
                if kb<=2000: decay.append({"kb":kb,"mean_r2":round(float(p[1]),4),"n":int(p[2])})
    decay.sort(key=lambda x:x['kb'])

chr_c=Counter(c['chr'] for c in clumps)
ld_data={"n_input":len(json.load(open("pbs_candidates.json"))),"n_loci":len(clumps),"clumps":clumps,"decay_curve":decay,"pbs_ld_pairs":ld_pairs[:300],"pbs_ld_pairs_total":len(ld_pairs),"clumps_with_proxies":len([c for c in clumps if c['n_total']>0]),"chr_loci_counts":{str(i):chr_c.get(str(i),0) for i in range(1,23)}}
json.dump(ld_data,open("ld_data.json","w"),indent=2)
print(f"  ld_data.json: {len(clumps)} loci, {len(ld_pairs)} pairs, {len(decay)} decay bins")
PYEOF2

echo ""
echo "===== COMPLETE ====="
echo "End: $(date)"
ls -lh *.json *.tsv 2>/dev/null
