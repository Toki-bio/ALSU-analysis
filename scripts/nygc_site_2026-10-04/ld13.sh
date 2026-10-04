#!/bin/bash
# LD decay, Uzbek (imputed HQ, 1,256) vs NYGC EUR/EAS/SAS/AFR on the same random SNP set.
# usage: ld13.sh <outdir-name> <chr list, e.g. "22" or "$(seq 1 22)"> <n_snps>
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
export TMPDIR=/staging/tmp/scratch
OUT=/staging/tmp/scratch/$1
CHRS="$2"
NSNP=$3
mkdir -p "$OUT"
cd "$OUT" || exit 1
F=/staging/ALSU-analysis/spring2026/full_expanded_cohort
U=$F/imputation_results/hq_filtered/UZB_imputed_HQ_qc
N=$F/nygc30x
K=$F/pbs_refqc_2026-10/nygc30x_1256
echo "start $(date) chrs=[$CHRS] nsnp=$NSNP"
# 1. Uzbek SNV list with MAF >= 0.05
plink --bfile $U --freq --out uzb_freq --threads 8 > /dev/null
awk 'NR==FNR{bp[$2]=$4; next} FNR>1 && length($3)==1 && length($4)==1 && $5>=0.05 && ($2 in bp) {print $1"\t"$2"\t"bp[$2]"\t"$3$4}' $U.bim uzb_freq.frq > cand.tsv
for C in $CHRS; do awk -v c=$C '$1==c' cand.tsv; done > cand_sel.tsv
wc -l cand_sel.tsv
# 2. reproducible random subset (seed 42), proportional to the candidate pool of the selected chromosomes
awk 'BEGIN{srand(42)} {print rand()"\t"$0}' cand_sel.tsv | sort -k1,1g | awk -v n=$NSNP 'NR<=n' | cut -f2- > snps.tsv
cut -f2 snps.tsv > snps.txt
wc -l snps.txt
# 3. Uzbek r2 (all chromosomes at once; pairs are within chromosomes only)
plink --bfile $U --extract snps.txt --maf 0.05 --r2 --ld-window-kb 2000 --ld-window 99999 --ld-window-r2 0 --out uzb --threads 8 > uzb.plink.log 2>&1
grep -i "variants remaining\|pass filters" uzb.plink.log || true
# 4. NYGC populations per chromosome
for C in $CHRS; do
  awk -F'\t' -v c=$C '$1==c {print "chr"$1"\t"$3-1"\t"$3"\t"$4}' snps.tsv | sort -k2,2n -u > reg$C.bed
  [ -s reg$C.bed ] || continue
  bcftools view -R reg$C.bed -m2 -M2 -v snps -Oz -o sub$C.vcf.gz $N/1kGP_high_coverage_Illumina.chr$C.filtered.SNV_INDEL_SV_phased_panel.vcf.gz 2> /dev/null
  for P in EUR EAS SAS AFR; do
    plink --vcf sub$C.vcf.gz --double-id --keep $K/keep_$P.txt --maf 0.05 --r2 --ld-window-kb 2000 --ld-window 99999 --ld-window-r2 0 --out r2_${P}_$C --threads 4 > r2_${P}_$C.plink.log 2>&1 || true
  done
  echo "chr$C done $(date)"
done
# 5. bin: distance bins in kb; mean r2 and pair counts
python3 - <<'PY'
import glob, collections, json
bins = [(0,10),(10,25),(25,50),(50,100),(100,200),(200,500),(500,1000),(1000,2000)]
def b(d):
    for i,(lo,hi) in enumerate(bins):
        if lo*1000 <= d < hi*1000: return i
    return None
res = {}
for P, files in {"UZB": ["uzb.ld"], "EUR": glob.glob("r2_EUR_*.ld"), "EAS": glob.glob("r2_EAS_*.ld"), "SAS": glob.glob("r2_SAS_*.ld"), "AFR": glob.glob("r2_AFR_*.ld")}.items():
    s = collections.defaultdict(float); n = collections.defaultdict(int)
    for fn in files:
        with open(fn) as fh:
            next(fh)
            for l in fh:
                f = l.split()
                if f[0] != f[3]: continue
                i = b(abs(int(f[4]) - int(f[1])))
                if i is None: continue
                r = float(f[6]) if f[6] != "nan" else None
                if r is None: continue
                s[i] += r; n[i] += 1
    res[P] = [[bins[i][0], bins[i][1], (s[i]/n[i] if n[i] else None), n[i]] for i in range(len(bins))]
json.dump(res, open("ld_bins.json", "w"))
for P, r in res.items(): print(P, [(x[0], round(x[2], 4) if x[2] is not None else None, x[3]) for x in r])
PY
echo "finished $(date)"
