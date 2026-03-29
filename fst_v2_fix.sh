#!/bin/bash
set -euo pipefail
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

GLOBAL=~/v2/global_admixture/global_v2_admix
OUT=~/v2/fst_v2

log "===== FIX: Reassign all non-1000G samples as UZB ====="

# Known 1000G superpopulations in pop_mapping
# Any sample NOT labeled EUR/EAS/SAS/AFR/AMR → UZB
POP_MAP=~/v2/pca/pop_mapping.txt

awk -v pm="$POP_MAP" '
BEGIN {
    while ((getline line < pm) > 0) {
        split(line, a, "\t")
        pop[a[1]] = a[2]
    }
}
{
    fid = $1; iid = $2
    p = "UZB"  # default
    if (iid in pop) {
        v = pop[iid]
        if (v == "EUR" || v == "EAS" || v == "SAS" || v == "AFR" || v == "AMR") p = v
    } else if (fid in pop) {
        v = pop[fid]
        if (v == "EUR" || v == "EAS" || v == "SAS" || v == "AFR" || v == "AMR") p = v
    }
    print fid, iid, p
}
' ${GLOBAL}.fam > $OUT/all_pop_labels_fixed.txt

log "Fixed population counts:"
awk '{print $3}' $OUT/all_pop_labels_fixed.txt | sort | uniq -c | sort -rn

# Recompute only UZB pairs
POPS="EUR EAS SAS AFR"
echo -e "Pop1\tPop2\tMeanFst\tWeightedFst\tValidMarkers" > $OUT/fst_v2_results_fixed.tsv

for pop2 in $POPS; do
    log "  Computing Fst: UZB vs $pop2"

    awk -v p2="$pop2" '$3 == "UZB" || $3 == p2 {print $1, $2}' \
      $OUT/all_pop_labels_fixed.txt > $OUT/keep_UZB_${pop2}_fix.txt

    awk -v p2="$pop2" '$3 == "UZB" || $3 == p2 {print $1, $2, $3}' \
      $OUT/all_pop_labels_fixed.txt > $OUT/within_UZB_${pop2}_fix.txt

    n1=$(awk '$3=="UZB"' $OUT/within_UZB_${pop2}_fix.txt | wc -l)
    n2=$(awk -v p="$pop2" '$3==p' $OUT/within_UZB_${pop2}_fix.txt | wc -l)
    log "    UZB: $n1 samples, $pop2: $n2 samples"

    plink --bfile $GLOBAL \
      --keep $OUT/keep_UZB_${pop2}_fix.txt \
      --within $OUT/within_UZB_${pop2}_fix.txt \
      --fst \
      --out $OUT/fst_UZB_vs_${pop2}_fix \
      --silent

    mean_fst=$(grep "Mean Fst estimate:" $OUT/fst_UZB_vs_${pop2}_fix.log | awk '{print $NF}')
    wt_fst=$(grep "Weighted Fst estimate:" $OUT/fst_UZB_vs_${pop2}_fix.log | awk '{print $NF}')
    valid=$(awk 'NR>1 && $5 != "nan" && $5 != "NA"' $OUT/fst_UZB_vs_${pop2}_fix.fst 2>/dev/null | wc -l)

    log "    Mean=$mean_fst, Weighted=$wt_fst, Valid=$valid"
    echo -e "UZB\t${pop2}\t${mean_fst}\t${wt_fst}\t${valid}" >> $OUT/fst_v2_results_fixed.tsv
done

# Append the correct reference-population pairs from first run
grep -v "^UZB" $OUT/fst_v2_results.tsv | grep -v "^Pop1" >> $OUT/fst_v2_results_fixed.tsv

log "=== CORRECTED FST V2 RESULTS ==="
cat $OUT/fst_v2_results_fixed.tsv
log "=== END ==="

# Package
base64 $OUT/fst_v2_results_fixed.tsv
