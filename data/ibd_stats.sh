#!/bin/bash
echo MARKER_START

# Get IBD pairs with PI_HAT > 0.125 (3rd-degree relatives)
echo "=== IBD_RELATED ==="
awk 'NR>1 && $10>0.125{print $1,$2,$3,$4,$7,$8,$9,$10}' /staging/ALSU-analysis/winter2025/PLINK_301125_0312/ConvSK_mind20_ibd.genome | wc -l
echo "=== IBD_BY_DEGREE ==="
awk 'NR>1{pi=$10; if(pi>0.98) dup++; else if(pi>0.354) d1++; else if(pi>0.177) d2++; else if(pi>0.0884) d3++} END{printf "Duplicates/MZ(>0.98)=%d First_degree(>0.354)=%d Second_degree(>0.177)=%d Third_degree(>0.0884)=%d\n",dup,d1,d2,d3}' /staging/ALSU-analysis/winter2025/PLINK_301125_0312/ConvSK_mind20_ibd.genome

# Get the top related pairs
echo "=== IBD_TOP20 ==="
awk 'NR>1{print}' /staging/ALSU-analysis/winter2025/PLINK_301125_0312/ConvSK_mind20_ibd.genome | sort -k10 -rn | head -20

# Get PI_HAT distribution
echo "=== PIHAT_DIST ==="
awk 'NR>1{pi=$10; if(pi<0.05)a++; else if(pi<0.10)b++; else if(pi<0.20)c++; else if(pi<0.50)d++; else e++} END{printf "pihat<0.05=%d 0.05-0.10=%d 0.10-0.20=%d 0.20-0.50=%d pihat>0.50=%d\n",a,b,c,d,e}' /staging/ALSU-analysis/winter2025/PLINK_301125_0312/ConvSK_mind20_ibd.genome

# Get ADMIXTURE Q-values for Uzbek-only K=2
echo "=== ADMIX_HEAD ==="
head -5 /staging/ALSU-analysis/admixture_analysis/UZB_for_admixture.2.Q
echo "=== ADMIX_N ==="
wc -l /staging/ALSU-analysis/admixture_analysis/UZB_for_admixture.2.Q
echo "=== FAM_HEAD ==="
head -5 /staging/ALSU-analysis/admixture_analysis/UZB_for_admixture.fam

# Merge FAM IDs with Q-values and compute stats
echo "=== Q2_STATS ==="
paste /staging/ALSU-analysis/admixture_analysis/UZB_for_admixture.fam /staging/ALSU-analysis/admixture_analysis/UZB_for_admixture.2.Q | awk '{q1=$7; q2=$8; sum1+=q1; sum2+=q2; count++} END{printf "N=%d mean_Q1=%.4f mean_Q2=%.4f\n",count,sum1/count,sum2/count}'

# ROH .hom.indiv as base64 (for embedding in HTML)
echo "=== INDIV_B64_START ==="
cat ~/ROH_analysis/UZB_ROH.hom.indiv | base64
echo "=== INDIV_B64_END ==="

echo MARKER_END
