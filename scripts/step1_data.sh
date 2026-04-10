#!/bin/bash
cd /staging/ALSU-analysis/spring2026

echo "=== HISTOGRAM ==="
awk 'NR>1{b=int(100*($6+0)); if(b>40)b=40; bins[b]++} END{for(i=0;i<=40;i++) printf "%d\t%.2f-%.2f\t%d\n",i,i/100.0,(i+1)/100.0,bins[i]+0}' ConvSK_raw_miss.imiss

echo ""
echo "=== REMOVED_99 ==="
while IFS=$'\t' read -r fid iid; do
  awk -v f="$fid" -v i="$iid" '$1==f && $2==i {print $1,$2,$6}' ConvSK_raw_miss.imiss
done < remove_miss20.txt

echo ""
echo "=== BASIC_STATS ==="
awk 'NR>1{s+=$6; ss+=$6*$6; n++; if($6+0>mx)mx=$6} END{printf "n=%d mean=%.6f sd=%.6f max=%.6f\n",n,s/n,sqrt(ss/n-(s/n)^2),mx}' ConvSK_raw_miss.imiss

echo ""
echo "=== RETAINED_STATS ==="
awk 'NR>1 && $6+0<=0.20{s+=$6; ss+=$6*$6; n++; if($6+0>mx)mx=$6; if(n==1||$6+0<mn)mn=$6} END{printf "n=%d mean=%.6f sd=%.6f min=%.6f max=%.6f\n",n,s/n,sqrt(ss/n-(s/n)^2),mn,mx}' ConvSK_raw_miss.imiss

echo ""
echo "=== MIND20_EXISTS ==="
ls -la ConvSK_mind20.bed ConvSK_mind20.bim ConvSK_mind20.fam 2>/dev/null || echo "NOT_YET_CREATED"

echo ""
echo "=== SEVEN_EXTRA ==="
diff <(awk '{print $2}' remove_miss20.txt | sort) <(awk '{print $2}' /staging/ALSU-analysis/winter2025/PLINK_301125_0312/remove_miss20.txt | sort) | grep '^<' | sed 's/^< //'

echo ""
echo "=== SEVEN_FMISS ==="
for s in 08-365 08-701 08-25 08-495 08-77 08-825 12-11; do
  awk -v id="$s" '$2==id{print $1,$2,$6}' ConvSK_raw_miss.imiss
done

echo ""
echo "=== REMOVE_FORMAT ==="
head -3 remove_miss20.txt | cat -A
echo ""
echo "=== IMISS_FORMAT ==="
head -3 ConvSK_raw_miss.imiss | cat -A
