#!/bin/bash
echo MARKER_START
echo "=== INDIV HEADER ==="
head -1 ~/ROH_analysis/UZB_ROH.hom.indiv
echo "=== STATS ==="
awk 'NR>1{sum_kb+=$5; sum_n+=$4; count++} END{printf "N=%d mean_n=%.1f mean_kb=%.1f mean_Mb=%.2f\n",count,sum_n/count,sum_kb/count,sum_kb/count/1000}' ~/ROH_analysis/UZB_ROH.hom.indiv
echo "=== PCTL ==="
awk 'NR>1{kb[NR-1]=$5; n[NR-1]=$4; count++} END{asort(kb); asort(n); printf "n_min=%d n_p25=%d n_med=%d n_p75=%d n_max=%d\n",n[1],n[int(count*0.25)],n[int(count*0.5)],n[int(count*0.75)],n[count]; printf "kb_min=%.0f kb_p25=%.0f kb_med=%.0f kb_p75=%.0f kb_max=%.0f\n",kb[1],kb[int(count*0.25)],kb[int(count*0.5)],kb[int(count*0.75)],kb[count]; printf "Froh_med=%.6f Froh_p75=%.6f Froh_max=%.6f\n",kb[int(count*0.5)]/2881000,kb[int(count*0.75)]/2881000,kb[count]/2881000}' ~/ROH_analysis/UZB_ROH.hom.indiv
echo "=== SIZE ==="
awk 'NR>1{kb=$8; if(kb<1000)s++; else if(kb<5000)m++; else l++} END{printf "lt1Mb=%d 1to5Mb=%d gt5Mb=%d tot=%d\n",s,m,l,s+m+l}' ~/ROH_analysis/UZB_ROH.hom
echo "=== TOP15 ==="
sort -k5 -rn ~/ROH_analysis/UZB_ROH.hom.indiv | head -16
echo "=== BOT5 ==="
sort -k5 -n ~/ROH_analysis/UZB_ROH.hom.indiv | head -6
echo "=== 5MB_N ==="
awk 'NR>1 && $8>=5000{seen[$2]=1} END{print length(seen)}' ~/ROH_analysis/UZB_ROH.hom
echo "=== 5MB_TOP ==="
awk 'NR>1 && $8>=5000{c[$2]++; s[$2]+=$8} END{for(k in c) printf "%s %d %.0f\n",k,c[k],s[k]}' ~/ROH_analysis/UZB_ROH.hom | sort -k3 -rn | head -10
echo "=== CHR_DIST ==="
awk 'NR>1{chr=$4; kb[chr]+=$8; n[chr]++} END{for(c in n) printf "chr%s n=%d total_kb=%.0f\n",c,n[c],kb[c]}' ~/ROH_analysis/UZB_ROH.hom | sort -t= -k2 -rn | head -22
echo MARKER_END
