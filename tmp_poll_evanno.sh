cd /staging/ALSU-analysis/admixture_analysis/evanno_runs
for K in 2 3 4 5 6 7 8; do
  dirs=$(find . -maxdepth 1 -type d -name "K${K}_rep*" | wc -l)
  logs=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/*.log" | wc -l)
  q=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/UZB_for_admixture.${K}.Q" | wc -l)
  p=$(find . -maxdepth 2 -type f -path "./K${K}_rep*/UZB_for_admixture.${K}.P" | wc -l)
  echo "K${K} dirs=${dirs} logs=${logs} Q=${q} P=${p}"
done
