echo "=== ADMIXTURE STATUS ==="
tail -1 ~/v2/reanalysis.log 2>/dev/null
echo "---"
ls -la ~/v2/admixture/*.Q 2>/dev/null
echo "=== CV ERRORS ==="
grep "CV error" ~/v2/admixture/admixture_K*.log 2>/dev/null
echo "=== LOG LIKELIHOODS ==="
grep "Loglikelihood" ~/v2/admixture/admixture_K*.log 2>/dev/null
echo "=== RUNNING ==="
ps aux | grep -i admixture | grep -v grep
echo "=== FAM HEAD ==="
wc -l ~/v2/admixture/UZB_v2_admix.fam 2>/dev/null
echo "=== Q FILES BASE64 ==="
for k in 2 3 4 5; do
  echo "--- K=$k ---"
  base64 ~/v2/admixture/UZB_v2_admix.${k}.Q 2>/dev/null
done
echo "=== DONE ==="