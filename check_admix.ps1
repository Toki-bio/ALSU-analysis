$out = & cmd /c 'plink.exe -batch -P 2222 copilot@127.0.0.1 "ps aux | grep admix | grep -v grep; echo ===LOG===; tail -30 /staging/ALSU-analysis/admixture_analysis/admixture_all.log 2>/dev/null; echo ===END==="' 2>&1
$out | Out-File -FilePath "c:\work\Surface\work\ALSU-analysis\admix_check_result.txt" -Encoding utf8
Write-Host "RESULT:"
Write-Host ($out -join "`n")
