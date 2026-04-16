@echo off
echo Checking DRAGEN ADMIXTURE status via tunnel (port 2222)...
echo.
echo === Running processes ===
plink.exe -batch -P 2222 copilot@127.0.0.1 "ps aux | grep -i admix | grep -v grep || echo No admixture process running"
echo.
echo === Log tail (last 30 lines) ===
plink.exe -batch -P 2222 copilot@127.0.0.1 "tail -30 /staging/ALSU-analysis/admixture_analysis/admixture_all.log 2>/dev/null || echo Log file not found"
echo.
echo === Admixture output files ===
plink.exe -batch -P 2222 copilot@127.0.0.1 "ls -lh /staging/ALSU-analysis/admixture_analysis/*.Q /staging/ALSU-analysis/admixture_analysis/*.P 2>/dev/null || echo No .Q/.P files found"
echo.
pause
