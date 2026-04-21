$ErrorActionPreference = 'Stop'
$keyPath = "$env:USERPROFILE\.ssh\id_ed25519.ppk"
$remoteDir = "/tmp/chr18_locus"
$localDir = "data\chr18_locus"
New-Item -ItemType Directory -Force -Path $localDir | Out-Null

$files = @("snps_top.tsv", "snps_saige.tsv", "ld_r2_ours.phased.vcor2.vars", "ld_r2_ours.phased.vcor2")

foreach ($f in $files) {
    $raw = cmd /c plink.exe -batch -P 2222 -i $keyPath copilot@127.0.0.1 "echo __B64_START__; base64 -w0 $remoteDir/$f; echo; echo __B64_END__"
    $text = ($raw -join "`n")
    $startIdx = $text.IndexOf("__B64_START__") + "__B64_START__".Length
    $endIdx = $text.IndexOf("__B64_END__")
    if ($startIdx -lt 0 -or $endIdx -lt 0) { Write-Host "$f : marker not found"; continue }
    $b64 = $text.Substring($startIdx, $endIdx - $startIdx) -replace '\s', ''
    [IO.File]::WriteAllBytes("$localDir\$f", [Convert]::FromBase64String($b64))
    $len = (Get-Item "$localDir\$f").Length
    Write-Host "$f : $len bytes"
}
