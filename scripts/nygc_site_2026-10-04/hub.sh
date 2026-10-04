set -u
export PATH=/staging/conda/envs/bioinfo/bin:$PATH
B=/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/cross_gwas2026_alsu_winter
mkdir -p /staging/tmp/scratch/hub_check; cd /staging/tmp/scratch/hub_check || exit 1
plink --bfile $B/merged --extract $B/cross_prune.prune.in --missing --het --out q --threads 8 > /dev/null 2>&1
plink --bfile $B/merged --extract $B/cross_prune.prune.in --freq --out f --threads 8 > /dev/null 2>&1
ls
head -2 q.imiss q.het
