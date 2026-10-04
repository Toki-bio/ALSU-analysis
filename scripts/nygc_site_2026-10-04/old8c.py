import subprocess, collections, sys, os
os.environ["PATH"] = "/staging/conda/envs/bioinfo/bin:" + os.environ["PATH"]
F = "/staging/ALSU-analysis/spring2026/full_expanded_cohort"
N = F + "/nygc30x/1kGP_high_coverage_Illumina.chr%s.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
K = F + "/pbs_refqc_2026-10/nygc30x_1256/keep_%s.txt"
sites = ["11:207698", "5:53879140", "3:133749168", "10:8512594", "11:20665570", "12:22967890", "12:125520190", "12:5664803"]
pops = {p: set(l.split()[1] for l in open(K % p)) for p in ["EUR", "EAS", "SAS", "AFR"]}
out = []
for s in sites:
    c, p = s.split(":")
    samp = subprocess.check_output(["bcftools", "query", "-l", N % c]).decode().split()
    rows = subprocess.check_output(["bcftools", "query", "-r", f"chr{c}:{p}-{p}", "-f", "%REF\t%ALT[\t%GT]\n", N % c]).decode().strip().split("\n")
    for r in rows:
        f = r.split("\t"); ref, alt, gts = f[0], f[1], f[2:]
        if len(ref) != 1 or len(alt) != 1: out.append((s, ref, alt, "non-SNV record")); continue
        res = {}
        for pop, ids in pops.items():
            n = a = 0
            for sm, g in zip(samp, gts):
                if sm in ids:
                    for x in g.replace("|", "/").split("/"):
                        if x in "01": n += 1; a += (x == "1")
            res[pop] = round(a / n, 4) if n else None
        out.append((s, ref, alt, res))
for o in out: print(o)
# Uzbek AF on the array data where present
bim = F + "/alsu_final.bim"
for s in sites:
    c, p = s.split(":")
    subprocess.run(["plink", "--bfile", F + "/alsu_final", "--snp", s, "--freq", "--out", "/staging/tmp/scratch/old8_" + s.replace(":", "_")], capture_output=True)
    fn = "/staging/tmp/scratch/old8_" + s.replace(":", "_") + ".frq"
    print(s, open(fn).read().split("\n")[1] if os.path.exists(fn) else "not in array set")
