import subprocess, os, math, json, itertools
os.chdir("/staging/tmp/scratch/fst22_check")
F = "/staging/ALSU-analysis/spring2026/full_expanded_cohort"
U = F + "/imputation_results/hq_filtered/UZB_imputed_HQ_qc"
K = F + "/pbs_refqc_2026-10/nygc30x_1256"
PL = "plink"


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# frequencies per population at chr22 SNVs
for P in ["EUR", "EAS", "SAS", "AFR"]:
    run(f"{PL} --bfile c22 --keep {K}/keep_{P}.txt --freq --out f_{P} --threads 8")
run(f"{PL} --bfile {U} --chr 22 --maf 0.01 --freq --out f_UZB --threads 8")


def readfrq(fn, bim=None):
    pos = {}
    if bim:
        for l in open(bim):
            p = l.split(); pos[p[1]] = p[3]
    d = {}
    for l in list(open(fn))[1:]:
        p = l.split()
        if p[4] == "NA" or p[5] == "0":
            continue
        key = pos[p[1]] if bim else p[1].split(":")[1] if ":" in p[1] else None
        d[key] = (p[2], p[3], float(p[4]), int(p[5]))   # A1, A2, MAF(of A1), NCHROBS
    return d


# NYGC bim ids
nb = {}
for l in open("c22.bim"):
    p = l.split(); nb[p[1]] = p[3]


def readfrq2(fn, idpos):
    d = {}
    for l in list(open(fn))[1:]:
        p = l.split()
        if p[4] == "NA" or p[5] == "0":
            continue
        d[idpos[p[1]]] = (p[2], p[3], float(p[4]), int(p[5]))
    return d


ub = {}
for l in open(U + ".bim"):
    p = l.split()
    if p[0] == "22":
        ub[p[1]] = p[3]
fr = {P: readfrq2(f"f_{P}.frq", nb) for P in ["EUR", "EAS", "SAS", "AFR"]}
fr["UZB"] = readfrq2("f_UZB.frq", ub)
for P, d in fr.items():
    print(P, len(d))


def aligned(a, b):
    """return (p_a, p_b) as frequency of a's A1 allele in both, or None"""
    A1, A2, pa, na = a
    B1, B2, pb, nb_ = b
    if len(A1) != 1 or len(A2) != 1 or len(B1) != 1 or len(B2) != 1:
        return None
    if A1 == B1 and A2 == B2:
        return pa, pb
    if A1 == B2 and A2 == B1:
        return pa, 1 - pb
    return None


def hudson(P1, P2):
    num = den = 0.0; n = 0
    for pos, a in fr[P1].items():
        b = fr[P2].get(pos)
        if b is None:
            continue
        al = aligned(a, b)
        if al is None:
            continue
        p1, p2 = al; n1, n2 = a[3], b[3]
        if n1 < 2 or n2 < 2:
            continue
        nu = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        de = p1 * (1 - p2) + p2 * (1 - p1)
        if de <= 0:
            continue
        num += nu; den += de; n += 1
    return num / den, n


out = {}
for a, b in itertools.combinations(["UZB", "EUR", "SAS", "EAS", "AFR"], 2):
    v, n = hudson(a, b)
    out[f"{a}_{b}"] = dict(hudson=v, n=n)
    print(a, b, round(v, 4), n)
json.dump(out, open("hudson_chr22.json", "w"), indent=1)
