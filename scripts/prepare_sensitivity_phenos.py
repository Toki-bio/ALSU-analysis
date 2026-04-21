"""
Build alternate phenotype files for GWAS sensitivity analyses.

Reads data/gwas_sample_map.txt (already merged genotype<->phenotype
classification) and emits four alternative pheno files:

  S1_strict3        CASE = >=3 strict losses (ESHRE definition).
                    CONTROL same as primary (>=1 live, 0 losses).
  S2_review_as_case CASE = primary CASE + REVIEW (loss at GA>=20).
                    CONTROL same.
  S3_gray_as_case   CASE = primary CASE + GRAY (any loss = >=1 sporadic).
                    CONTROL same.
  S4_probable       CASE = primary CASE + CASE_PROBABLE (>=2 total losses
                    when GA missing). CONTROL same + CONTROL_PROBABLE.

Output format matches data/gwas_pheno.txt:
  #FID  IID  RPL    (1=case, 0=control, NA=excluded)
"""
import os, sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAP = os.path.join(ROOT, 'data', 'gwas_sample_map.txt')
OUT_DIR = os.path.join(ROOT, 'data', 'sensitivity')
os.makedirs(OUT_DIR, exist_ok=True)


def parse_map():
    rows = []
    with open(MAP) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#') or not line:
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            fid, iid, _, cls, details = parts[:5]
            # details: "loss=N live=N preg=N"
            d = {}
            for tok in details.split():
                if '=' in tok:
                    k, v = tok.split('=', 1)
                    try:
                        d[k] = int(v)
                    except ValueError:
                        pass
            rows.append({'fid': fid, 'iid': iid, 'cls': cls,
                         'loss': d.get('loss', 0), 'live': d.get('live', 0),
                         'preg': d.get('preg', 0)})
    return rows


def write_pheno(rows, name, classifier):
    """classifier(row) -> 1, 0, or None (None=NA/excluded)"""
    counts = Counter()
    out = os.path.join(OUT_DIR, f'gwas_pheno_{name}.txt')
    with open(out, 'w') as f:
        f.write('#FID\tIID\tRPL\n')
        for r in rows:
            v = classifier(r)
            if v is None:
                rpl = 'NA'
                counts['NA'] += 1
            else:
                rpl = str(v)
                counts['CASE' if v == 1 else 'CONTROL'] += 1
            f.write(f"{r['fid']}\t{r['iid']}\t{rpl}\n")
    return out, counts


def main():
    rows = parse_map()
    print(f'Parsed {len(rows)} sample rows')
    print('Class distribution:', Counter(r['cls'] for r in rows).most_common())

    # ── S1: strict ESHRE >=3 strict losses (CASE only kept if loss>=3) ──
    def s1(r):
        if r['cls'] == 'CASE' and r['loss'] >= 3:
            return 1
        if r['cls'] == 'CONTROL':
            return 0
        return None

    # ── S2: include REVIEW (GA>=20 losses likely miscoded stillbirths) ──
    def s2(r):
        if r['cls'] == 'CASE':
            return 1
        if r['cls'] == 'REVIEW':
            return 1
        if r['cls'] == 'CONTROL':
            return 0
        return None

    # ── S3: include GRAY (any loss) — broadest case definition ──
    def s3(r):
        if r['cls'] == 'CASE' or r['cls'] == 'GRAY':
            return 1
        if r['cls'] == 'CONTROL':
            return 0
        return None

    # ── S4: include PROBABLE (cases with missing GA + control with current preg) ──
    def s4(r):
        if r['cls'] in ('CASE', 'CASE_PROBABLE'):
            return 1
        if r['cls'] in ('CONTROL', 'CONTROL_PROBABLE'):
            return 0
        return None

    for name, fn in [('S1_strict3', s1),
                     ('S2_review_as_case', s2),
                     ('S3_gray_as_case', s3),
                     ('S4_probable', s4)]:
        path, c = write_pheno(rows, name, fn)
        print(f'  {name:22s}  CASE={c["CASE"]:4d}  CONTROL={c["CONTROL"]:4d}  '
              f'NA={c["NA"]:4d}  ratio={c["CASE"]/max(1,c["CONTROL"]):.3f}  -> {path}')


if __name__ == '__main__':
    sys.exit(main())
