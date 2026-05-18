#!/usr/bin/env python3
ids_48 = [
'03-328','03-72','03-73','03-85','07-42','07-79','07-96',
'08-290','08-291','08-292','08-296','08-297','08-298','08-305',
'08-306','08-309','08-310','08-311','08-314','08-317','08-321',
'08-322','08-323','08-326','08-331','08-332','08-336','08-341',
'08-352','08-356','08-360','08-374','08-386','08-390','08-391',
'08-392','08-399','08-407','08-410','08-411','08-416','08-419',
'08-421','08-430','08-433','08-487','08-490','08-505'
]

convsk = set()
with open('/staging/ALSU-analysis/spring2026/ConvSK_raw.fam') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            convsk.add(parts[1])

in_convsk = [x for x in ids_48 if x in convsk]
not_in_convsk = [x for x in ids_48 if x not in convsk]
print(f'Total 48 rescan: {len(ids_48)}')
print(f'IN ConvSK_raw.fam ({len(in_convsk)}): {in_convsk}')
print(f'NOT in ConvSK_raw.fam ({len(not_in_convsk)}): {not_in_convsk}')

# Also check gwas2026_raw.fam
gwas = set()
try:
    with open('/staging/ALSU-analysis/admixture_analysis/temp/dragen_array_test_20260513/gwas2026/plink/gwas2026_raw.fam') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                gwas.add(parts[1])
    in_gwas = [x for x in ids_48 if x in gwas]
    print(f'IN gwas2026_raw.fam ({len(in_gwas)}): {in_gwas}')
except Exception as e:
    print(f'gwas check error: {e}')
