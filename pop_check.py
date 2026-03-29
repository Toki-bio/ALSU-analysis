import collections
lines = open('pop_mapping.txt').readlines()
print(f"Total lines: {len(lines)}")
print(f"Last 5:")
for l in lines[-5:]:
    print(f"  {l.strip()!r}")
pops = collections.Counter()
for l in lines:
    parts = l.strip().split('\t')
    if len(parts) >= 2:
        pops[parts[1]] += 1
for p, c in sorted(pops.items()):
    print(f"  {p}: {c}")
