import json, os

# Load old (has PBS data) and new (has K=2-8)
old = json.load(open('data/admix_data.json'))
new = json.load(open('data/admix_data_full.json'))

# Merge PBS into new
new['pbs_top'] = old['pbs_top']
new['pbs_header'] = old['pbs_header']

# Write back
json.dump(new, open('data/admix_data.json', 'w'))
sz = os.path.getsize('data/admix_data.json')
print(f'Updated admix_data.json: {sz:,} bytes')
print(f'Keys: {list(new.keys())}')
