import json

data = json.load(open('data/admix_data.json'))

# Print as compact JS
print('    const ADMIX_DATA = ' + json.dumps(data, separators=(',', ':')) + ';')
