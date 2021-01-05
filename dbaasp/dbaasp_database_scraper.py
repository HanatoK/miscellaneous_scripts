#!/usr/bin/env python3
import sys
import json
from urllib.request import urlopen
from urllib.parse import quote

def get_query_url(offset, limit, complexity, targetGroup):
    # DBAASP server
    server_base_url = r'https://dbaasp.org'
    query_url = (server_base_url
              + '/peptides?'
              + f'complexity.value={complexity}&'
              + f'limit={int(limit)}&'
              + f'offset={int(offset)}&'
              + f'targetGroup.value={quote(targetGroup)}')
    return query_url


complexity = 'multimer'
target_group = 'Gram+'
start = 0
batch_count = 100
query_url = get_query_url(start, batch_count, complexity, target_group)
print(f'Request URL: {query_url}')

# Get the first batch of sequences from the database
# Also read the total number of sequences
total_count = 0
output_filename = f'database/{complexity}_{target_group}_{str(0).zfill(5)}.dat'
with urlopen(query_url) as finput, open(output_filename, 'w') as foutput:
    json_string = finput.read().decode('utf-8')
    json_reader = json.loads(json_string)
    foutput.write(json.dumps(json_reader, indent=4))
    total_count = json_reader["totalCount"]
    print(f'Total number of sequences: {total_count}')
    print(f'Saved to {output_filename}')

# Get all the sequences
num_batches = int((total_count - total_count % batch_count) / batch_count)
for i in range(1, num_batches + 1):
    start = i * batch_count
    query_url = get_query_url(start, batch_count, complexity, target_group)
    print(f'Request URL: {query_url}')
    output_filename = f'database/{complexity}_{target_group}_{str(i).zfill(5)}.dat'
    with urlopen(query_url) as finput, open(output_filename, 'w') as foutput:
        json_string = finput.read().decode('utf-8')
        json_reader = json.loads(json_string)
        foutput.write(json.dumps(json_reader, indent=4))
        print(f'Saved to {output_filename}')
