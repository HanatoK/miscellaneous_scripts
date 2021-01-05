#!/usr/bin/env python3
import sys
import json
import glob
import os
from urllib.request import urlopen
from urllib.parse import quote
from urllib.error import HTTPError
from urllib.error import URLError
from http.client import HTTPException
from socket import timeout

max_retry = 100

def get_query_url(peptide_id):
    # DBAASP server
    server_base_url = r'https://dbaasp.org'
    query_url = server_base_url + f'/peptides/{peptide_id}'
    return query_url

def scrape_json(json_filename, output_prefix):
    print(f'Open {json_filename}')
    # Open the json file database
    with open(json_filename, 'r') as fdatabase:
        data = json.load(fdatabase)
        for item in data['data']:
            # Query every peptide id
            peptide_id = item['id']
            query_url = get_query_url(peptide_id)
            print(f'Request {query_url}')
            output_filename = f'{output_prefix}_{peptide_id}.dat'
            if os.path.exists(output_filename):
                print(f'Skip {output_filename}')
                continue
            retry = 0
            while True:
                try:
                    json_string = urlopen(query_url).read().decode('utf-8')
                    break
                except (HTTPError, URLError, HTTPException, timeout) as e:
                    print(f'Failed to read URL: {query_url}')
                    print(f'Error code: {e}')
                    if retry < max_retry:
                        retry = retry + 1
                        print(f'Retrying... ({retry})')
                        continue
                    else:
                        raise RuntimeError('Maximum retries reached.')
            with open(output_filename, 'w') as foutput:
                json_reader = json.loads(json_string)
                foutput.write(json.dumps(json_reader, indent=4))
                print(f'Saved to {output_filename}')

output_prefix = 'peptides_info/monomer'
# database files list
file_list = sorted(glob.glob('database/monomer_Gram-_[0-9][0-9][0-9][0-9][0-9].dat'))
for filename in file_list:
    scrape_json(filename, output_prefix)
