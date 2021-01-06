#!/usr/bin/env python3
import glob
from natsort import natsorted
import sqlite3
import json
from sqlite3 import Error

def property_list_to_dict(list_in):
    result = {}
    for item in list_in:
        result[item['name']] = float(item['value'])
    return result

def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return conn

def create_table(conn):
    sql = '''CREATE TABLE IF NOT EXISTS monomers (
                 peptideId INTEGER PRIMARY KEY,
                 name TEXT,
                 sequence TEXT NOT NULL,
                 synthesisType TEXT NOT NULL,
                 targetGroups TEXT,
                 PDB BLOB,
                 netCharge REAL,
                 targetObjects TEXT
    );'''
    try:
        cur = conn.cursor()
        cur.execute(sql)
    except Error as e:
        print(e)

def check_record(conn, peptideId):
    # check if the current id exists
    id_query = 'SELECT peptideId FROM monomers WHERE peptideId=?'
    cur = conn.cursor()
    cur.execute(id_query, (peptideId,))
    data = cur.fetchone()
    if data is None:
        return False
    else:
        return True

def create_monomer(conn, json_obj):
    cur = conn.cursor()
    peptideId = json_obj['id']
    if check_record(conn, peptideId):
        print(f'Monomer id {peptideId} exists. Skip.')
        return cur.lastrowid
    name = json_obj['name']
    sequence = json_obj['sequence']
    if sequence is None:
        print(f'No sequence for {peptideId}. Skip.')
        return cur.lastrowid
    synthesisType = json_obj['synthesisType']['name']
    if synthesisType is None:
        print(f'No synthesisType for {peptideId}. Skip.')
        return cur.lastrowid
    targetGroups = None
    if json_obj['targetGroups']:
        targetGroups = ''
        for item in json_obj['targetGroups']:
            targetGroups += '"' + item['name'] + '",'
        targetGroups = targetGroups.strip(',')
    PDB = None
    netCharge = None
    if json_obj['physicoChemicalProperties']:
        property_dict = property_list_to_dict(json_obj['physicoChemicalProperties'])
        if 'Net Charge' in property_dict.keys():
            netCharge = property_dict['Net Charge']
    targetObjects = None
    if json_obj['targetObjects']:
        targetObjects = ''
        for item in json_obj['targetObjects']:
            targetObjects += '"' + item['name'] + '",'
        targetObjects = targetObjects.strip(',')
    sql = '''INSERT INTO monomers(peptideId, name, sequence, synthesisType, targetGroups, PDB, netCharge, targetObjects)
             VALUES(?,?,?,?,?,?,?,?) '''
    monomer = (peptideId, name, sequence, synthesisType, targetGroups, PDB, netCharge, targetObjects)
    print(f'Insert {monomer}')
    cur.execute(sql, monomer)
    conn.commit()
    return cur.lastrowid

def main():
    database = 'peptides.db'
    conn = create_connection(database)
    create_table(conn)
    with conn:
        file_list = natsorted(glob.glob('peptides_info/monomer_*.dat'))
        for filename in file_list:
            with open(filename, 'r') as finput:
                json_obj = json.load(finput)
                create_monomer(conn, json_obj)

if __name__ == '__main__':
    main()
