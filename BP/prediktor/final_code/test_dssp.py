#!/usr/bin/python

import Bio
import sys
from Bio.PDB import *
import sys
import urllib2
import csv

encode_sec_struc = {
    	'C':'1','H':'2','S':'3'
    }
encode_asa = {
    	'burried':'10','partially_burried':'11','exposed':'12'
    }

#get pdb file for entered pdb id
def get_pdb(pdb_id):
    """download PDB file"""
    protein_file = pdb_id + ".pdb"
    with open(protein_file, "w") as pdb_output:
        url_pdb = "https://files.rcsb.org/download/%s.pdb" %pdb_id
        try:
            handle = urllib2.urlopen(url_pdb)
        except URLError as error:
            print(error.reason)
            sys.exit(1)
        pdb_output.write(handle.read())
        pdb_output.close()

def PDB_parse(name,chain_type):
    """get start position of residue from PDB file"""
    parser = PDBParser()
    structure = parser.get_structure(name,name+".pdb")
    model = structure[0]
    try:
        chain = model[chain_type]
    except KeyError as error:
        print(error.reason)
    residue_list = Selection.unfold_entities(chain,'A')
    residue_start = residue_list[0].get_full_id()[3][1]
    return residue_start

def compute_ASA(pdb_id, residue_start, index, wild_type):
    """compute accessible surface area and secondary structure"""
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, pdb_id + ".pdb")
    model = structure[0]
    dssp = DSSP(model, pdb_id + ".pdb", dssp='/home/juraj/dssp-2.0.4-linux-amd64')
    # index is position of mutation
    residue_start = int(residue_start)
    index = int(index)
    asa_keys = list(dssp.keys())
    print(asa_keys)

    if residue_start == 1:
        position = index - residue_start +1
    elif residue_start == 0:
        position = index-1
    elif residue_start < 0:
        position = index - abs(residue_start)
    else:
        position = index + residue_start -1
    for key in asa_keys:
        if(key[1][1] == position):
            if(dssp[key][1] == wild_type):
                asa_key = key
            else:
                position = index
                for key in asa_keys:
                    if(key[1][1] == position):
                        asa_key = key

    print(asa_key)
    print(dssp[asa_key])
    asa = str(dssp[asa_key][3])
    sec_structure = str(dssp[asa_key][2])
    sec_struc_id = sec_structure
    return sec_struc_id,asa

def sec_struc_code(id):
    """encode secondary structure"""
    if(id == 'H' or id == 'G' or id == 'I'):
    	sec_struc_res = 'H'
    elif(id == 'E' or id == 'B'):
    	sec_struc_res = 'S'
    elif(id == 'T' or id == 'S' or id == '-'):
    	sec_struc_res = 'C'
    return sec_struc_res

def divide_asa(asa):
    """encode asa values"""
    encoded_asa = float(asa)
    if(encoded_asa < 0.25 ):
        return "burried"
    elif (encoded_asa >= 0.25 and encoded_asa <= 0.5):
    	return "partially_burried"
    elif (encoded_asa > 0.5):
    	return "exposed"

handle = open('compiled_results_edit1.txt','r')
test_file = csv.reader(handle, delimiter=',')
asa_file = open('asa_results.txt','w')
sec_struc_file = open('structure_results.txt','w')

i = 0
for record in test_file:
    if i==0:
        i = 1
        continue
    pdb_id = record[0]
    chain = record[1]
    mutpos = record[2]
    wild_type = mutpos[0]
    mutated = mutpos[-1]
    position = mutpos[1:len(mutpos)-1]
    get_pdb(pdb_id)
    residue_start = PDB_parse(pdb_id,chain)
    structure, asa = compute_ASA(pdb_id,residue_start,position,wild_type)

    asa_code = divide_asa(asa)
    final_asa = encode_asa[asa_code]
    asa_file.write(final_asa+'\n')

    struc_code = sec_struc_code(structure)
    final_structure = encode_sec_struc[struc_code]
    sec_struc_file.write(final_structure+'\n')

asa_file.close()
sec_struc_file.close()

