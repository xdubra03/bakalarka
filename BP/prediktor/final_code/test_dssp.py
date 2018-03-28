#!/usr/bin/python

import Bio
import sys
from Bio.PDB import *


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

sec_struc, asa = compute_ASA(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

print('SeC struc' + str(sec_struc))
print('asa: ' + str(asa))