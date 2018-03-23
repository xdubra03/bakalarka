#!/usr/bin/python

import sys
import subprocess
import pandas as pd
import csv

dataset = open('dataset_S350_cleaned2.csv', 'r')
test_file = csv.reader(dataset, delimiter=',')

pdb_id = ''
wild_type = ''
mutant = ''
position = 0
chain = ''
predicted_class = []


for record in test_file:
    pdb_id = record[0]
    wild_type = record[1]
    position = record[2]
    mutant = record[3]
    chain = record[8]
    print(pdb_id,wild_type,position,mutant,chain)