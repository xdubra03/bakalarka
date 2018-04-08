#!/usr/bin/python

import sys
import subprocess
import pandas as pd
import csv

dataset = open('dataset_S350.csv', 'r')
test_file = csv.reader(dataset, delimiter=',')
i = 0
classes = []
for line in test_file:
    if(i == 0):
        i = 1
        continue
    if(float(line[11]) >= 0):
        classes.append(1)
    else:
        classes.append(-1)

print(len(classes))
f = pd.read_csv('dataset_S350.csv')
f['class'] = classes
f.to_csv('dataset_S350.csv')

"""
pdb_id = ''
wild_type = ''
mutant = ''
position = 0
chain = ''
predicted_class = []
i = 0
j = 0
for record in test_file:
    if i==0:
        i = 1
        continue
    pdb_id = record[0]
    chain_id = record[0]
    pdb_id = pdb_id[0:-1]
    wild_type = record[1]
    position = record[2]
    mutant = record[3]
    chain = chain_id[-1]
    predicted_value = subprocess.check_output(['./stability_predictor.py', '%s' %pdb_id, '%s' %chain , '%s' %position , '%s' %wild_type, '%s' %mutant])
    j +=1
    predicted_class.append(predicted_value)
    print('Protein' + str(j) + 'done')

save_file = pd.read_csv('dataset_S350.csv')
save_file['predicted3'] = predicted_class
save_file.to_csv('dataset_S350.csv')
"""