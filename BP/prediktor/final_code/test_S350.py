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
i = 0
j = 0
for record in test_file:
    if i==0:
        i = 1
        continue
    pdb_id = record[0]
    wild_type = record[1]
    position = record[2]
    mutant = record[3]
    chain = record[8]
    predicted_value = subprocess.check_output(['./stability_predictor.py', '%s' %pdb_id, '%s' %chain , '%s' %position , '%s' %wild_type, '%s' %mutant])
    predicted_class.append(predicted_value)
    j +=1
    print(predicted_value)
    print('Protein' + str(j) + 'done')

save_file = pd.read_csv('dataset_S350_cleaned2_new.csv')
save_file['predicted2'] = predicted_class
save_file.to_csv('dataset_S350_cleaned2_new.csv')
