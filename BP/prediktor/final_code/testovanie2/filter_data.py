#!/usr/bin/python

import pandas as pd
import csv


"""csv_data = open('dataset_S350.csv','r')
csv_write = open('dataset_S350_cleaned.csv','w')
dataset = csv.reader(csv_data, delimiter = ',')
filtered_data = csv.writer(csv_write, delimiter = ',')
mutations = list()


i = 0
for row in dataset:
    if i == 0:
        filtered_data.writerow(row)
        i = 1
        continue
    if float(row[6]) >= 0.5 or float(row[6]) <= -0.5:
        if(float(row[6]) >= 0.5):
            mutations.append(-1)
        else:
            mutations.append(1)
        filtered_data.writerow(row)
    else:
        continue

csv_write.close()

fdata = pd.read_csv('dataset_S350_cleaned.csv')
print(len(mutations))
fdata['class'] = mutations
fdata.to_csv('dataset_S350_cleaned.csv',index=False)
"""

def count_neg(name):
	csv_file = open(name,'r')
	data_frame = csv.reader(csv_file,delimiter = ',')
	res = list()
	countNeg = 0
	countPos = 0
	for row in data_frame:
		print(row)
		if(row[7] == '-1'):
			countNeg +=1
		else:
			countPos +=1
	print('POSITIVES: ' +str(countPos))
	print('NEGATIVES:' + str(countNeg))

count_neg('dataset_S350_cleaned.csv')