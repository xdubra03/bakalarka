#!/usr/bin/python3.6

import pandas as pd
import csv

class_column = []

file1 = csv.reader(open('dataset_S350.csv'))
for line in file1:
    class_column.append(line[13])

class_column.pop(0)

file2 = pd.read_csv('records.csv')
file2['class'] = class_column
file2.to_csv('records.csv')
