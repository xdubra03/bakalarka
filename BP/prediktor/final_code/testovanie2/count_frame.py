#!/usr/bin/python

from __future__ import division
import numpy as np
import pandas as pd
import random
import sys
import csv
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import classification_report

def count_frame(file):
	data_frame = pd.read_csv(file)
	data_frame_len = len(data_frame)
	res_frame = data_frame[data_frame['class'] == data_frame['predicted1']]
	count = len(res_frame)
	acc = count / data_frame_len
	print("Accuracy: " + str(acc))

def compute_mcc(file):
	data_frame = pd.read_csv(file)
	testarr = data_frame['class'].values
	trainarr = data_frame['predicted1'].values
	mcc = matthews_corrcoef(testarr, trainarr)
	target_names = ['-1', '1']
	print(classification_report(testarr,trainarr, target_names=target_names))
	print('Correlation: ' +str(mcc))


compute_mcc(sys.argv[1])
count_frame(sys.argv[1])
