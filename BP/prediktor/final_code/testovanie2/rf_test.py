#!/usr/bin/python
from __future__ import division
import numpy as np
import pandas as pd
import sys
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
import csv
from sklearn.metrics import matthews_corrcoef

#train_file = sys.argv[1]
#train_frame = pd.read_csv(sys.argv[1])

names = ['train100_100.csv',
'train100_100_1.csv',
'train100_100_2.csv',
'train100_150.csv',
'train100_150_1.csv',
'train100_150_2.csv',
'train130_100.csv',
'train130_100_1.csv',
'train130_100_2.csv',
'train150_100.csv',
'train150_100_1.csv',
'train150_100_2.csv',
'train150_150.csv',
'train150_150_1.csv',
'train150_150_2.csv',
'train150_200.csv',
'train150_200_1.csv',
'train150_200_2.csv',
'train200_150.csv',
'train200_150_1.csv',
'train200_150_2.csv',
'train200_200.csv',
'train200_200_1.csv',
'train200_200_2.csv',
'train200_250.csv',
'train200_250_1.csv',
'train200_250_2.csv',
'train250_200.csv',
'train250_200_1.csv',
'train250_200_2.csv',
'train250_250.csv',
'train250_250_1.csv',
'train250_250_2.csv',
'train300_300.csv',
'train300_300_1.csv',
'train300_300_2.csv']

test_file = sys.argv[1]
test_frame = pd.read_csv(test_file)

cols = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
cols1 = ['correlation','conservation','polaritychange','hydroindexchange','secondarystruc','asa','sizechange']
cols2 = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
cols3 = ['correlation','conservation','polaritychange','chargechange','secondarystruc','sizechange']

colsRes = ['class']

for name in names:
    train_frame = pd.read_csv(name)
    #dataset 150 150 vsetky parametre
    trainArr = train_frame.as_matrix(cols3)
    trainRes = train_frame.as_matrix(colsRes)
    trainRes = trainRes.ravel()

    testArr = test_frame.as_matrix(cols3)
    testRes = test_frame.as_matrix(colsRes)
    testRes = testRes.ravel()

    test_class = test_frame[['class']]

    rf = RandomForestClassifier(max_features=0.3,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
    rf.fit(trainArr,trainRes)
    result = rf.predict(testArr)

    predicted_class = result

    mcc = matthews_corrcoef(test_class, predicted_class)
    print(name +": "+ str(mcc))

#with open('majority_voting4_new.csv', 'w') as f:
#	writer = csv.writer(f, delimiter=',')
#	writer.writerows(zip(result))

#classifier = svm.SVC(kernel = 'linear',class_weight={1: .5, -1: .5 })
#classifier.fit(trainArr0, trainRes0)
#results = classifier.predict(testArr)
#rf = RandomForestClassifier(max_features=0.3,n_estimators=400,n_jobs=1,min_samples_leaf=50)
#rf.fit(trainArr0,trainRes0)
#results = rf.predict(testArr)
#test_frame['predicted1'] = results
#test_frame.to_csv(sys.argv[1])
#print(test_frame)
