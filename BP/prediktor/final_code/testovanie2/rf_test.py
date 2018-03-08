#!/usr/bin/python3.6
from __future__ import division
import numpy as np
import pandas as pd
import sys
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
import csv
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import classification_report
#train_file = sys.argv[1]
#train_frame = pd.read_csv(sys.argv[1])

names = ['train100_100.csv',
'train100_100_1.csv',
'train100_100_2.csv',
'train150_150.csv',
'train150_150_1.csv',
'train150_150_2.csv',
'train130_100.csv',
'train130_100_1.csv',
'train130_100_2.csv',
'train100_150.csv',
'train100_150_1.csv',
'train100_150_2.csv',
'train150_100.csv',
'train150_100_1.csv',
'train150_100_2.csv',
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

#names = ['train_data_cleaned_new.csv']

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
    trainArr = train_frame.as_matrix(cols1)
    trainRes = train_frame.as_matrix(colsRes)
    trainRes = trainRes.ravel()

    testArr = test_frame.as_matrix(cols1)
    testRes = test_frame.as_matrix(colsRes)
    testRes = testRes.ravel()

    test_class = test_frame[['class']]
    rf = RandomForestClassifier(max_features=0.3,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
    rf.fit(trainArr,trainRes)
    result = rf.predict(testArr)

    """classifier = svm.SVC(kernel = 'linear',class_weight={1: .2, -1: .8 })
    classifier.fit(trainArr, trainRes)
    result = classifier.predict(testArr)
    """
    predicted_class = result
    #predicted_class = list(predicted_class)
    #predicted_class.remove(-1)
   # predicted_class.append(1)
    #print(predicted_class)
    #frame = pd.read_csv(test_file)
    #frame['predicted1'] = predicted_class
    #frame.to_csv(test_file)
    #data_frame = pd.read_csv(test_file)
    #new_frame = data_frame[(data_frame['class'] == 1) & (data_frame['predicted1'] == 1)]
    #neg_frame = data_frame[(data_frame['class'] == -1) & (data_frame['predicted1'] == -1)]
	#print(len(new_frame))
	#print(len(neg_frame))
    #count = len(new_frame) + len(neg_frame)
    #print(count)
    #r = count / 250
    #print(str(r))
    #testarr = data_frame['class'].values
    #trainarr = data_frame['predicted1'].values
    #target_names = ['1', '-1']
    #print(classification_report(testarr,trainarr, target_names=target_names))
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
