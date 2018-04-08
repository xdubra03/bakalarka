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



#names = ['train_data_cleaned_new.csv']
weights = [{1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .57, -1: .43 },
            {1: .57, -1: .43 },
            {1: .57, -1: .43 },
            {1: .4, -1: .6 },
            {1: .4, -1: .6 },
            {1: .4, -1: .6 },
            {1: .6, -1: .4 },
            {1: .6, -1: .4 },
            {1: .6, -1: .4 },
            {1: .43, -1: .57 },
            {1: .43, -1: .57 },
            {1: .43, -1: .57 },
            {1: .57, -1: .43 },
            {1: .57, -1: .43 },
            {1: .57, -1: .43 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .45, -1: .55 },
            {1: .45, -1: .55 },
            {1: .45, -1: .55 },
            {1: .55, -1: .45 },
            {1: .55, -1: .45 },
            {1: .55, -1: .45 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            {1: .5, -1: .5 },
            ]

test_file = ''
test_frame = pd.read_csv(test_file)

cols = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
cols1 = ['correlation','conservation','polaritychange','hydroindexchange','secondarystruc','asa','sizechange']
cols2 = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
cols3 = ['correlation','conservation','polaritychange','chargechange','secondarystruc','sizechange']


cm1 = ['correlation','conservation','polaritychange','chargechange']
cm2 = ['correlation','conservation','polaritychange','hydroindexchange']
cm3 = ['correlation','conservation','polaritychange','secondarystruc']
cm4 = ['correlation','conservation','polaritychange','asa']
cm5 = ['correlation','conservation','polaritychange','sizechange']

colsRes = ['class']


train_frame = pd.read_csv(sys.argv[1])
#dataset 150 150 vsetky parametre
trainArr = train_frame.as_matrix(cm1)
trainRes = train_frame.as_matrix(colsRes)
trainRes = trainRes.ravel()

testArr = test_frame.as_matrix(cm1)
testRes = test_frame.as_matrix(colsRes)
testRes = testRes.ravel()

    test_class = test_frame[['class']]
    rf = RandomForestClassifier(max_features='auto',n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
    rf.fit(trainArr,trainRes)
    result = rf.predict(testArr)

    """classifier = svm.SVC(kernel = 'poly',degree=4,class_weight='balanced')
    classifier.fit(trainArr, trainRes)
    result = classifier.predict(testArr)
    """
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
