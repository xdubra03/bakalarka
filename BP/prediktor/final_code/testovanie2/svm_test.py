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
train_frame = pd.read_csv(sys.argv[1])

test_file = sys.argv[2]
test_frame = pd.read_csv(test_file)

cols = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
cols1 = ['correlation','conservation','polaritychange','hydroindexchange','secondarystruc','asa','sizechange']
cols2 = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
cols3 = ['correlation','conservation','polaritychange','chargechange','secondarystruc','sizechange']

colsRes = ['class']


trainArr = train_frame.as_matrix(cols2)
trainRes = train_frame.as_matrix(colsRes)
trainRes = trainRes.ravel()

testArr = test_frame.as_matrix(cols2)
testRes = test_frame.as_matrix(colsRes)
testRes = testRes.ravel()

test_class = test_frame[['class']]
#correct = 0
classifier = svm.SVC(kernel = 'linear',class_weight={1: .5, -1: .5 })
classifier.fit(trainArr, trainRes)
results = classifier.predict(testArr)

classifier = svm.SVC(kernel = 'linear',class_weight={1: .45, -1: .55 })
classifier.fit(trainArr, trainRes)
results2 = classifier.predict(testArr)

predicted_class = results2

mcc = matthews_corrcoef(test_class, predicted_class)
print(mcc)



#with open('majority_voting4_new.csv', 'w') as f:
#	writer = csv.writer(f, delimiter=',')
#	writer.writerows(zip(results))
