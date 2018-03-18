#!/usr/bin/python
from __future__ import division
import numpy as np
import pandas as pd
import sys
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
import csv
from sklearn.metrics import matthews_corrcoef
from sklearn.model_selection import train_test_split


# train_file = sys.argv[1]
train_frame = pd.read_csv(sys.argv[1])

test_file = sys.argv[2]
test_frame = pd.read_csv(test_file)

cols = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
cols1 = ['correlation','conservation','polaritychange','hydroindexchange','secondarystruc','asa','sizechange']
cols2 = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
cols3 = ['correlation','conservation','polaritychange','chargechange','secondarystruc','sizechange']

colsRes = ['class']

X = train_frame.iloc[:,0:7].values
y = train_frame.iloc[:, 8].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.17, random_state = 0)

trainArr = train_frame.as_matrix(cols)
trainRes = train_frame.as_matrix(colsRes)
trainRes = trainRes.ravel()

testArr = test_frame.as_matrix(cols)
testRes = test_frame.as_matrix(colsRes)
testRes = testRes.ravel()

test_class = test_frame[['class']]
#correct = 0
classifier = svm.SVC(kernel = 'poly',class_weight='balanced')
classifier.fit(trainArr, trainRes)
results = classifier.predict(testArr)

classifier = svm.SVC(kernel = 'linear',class_weight={1: .45, -1: .55 })
classifier.fit(trainArr, trainRes)
results2 = classifier.predict(testArr)

predicted_class = results

mcc = matthews_corrcoef(test_class, predicted_class)
print(mcc)



#with open('majority_voting4_new.csv', 'w') as f:
#	writer = csv.writer(f, delimiter=',')
#	writer.writerows(zip(results))
