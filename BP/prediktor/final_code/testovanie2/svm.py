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
train_frame0 = pd.read_csv("train200_250_2.csv")
train_frame1 = pd.read_csv("train250_250.csv")
train_frame2 = pd.read_csv("train250_200.csv")
train_frame3 = pd.read_csv("train300_300.csv")
train_frame4 = pd.read_csv("train300_300_1.csv")
train_frame5 = pd.read_csv("train300_300_2.csv")


test_file = sys.argv[1]
test_frame = pd.read_csv(test_file)



#test_file1 = "RF_data"
#test_frame1 = pd.read_csv(test_file)


cols = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
cols1 = ['correlation','conservation','polaritychange','hydroindexchange','secondarystruc','asa','sizechange']
cols2 = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
cols3 = ['correlation','conservation','polaritychange','chargechange','secondarystruc','sizechange']

colsRes = ['class']
#dataset 261,300 vsetky parametre
#dataset 150 150 vsetky parametre

trainArr0 = train_frame0.as_matrix(cols)
trainRes0 = train_frame0.as_matrix(colsRes)
trainRes0 = trainRes0.ravel()

#dataset 100 150 vsetky parametre
trainArr1 = train_frame1.as_matrix(cols1)
trainRes1 = train_frame1.as_matrix(colsRes)
trainRes1 = trainRes1.ravel()

#dataset 130 100 vsetky parametre
trainArr2 = train_frame2.as_matrix(cols)
trainRes2 = train_frame2.as_matrix(colsRes)
trainRes2 = trainRes2.ravel()

#dataset 150 200 bez hydroindexu
trainArr3 = train_frame3.as_matrix(cols)
trainRes3 = train_frame3.as_matrix(colsRes)
trainRes3 = trainRes3.ravel()

#datset 261 300 bez asa
trainArr4 = train_frame4.as_matrix(cols)
trainRes4 = train_frame4.as_matrix(colsRes)
trainRes4 = trainRes4.ravel()

#200 200 bez asa/hydro
trainArr5 = train_frame5.as_matrix(cols)
trainRes5 = train_frame5.as_matrix(colsRes)
trainRes5 = trainRes5.ravel()

#daatset 250 250 bez asa/hydro
trainArr6 = train_frame0.as_matrix(cols3)
trainRes6 = train_frame0.as_matrix(colsRes)
trainRes6 = trainRes6.ravel()

#daatset 250 250
trainArr8 = train_frame1.as_matrix(cols)
trainRes8 = train_frame1.as_matrix(colsRes)
trainRes8 = trainRes8.ravel()

#261 300 bez hydro/asa
trainArr7 = train_frame3.as_matrix(cols)
trainRes7 = train_frame3.as_matrix(colsRes)
trainRes7 = trainRes7.ravel()

testArr = test_frame.as_matrix(cols)
testRes = test_frame.as_matrix(colsRes)
testRes = testRes.ravel()

testArr1 = test_frame.as_matrix(cols2)
testRes1 = test_frame.as_matrix(colsRes)
testRes1 = testRes1.ravel()

testArr2 = test_frame.as_matrix(cols3)
testRes2 = test_frame.as_matrix(colsRes)
testRes2 = testRes2.ravel()

testArr3 = test_frame.as_matrix(cols1)
testRes3 = test_frame.as_matrix(colsRes)
testRes3 = testRes3.ravel()

#200 250  SVM
classifier = svm.SVC(kernel = 'linear',class_weight={1: .45, -1: .55 })
classifier.fit(trainArr0, trainRes0)
results0 = classifier.predict(testArr)
#{1: .47, -1: .53 }
#250 250 SVM
classifier = svm.SVC(kernel = 'linear',class_weight={1: .5, -1: .5 })
classifier.fit(trainArr1, trainRes1)
results1 = classifier.predict(testArr3)
#200 250 bez asa SVM
classifier = svm.SVC(kernel = 'linear',class_weight={1: .45, -1: .55 })
classifier.fit(trainArr6, trainRes6)
results2 = classifier.predict(testArr2)
#res_frame['predicted1'] = results
#res_frame.to_csv("majority_voting.csv")
#print(results)

rf = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
rf.fit(trainArr2,trainRes2)
result0 = rf.predict(testArr)

rf = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
rf.fit(trainArr3,trainRes3)
result2 = rf.predict(testArr)

rf = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
rf.fit(trainArr4,trainRes4)
result3 = rf.predict(testArr)

rf = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
rf.fit(trainArr5,trainRes5)
result4 = rf.predict(testArr)

#rf = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
#rf.fit(trainArr7,trainRes7)
#result5 = rf.predict(testArr3)

#rf = RandomForestClassifier(max_features=0.3,n_estimators=400,n_jobs=1,min_samples_leaf=50)
#rf.fit(trainArr8,trainRes8)
#result7 = rf.predict(testArr)

with open('majority_voting4_new.csv', 'w') as f:
	writer = csv.writer(f, delimiter=',')
	writer.writerows(zip(results0,results1,results2,result0,result2,result3,result4))

csv_file = open("majority_voting4_new.csv",'r')
data_frame = csv.reader(csv_file,delimiter = ',')
res = list()
countNeg = 0
countPos = 0
for row in data_frame:
	for item in row:
		if(item == '-1'):
			countNeg +=1
		else:
			countPos +=1
	if(countPos > countNeg):
		res.append(1)
	elif(countPos < countNeg):
		res.append(-1)
	else:
		res.append(1)

	countNeg = 0
	countPos = 0

frame = pd.read_csv(test_file)
frame['predicted1'] = res
frame.to_csv(test_file)

data_frame = pd.read_csv(test_file)
new_frame = data_frame[(data_frame['class'] == 1) & (data_frame['predicted1'] == 1)]
neg_frame = data_frame[(data_frame['class'] == -1) & (data_frame['predicted1'] == -1)]
#print(len(new_frame))
#print(len(neg_frame))
count = len(new_frame) + len(neg_frame)
print(count)
r = count / 250
print(str(r))


test_class = test_frame[['class']]
mcc = matthews_corrcoef(test_class, res)
print(str(mcc))

#classifier = svm.SVC(kernel = 'linear',class_weight={1: .5, -1: .5 })
#classifier.fit(trainArr0, trainRes0)
#results = classifier.predict(testArr)
#rf = RandomForestClassifier(max_features=0.3,n_estimators=400,n_jobs=1,min_samples_leaf=50)
#rf.fit(trainArr0,trainRes0)
#results = rf.predict(testArr)
#test_frame['predicted1'] = results
#test_frame.to_csv(sys.argv[1])
#print(test_frame)
