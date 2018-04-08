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


test_file = 'records.csv'
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
cm6 = ['correlation','conservation','chargechange','secondarystruc']
cm7 = ['correlation','conservation','hydroindexchange','secondarystruc','asa']
cm8 = ['conservation','polaritychange','hydroindexchange','secondarystruc','asa']
cm9 = ['chargechange']
cm10 = ['hydroindexchange']
cm11 = ['secondarystruc']
cm12 = ['asa']
cm13 = ['sizechange']


colsRes = ['class']

test_files = ['corr_cons_polar_charge_class.csv ',
              'corr_cons_polar_hydro_class.csv',
              'corr_cons_polar_secondary_class.csv',
              'corr_cons_polar_asa_class.csv',
              'corr_cons_polar_size_class.csv',
              ]

#dataset 150 150 vsetky parametre
train_file = 'train_data_cleaned_new.csv'
train_frame = pd.read_csv(train_file)
################################
trainArr = train_frame.as_matrix(cols)
trainRes = train_frame.as_matrix(colsRes)
trainRes = trainRes.ravel()

testArr = test_frame.as_matrix(cols)
testRes = test_frame.as_matrix(colsRes)
testRes = testRes.ravel()
######################################

# na SVM aj RF
###################################
trainArr0 = train_frame.as_matrix(cols3)
trainRes0 = train_frame.as_matrix(colsRes)
trainRes0 = trainRes0.ravel()

testArr0 = test_frame.as_matrix(cols3)
testRes0 = test_frame.as_matrix(colsRes)
testRes0 = testRes0.ravel()
##########################################
# na SVM aj RF vsetky parametre
#########################################
trainArr1 = train_frame.as_matrix(cols)
trainRes1 = train_frame.as_matrix(colsRes)
trainRes1 = trainRes1.ravel()

testArr1 = test_frame.as_matrix(cols)
testRes1 = test_frame.as_matrix(colsRes)
testRes1 = testRes1.ravel()
##########################################
# na SVM
#########################################

trainArr2 = train_frame.as_matrix(cm6)
trainRes2 = train_frame.as_matrix(colsRes)
trainRes2= trainRes2.ravel()

testArr2 = test_frame.as_matrix(cm6)
testRes2 = test_frame.as_matrix(colsRes)
testRes2 = testRes2.ravel()
##########################################
# na SVM
##########################################

trainArr3 = train_frame.as_matrix(cm5)
trainRes3 = train_frame.as_matrix(colsRes)
trainRes3 = trainRes3.ravel()

testArr3 = test_frame.as_matrix(cm5)
testRes3 = test_frame.as_matrix(colsRes)
testRes3 = testRes3.ravel()
##########################################
# na RF
##########################################

trainArr4 = train_frame.as_matrix(cols1)
trainRes4 = train_frame.as_matrix(colsRes)
trainRes4 = trainRes4.ravel()

testArr4 = test_frame.as_matrix(cols1)
testRes4 = test_frame.as_matrix(colsRes)
testRes4 = testRes4.ravel()
########################################
trainArr5 = train_frame.as_matrix(cm5)
trainRes5 = train_frame.as_matrix(colsRes)
trainRes5 = trainRes5.ravel()

testArr5 = test_frame.as_matrix(cm5)
testRes5 = test_frame.as_matrix(colsRes)
testRes5 = testRes5.ravel()


test_class = test_frame[['class']]


#rf = RandomForestClassifier(max_features='auto',n_estimators=100,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
#rf.fit(trainArr,trainRes)
#result = rf.predict(testArr)

classifier = svm.SVC(kernel = 'linear',class_weight='balanced')
classifier.fit(trainArr, trainRes)
result = classifier.predict(testArr)



"""
rf = RandomForestClassifier(max_features='auto',n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
rf.fit(trainArr0,trainRes0)
result4 = rf.predict(testArr0)

rf = RandomForestClassifier(max_features='auto',n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
rf.fit(trainArr1,trainRes1)
result5 = rf.predict(testArr1)

rf = RandomForestClassifier(max_features='auto',n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
rf.fit(trainArr4,trainRes4)
result6 = rf.predict(testArr4)

classifier = svm.SVC(kernel = 'rbf',degree=3,class_weight='balanced')
classifier.fit(trainArr0, trainRes0)
result0 = classifier.predict(testArr0)


classifier = svm.SVC(kernel = 'poly',degree=4,class_weight='balanced')
classifier.fit(trainArr1, trainRes1)
result1 = classifier.predict(testArr1)


classifier = svm.SVC(kernel = 'poly',degree=4,class_weight='balanced')
classifier.fit(trainArr2, trainRes2)
result2 = classifier.predict(testArr2)


classifier = svm.SVC(kernel = 'rbf',degree=3,class_weight='balanced')
classifier.fit(trainArr3, trainRes3)
result3 = classifier.predict(testArr3)

"""


predicted_class = result
mcc = matthews_corrcoef(test_class, predicted_class)
print("File : "+ str(mcc))


#with open('majority_voting4_new.csv', 'w') as f:
#        writer = csv.writer(f, delimiter=',')
#        writer.writerows(zip(result0,result1,result2,result3,result4,result5,result6))

#zatial korelacia 0.476