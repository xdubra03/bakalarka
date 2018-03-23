#!/usr/bin/python
import pandas as pd
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import pickle

class Ensemble():
    """ensemble system class"""
    def __init__(self):
        self.X_frame1 = pd.read_csv("train200_250_2.csv")
        self.X_frame2 = pd.read_csv("train250_250.csv")
        self.X_frame3 = pd.read_csv("train250_200.csv")
        self.X_frame4 = pd.read_csv("train300_300.csv")
        self.X_frame5 = pd.read_csv("train300_300_1.csv")
        self.X_frame6 = pd.read_csv("train300_300_2.csv")
        #self.X_test = pd.read_csv("mutation_record.csv")

    def predict(self):
        """predict class of mutation"""

        columns = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
        columns_asa = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
        class_column = ['class']
        predicted_class = list()
        # prepare X_frame1
        X_trainArray1 = self.X_frame1.as_matrix(columns)
        Y_trainResults1 = self.X_frame1.as_matrix(class_column)
        Y_trainResults1 = Y_trainResults1.ravel()
        #prepare X_frame2
        X_trainArray2 = self.X_frame2.as_matrix(columns)
        Y_trainResults2 = self.X_frame2.as_matrix(class_column)
        Y_trainResults2 = Y_trainResults2.ravel()
        #prapare X_frame3
        X_trainArray3 = self.X_frame3.as_matrix(columns)
        Y_trainResults3 = self.X_frame3.as_matrix(class_column)
        Y_trainResults3 = Y_trainResults3.ravel()
        #prepare X_frame4
        X_trainArray4 = self.X_frame4.as_matrix(columns)
        Y_trainResults4 = self.X_frame4.as_matrix(class_column)
        Y_trainResults4 = Y_trainResults4.ravel()
        #prepare X_frame5
        X_trainArray5 = self.X_frame5.as_matrix(columns)
        Y_trainResults5 = self.X_frame5.as_matrix(class_column)
        Y_trainResults5 = Y_trainResults4.ravel()
        #prepare X_frame6
        X_trainArray6 = self.X_frame6.as_matrix(columns)
        Y_trainResults6 = self.X_frame6.as_matrix(class_column)
        Y_trainResults6 = Y_trainResults6.ravel()
        #prepare X_frame1 without asa
        X_trainArray7 = self.X_frame1.as_matrix(columns_asa)
        Y_trainResults7 = self.X_frame1.as_matrix(class_column)
        Y_trainRes7 = Y_trainResults7.ravel()
        #prepare frame for test record with all columns
        #X_testArray = self.X_test.as_matrix(columns)
        #Y_testResults = self.X_test.as_matrix(class_column)
        #Y_testResults = Y_testResults.ravel()
        #prepare frame for test record without asa column
        #X_testArray1 = self.X_test.as_matrix(columns_asa)
        #Y_testResults1 = self.X_test.as_matrix(class_column)
        #Y_testResults1 = Y_testResults1.ravel()

        #train classifiers on specific subset of dataset, 3 SVM classifiers and 4 Random Forest
        #200 250  SVM
        classifier1 = svm.SVC(kernel = 'poly',degree=4,class_weight={1: .5, -1: .5 })
        classifier1.fit(X_trainArray4, Y_trainResults4)
        filename = 'svm_model4.sav'
        pickle.dump(classifier1, open(filename, 'wb'))

        #svm_results1 = classifier1.predict(X_testArray)
        #predicted_class.append(svm_results1)
        #250 250 SVM
        classifier2 = svm.SVC(kernel = 'poly',degree=4,class_weight={1: .5, -1: .5 })
        classifier2.fit(X_trainArray5, Y_trainResults5)
        filename = 'svm_model5.sav'
        pickle.dump(classifier2, open(filename, 'wb'))


        #svm_results2 = classifier2.predict(X_testArray)
        #predicted_class.append(svm_results2)
        #200 250 bez asa SVM
        classifier3 = svm.SVC(kernel = 'poly',degree=4,class_weight={1: .5, -1: .5 })
        classifier3.fit(X_trainArray6, Y_trainResults6)
        filename = 'svm_model6.sav'
        pickle.dump(classifier3, open(filename, 'wb'))

        #svm_results3 = classifier3.predict(X_testArray1)
        #predicted_class.append(svm_results3)

        #rf1 = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        #rf1.fit(X_trainArray3,Y_trainResults3)
        #ilename = 'rf_model1.sav'
        #pickle.dump(rf1, open(filename, 'wb'))

        #rf_result1 = rf1.predict(X_testArray)
        #predicted_class.append(rf_result1)

        #rf2 = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        #rf2.fit(X_trainArray4,Y_trainResults4)
        #filename = 'rf_model2.sav'
        #pickle.dump(rf2, open(filename, 'wb'))

        #rf_result2 = rf2.predict(X_testArray)
        #predicted_class.append(rf_result2)

        ##rf3 = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        #rf3.fit(X_trainArray5,Y_trainResults5)
        #filename = 'rf_model3.sav'
        #pickle.dump(rf3, open(filename, 'wb'))

        #rf_result3 = rf3.predict(X_testArray)
        #predicted_class.append(rf_result3)

        #rf4 = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        #rf4.fit(X_trainArray6,Y_trainResults6)
        #filename = 'rf_model4.sav'
        #pickle.dump(rf4, open(filename, 'wb'))

        #rf_result4 = rf4.predict(X_testArray)
        #predicted_class.append(rf_result4)

        """""#get final result by majority voting
        countNeg = 0
        countPos = 0
        final_class = 1
        for item in predicted_class:
            if(item == -1):
                countNeg +=1
            else:
                countPos +=1
        if(countPos > countNeg):
            final_class = 1
        elif(countPos < countNeg):
            final_class = -1

        return final_class"""

p = Ensemble()
p.predict()