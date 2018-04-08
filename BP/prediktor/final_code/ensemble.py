#!/usr/bin/python3.6
import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
import csv
import pickle

class Ensemble():
    """ensemble system class"""
    def __init__(self):
        self.X_test = pd.read_csv("records.csv")

    def predict(self):
        """predict class of mutation"""

        columns = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
        columns_asa = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
        class_column = ['class']
        predicted_class = list()
        #prepare frame for test record with all columns
        X_testArray = self.X_test.as_matrix(columns)
        Y_testResults = self.X_test.as_matrix(class_column)
        Y_testResults = Y_testResults.ravel()
        #prepare frame for test record without asa column
        X_testArray1 = self.X_test.as_matrix(columns_asa)
        Y_testResults1 = self.X_test.as_matrix(class_column)
        Y_testResults1 = Y_testResults1.ravel()

        #train classifiers on specific subset of dataset, 3 SVM classifiers and 4 Random Forest
        #200 250  SVM
        classifier = pickle.load(open('./models/svm_model10.sav', 'rb'))
        svm_results1 = classifier.predict(X_testArray)
        predicted_class.append(svm_results1)

        #250 250 SVM
        classifier = pickle.load(open('./models/svm_model11.sav', 'rb'))
        svm_results2 = classifier.predict(X_testArray)
        predicted_class.append(svm_results2)

        #200 250 bez asa SVM
        classifier = pickle.load(open('./models/svm_model12.sav', 'rb'))
        svm_results3 = classifier.predict(X_testArray)
        predicted_class.append(svm_results3)

        classifier = pickle.load(open('./models/svm_model13.sav', 'rb'))
        svm_results4 = classifier.predict(X_testArray)
        predicted_class.append(svm_results4)

        rf = pickle.load(open('./models/rf_model10.sav', 'rb'))
        rf_result1 = rf.predict(X_testArray)
        predicted_class.append(rf_result1)

        rf = pickle.load(open('./models/rf_model11.sav', 'rb'))
        rf_result2 = rf.predict(X_testArray)
        predicted_class.append(rf_result2)

        rf = pickle.load(open('./models/rf_model12.sav', 'rb'))
        rf_result3 = rf.predict(X_testArray)
        predicted_class.append(rf_result3)

        with open('majority_voting4_new.csv', 'w') as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerows(zip(svm_results1,svm_results2,svm_results3,svm_results4,rf_result1,rf_result2,rf_result3))

e = Ensemble()
e.predict()