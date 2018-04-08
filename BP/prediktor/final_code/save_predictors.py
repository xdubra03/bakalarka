#!/usr/bin/python
import pandas as pd
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
import pickle

class Ensemble():
    """ensemble system class"""
    def __init__(self):
        self.X_frame1 = pd.read_csv("train200_200_2.csv")
        self.X_frame2 = pd.read_csv("train200_250_2.csv")
        self.X_frame3 = pd.read_csv("train250_200_2.csv")
        self.X_frame4 = pd.read_csv("train250_250_2.csv")
        self.X_frame5 = pd.read_csv("train250_200.csv")
        self.X_frame6 = pd.read_csv("train250_250_1.csv")

    def create_models(self):
        """create models of predictors"""

        columns = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
        columns_asa = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
        class_column = ['class']
        # prepare X_frame1
        x_trainArray1 = self.X_frame1.as_matrix(columns_asa)
        X_trainResults1 = self.X_frame1.as_matrix(class_column)
        Y_trainResults1 = Y_trainResults1.ravel()
        # prepare X_frame2
        X_trainArray2 = self.X_frame2.as_matrix(columns)
        Y_trainResults2 = self.X_frame2.as_matrix(class_column)
        Y_trainResults2 = Y_trainResults2.ravel()
        # prapare X_frame3
        X_trainArray3 = self.X_frame3.as_matrix(columns)
        Y_trainResults3 = self.X_frame3.as_matrix(class_column)
        Y_trainResults3 = Y_trainResults3.ravel()
        # prepare X_frame4
        X_trainArray4 = self.X_frame4.as_matrix(columns)
        Y_trainResults4 = self.X_frame4.as_matrix(class_column)
        Y_trainResults4 = Y_trainResults4.ravel()
        # prepare X_frame5
        X_trainArray5 = self.X_frame5.as_matrix(columns)
        Y_trainResults5 = self.X_frame5.as_matrix(class_column)
        Y_trainResults5 = Y_trainResults4.ravel()
        # prepare X_frame6
        X_trainArray6 = self.X_frame6.as_matrix(columns)
        Y_trainResults6 = self.X_frame6.as_matrix(class_column)
        Y_trainResults6 = Y_trainResults6.ravel()

        X_trainArray7 = self.X_frame3.as_matrix(column)
        Y_trainResults7 = self.X_frame3.as_matrix(class_column)
        Y_trainRes7 = Y_trainResults7.ravel()


        # train classifiers on specific subset of dataset, 3 SVM classifiers and 4 Random Forest
        # 200 250  SVM
        classifier1 = svm.SVC(kernel = 'poly',degree=4,class_weight='balanced')
        classifier1.fit(X_trainArray1, Y_trainResults1)
        filename = 'svm_model1.sav'
        pickle.dump(classifier1, open(filename, 'wb'))

        classifier2 = svm.SVC(kernel = 'poly',degree=4,class_weight='balanced')
        classifier2.fit(X_trainArray2, Y_trainResults2)
        filename = 'svm_model2.sav'
        pickle.dump(classifier2, open(filename, 'wb'))

        classifier3 = svm.SVC(kernel = 'poly',degree=4,class_weight='balanced')
        classifier3.fit(X_trainArray3, Y_trainResults3)
        filename = 'svm_model3.sav'
        pickle.dump(classifier3, open(filename, 'wb'))

        classifier4 = svm.SVC(kernel='poly', degree=4, class_weight='balanced')
        classifier4.fit(X_trainArray4, Y_trainResults4)
        filename = 'svm_model4.sav'
        pickle.dump(classifier4, open(filename, 'wb'))

        rf1 = RandomForestClassifier(max_features='auto',n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        rf1.fit(X_trainArray5,Y_trainResults5)
        filename = 'rf_model1.sav'
        pickle.dump(rf1, open(filename, 'wb'))

        rf2 = RandomForestClassifier(max_features='auto',n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        rf2.fit(X_trainArray3,Y_trainResults3)
        filename = 'rf_model2.sav'
        pickle.dump(rf2, open(filename, 'wb'))

        rf3 = RandomForestClassifier(max_features='auto',n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        rf3.fit(X_trainArray6,Y_trainResults6)
        filename = 'rf_model3.sav'
        pickle.dump(rf3, open(filename, 'wb'))



p = Ensemble()
p.create_models()