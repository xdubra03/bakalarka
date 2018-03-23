#!/usr/bin/python
"""main module of predictor"""

from __future__ import division

__version__ = '1.0'
__author___ = 'Juraj Ondrej Dubrava'

import os
import shutil
import sys
import subprocess
import csv
import pandas as pd
import numpy as np
import pickle

import Bio
from Bio.PDB import *
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier

from mutation_record import *
from blosum import *

class Ensemble():
    """ensemble system class"""
    def __init__(self):
        #self.X_frame1 = pd.read_csv("train200_250_2.csv")
        #self.X_frame2 = pd.read_csv("train250_250.csv")
        #self.X_frame3 = pd.read_csv("train250_200.csv")
        #self.X_frame4 = pd.read_csv("train300_300.csv")
        #self.X_frame5 = pd.read_csv("train300_300_1.csv")
        #self.X_frame6 = pd.read_csv("train300_300_2.csv")
        self.X_test = pd.read_csv("mutation_record.csv")

    def predict(self):
        """predict class of mutation"""

        columns = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
        columns_asa = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
        class_column = ['class']
        predicted_class = list()
        # prepare X_frame1
        """X_trainArray1 = self.X_frame1.as_matrix(columns)
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
        """
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
        #classifier = svm.SVC(kernel = 'linear',class_weight={1: .45, -1: .55 })
        #classifier.fit(X_trainArray1, Y_trainResults1)
        classifier = pickle.load(open('svm_model10.sav', 'rb'))
        svm_results1 = classifier.predict(X_testArray)
        predicted_class.append(svm_results1)
        
        #250 250 SVM
        #classifier = svm.SVC(kernel = 'linear',class_weight={1: .5, -1: .5 })
        #classifier.fit(X_trainArray2, Y_trainResults2)
        classifier = pickle.load(open('svm_model11.sav', 'rb'))
        svm_results2 = classifier.predict(X_testArray)
        predicted_class.append(svm_results2)

        #200 250 bez asa SVM
        classifier = pickle.load(open('svm_model12.sav', 'rb'))
        svm_results3 = classifier.predict(X_testArray)
        predicted_class.append(svm_results3)

        classifier = pickle.load(open('svm_model13.sav', 'rb'))
        svm_results4 = classifier.predict(X_testArray)
        predicted_class.append(svm_results4)

        #rf = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        #rf.fit(X_trainArray3,Y_trainResults3)
        rf = pickle.load(open('rf_model10.sav', 'rb'))
        rf_result1 = rf.predict(X_testArray)
        predicted_class.append(rf_result1)

        #rf = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        #rf.fit(X_trainArray4,Y_trainResults4)
        rf = pickle.load(open('rf_model11.sav', 'rb'))
        rf_result2 = rf.predict(X_testArray)
        predicted_class.append(rf_result2)

        #rf = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        #rf.fit(X_trainArray5,Y_trainResults5)
        rf = pickle.load(open('rf_model12.sav', 'rb'))
        rf_result3 = rf.predict(X_testArray)
        predicted_class.append(rf_result3)

        #rf = RandomForestClassifier(max_features=0.4,n_estimators=1000,n_jobs=1,min_samples_leaf=50,class_weight="balanced")
        #rf.fit(X_trainArray6,Y_trainResults6)
        #rf = pickle.load(open('rf_model4.sav', 'rb'))
        #rf_result4 = rf.predict(X_testArray)
        #predicted_class.append(rf_result4)

        #get final result by majority voting
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

        return final_class


class StabilityPredictor():
    """main predictor class, creates new predictor object"""
    def __init__(self, arg):
        self.pdb_id = ''
        self.chain = ''
        self.wild_type = ''
        self.mutant = ''
        self.position = ''

    def getPdb(self):
        return self.pdb_id
    def getChain(self):
        return self.chain
    def getWildType(self):
        return self.wild_type
    def getMutant(self):
        return self.mutant
    def getPosition(self):
        return self.position

    #parse command line arguments
    def parseArguments(self):
        if(len(sys.argv) > 6):
            sys.stderr.write("Too many arguments.\n")
            sys.exit(1)
        self.pdb_id = sys.argv[1]
        self.chain = sys.argv[2]
        self.position = int(sys.argv[3])
        self.wild_type = sys.argv[4]
        self.mutant = sys.argv[5]

    #prepare mutation record with mutation features
    def createMutationRecord(self):
        record = MutationRecord(self.pdb_id,self.chain,self.position,self.wild_type,self.mutant)
        record.create_record()
    #predict class of mutation
    def predict(self):
        ensemble = Ensemble()
        stability_class = ensemble.predict()
        return stability_class

        """if(stability_class == -1):
            return "destabilizing"
        else:
            return "stabilizing"""

def main():
    predictor = StabilityPredictor(object)
    predictor.parseArguments()
    predictor.createMutationRecord()
    result_class = predictor.predict()
    print("Result "+str(result_class))
    return result_class
    """print("------------RESULTS------------")
    print("PDB ID: " + str(predictor.getPdb()))
    print("CHAIN TYPE: " + str(predictor.getChain()))
    print("WILD TYPE ACID: " + str(predictor.getWildType()))
    print("MUTANT ACID: " + str(predictor.getMutant()))
    print("POSITION OF MUTATION: " + str(predictor.getPosition()))
    print("TYPE OF MUTATION: " + str(result_class))"""

if __name__ == '__main__':
    main()
