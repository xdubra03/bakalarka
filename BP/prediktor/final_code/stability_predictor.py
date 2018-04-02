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
        self.X_test = pd.read_csv("mutation_record.csv")

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

    def print_help(self):
        sys.stderr.write('Predictor usage:\n')
        sys.stderr.write('------------------')
        sys.stderr.write('\n')
        sys.stderr.write('./stability_predictor PDB_ID CHAIN POSITION WILD_TYPE MUTANT\n')
        sys.stderr.write('\n')
        sys.stderr.write('PDB_ID PDB identifier/name of protein sequence\n')
        sys.stderr.write('CHAIN chain type, e.g. A\n')
        sys.stderr.write('POSITION position of mutation in protein sequence, at firt you have to renumber position as the residue would start at position 1\n')
        sys.stderr.write('WILD_TYPE wild type amino acid\n')
        sys.stderr.write('MUTANT mutated amino acid\n')
    #parse command line arguments
    def parseArguments(self):
        if('-h' in sys.argv or '--help' in sys.argv):
             self.print_help()
             sys.exit(0)
        elif(len(sys.argv) > 6):
            sys.stderr.write("Too many arguments.\n")
            sys.exit(1)

        self.pdb_id = sys.argv[1]
        self.chain = sys.argv[2]
        self.position = int(sys.argv[3])
        self.wild_type = sys.argv[4]
        self.mutant = sys.argv[5]

        self.pdb_id = self.pdb_id.upper()
        self.chain = self.chain.upper()
        self.wild_type = self.wild_type.upper()
        self.mutant = self.mutant.upper()

    #prepare mutation record with mutation features
    def createMutationRecord(self):
        record = MutationRecord(self.pdb_id,self.chain,self.position,self.wild_type,self.mutant)
        record.create_record()
    #predict class of mutation
    def predict(self):
        ensemble = Ensemble()
        stability_class = ensemble.predict()
        return stability_class

        if(stability_class == -1):
            return "destabilizing"
        else:
            return "stabilizing"

def main():
    predictor = StabilityPredictor(object)
    predictor.parseArguments()
    predictor.createMutationRecord()
    result_class = predictor.predict()
    print("------------RESULTS------------")
    print("PDB ID: " + str(predictor.getPdb()))
    print("CHAIN TYPE: " + str(predictor.getChain()))
    print("WILD TYPE ACID: " + str(predictor.getWildType()))
    print("MUTANT ACID: " + str(predictor.getMutant()))
    print("POSITION OF MUTATION: " + str(predictor.getPosition()))
    print("TYPE OF MUTATION: " + str(result_class))

if __name__ == '__main__':
    main()
