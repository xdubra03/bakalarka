import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
import csv

class Ensemble():

    def __init__(self):
        X_frame1 = pd.read_csv("train200_250_2.csv")
        X_frame2 = pd.read_csv("train250_250.csv")
        X_frame3 = pd.read_csv("train250_200.csv")
        X_frame4 = pd.read_csv("train300_300.csv")
        X_frame5 = pd.read_csv("train300_300_1.csv")
        X_frame6 = pd.read_csv("train300_300_2.csv")

    def predict(self):
        columns = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','asa','sizechange']
        columns_asa = ['correlation','conservation','polaritychange','chargechange','hydroindexchange','secondarystruc','sizechange']
        result_column = ['class']
        # prepare X_frame1
        X_trainArray0 = X_frame1.as_matrix(columns)
        Y_trainResults0 = X_frame1.as_matrix(result_column)
        Y_trainResults0 = X_trainResults0.ravel()
        #prepare X_frame2
        X_trainArray1 = X_frame2.as_matrix(columns)
        Y_trainResults1 = X_frame2.as_matrix(result_column)
        Y_trainResults1 = trainRes1.ravel()
        #prapare X_fram3
        trainArr2 = train_frame3.as_matrix(cols)
        trainRes2 = train_frame3.as_matrix(colsRes)
        trainRes2 = trainRes2.ravel()

        #dataset 150 200 bez hydroindexu
        trainArr3 = train_frame4.as_matrix(cols)
        trainRes3 = train_frame4.as_matrix(colsRes)
        trainRes3 = trainRes3.ravel()

        #datset 261 300 bez asa
        trainArr4 = train_frame5.as_matrix(cols)
        trainRes4 = train_frame5.as_matrix(colsRes)
        trainRes4 = trainRes4.ravel()

        #200 200 bez asa/hydro
        trainArr5 = train_frame6.as_matrix(cols)
        trainRes5 = train_frame6.as_matrix(colsRes)
        trainRes5 = trainRes5.ravel()

        #daatset 250 250 bez asa/hydro
        trainArr6 = train_frame1.as_matrix(cols2)
        trainRes6 = train_frame1.as_matrix(colsRes)
        trainRes6 = trainRes6.ravel()
