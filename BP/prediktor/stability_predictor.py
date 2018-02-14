#!/usr/bin/python
#imports

from __future__ import division
#from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import Bio
from Bio.PDB import *
import urllib2
import os
import shutil
import sys
import subprocess
import copy
from mutation_record import *
from blosum import *
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline


class StabilityPredictor():
    """main predictor class, creates new predictor object."""
    def __init__(self, arg):
        self.pdb_id = ''
        self.chain = ''
        self.wild_type = ''
        self.mutant = ''
        self.position = ''

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

def main():
    predictor = StabilityPredictor(object)
    predictor.parseArguments()
    predictor.createMutationRecord()

if __name__ == '__main__':
    main()
