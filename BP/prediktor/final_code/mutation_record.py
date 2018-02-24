#this is module for preparing mutation record, e.g. all features in record

from __future__ import division
import math
import pandas as pd
import Bio
from Bio.PDB import *
import urllib2
import os
import shutil
import sys
import subprocess
from ete3 import Tree
from blosum import *
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline

class Etree(Tree):
    """class for creating phylogenetic tree and computing conservation of mutated position"""
    _names = []
    alignements = dict()
    _identificators = []
    _IDs = dict()
    _idArray = dict()

    #_matrix = BlosumMatrix('./tools/blosum62.txt')
    #TODO predat meno do funkcie
    #clustal_file = sys.argv[1]
    #_probability_matrix = ProbabilityMatrix(clustal_file+'clustal.fasta',_matrix)

     #parse PDB file to find stating position
    def PDB_parse(self,name):
		p = PDBParser()
		structure = p.get_structure(name,name+".pdb")
		model = structure[0]
		#pridat try na jednotlive chainy
		try:
			chain = model['A']
		except KeyError as error:
			try:
				chain = model['B']
			except KeyError as error:
				try:
					chain = model['C']
				except KeyError as error:
					try:
						chain = model['I']
					except KeyError as error:
						try:
							chain = model['X']
						except KeyError as error:
							print("Cannot find this type of chain.")
							sys.exit(1)
						else:
							pass
					else:
						pass
				else:
					pass
			else:
				pass
		else:
			pass
		#always returns position of first chain which could no be correct
		residue_list = Selection.unfold_entities(chain,'A')
		#print(residue_list[0].get_full_id()[3][1])
		residue_start = residue_list[0].get_full_id()[3][1]
		return residue_start

    #compute conservation score on entered position
    def compute_conservation(self,file,residue_start,index,weightsArray,acid1):
		count_mezera = 0
		count_basic_acid = 0
		count_mutated_acid = 0
		count_else=0
		all_count =0
		count_pos=0
		start_position = 1 #meni sa
		pos = 0
		handle = open(file,"r")

		lines = iter(handle.readlines())
		for line in lines:
			if(line.startswith('>')):
				continue
			else:
				for word in line.split():
	 			#if(word[0] == '-'):
	 			#	break
		 			#if(word[0] == 'M'):
		 			#		count_pos -=1#-residue_start+1

		 			if(residue_start > len(word)):
		 				print(residue_start)
		 				print(index)
		 				count_pos = residue_start
		 				print(count_pos)
		 				for i in range(0,len(word),1):
		 					if(word[i] != '-'):
		 						count_pos +=1
		 						if(count_pos == residue_start+index):
		 							pos = i
		 							print(word[i])
		 							break
		 			else:
		 				print(residue_start)
		 				print(index)
		 				count_pos = residue_start
		 				if(residue_start < 0):
		 					chain_res = index#+residue_start #+ abs(residue_start) + abs(residue_start) -1
		 				elif (residue_start == 1):
		 					chain_res= index+residue_start
		 				else:
		 					chain_res= index+residue_start+2

		 				for i in range(0,len(word),1):
		 					if(word[i] != '-'):
		 						count_pos +=1
		 						if(count_pos == chain_res):
		 							pos = i
		 							#print(pos)
		 							print(word[i])
		 							break
	 		break
	 	print(pos)
	 	conservation_value = 0
	 	base_acid = 0
	 	weights = 0
	 	for name in self._names:
	 		sequence = self._idArray[name]
	 		acid = sequence[pos]
	 		if(acid == acid1):
	 			base_acid = 1
	 		else:
	 			base_acid= 0
	 		weights += weightsArray[name]
	 		conservation_value += weightsArray[name] * base_acid

	 	accuracy = conservation_value/ weights
	 	return accuracy

    def create_ID_table(self):
		"""create table where key is node name and value is sequence to speed up lookup"""
		for name in self._names:
			key1 = self._IDs.get(name)
			seq1 = self.alignements[key1]
			self._idArray[name] = seq1

    def create_alignement_table(self,file):
		"""creates lookup table for sequence names and sequences"""
		with open(file,'r') as f:
			lines = iter(f.readlines())
			for line in lines:
				if(line.startswith('>')):
					name = line.strip('>').strip('\n')
					sequence = lines.next().strip('\n')
					self.alignements[name] = sequence


    def create_names_table(self,file):
		"""create lookup table for complete sequence ID according to its abbrevation"""
		with open(file,'r') as f:
			lines = iter(f.readlines())
			for line in lines:
				if(line.startswith('>')):
					self._identificators.append(line.strip('>').strip('\n'))

		for item in self._identificators:
			for name in self._names:
				if(name in item):
					self._IDs[name] = item


    def get_table_value(self,value):
		"""get value from alignements table"""
		return self.alignements[value]


    def get_names(self):
		"""get all leaf names in the tree and stores them in _names array"""
		for leaf in self:
			if(leaf.is_leaf()):
				self._names.append(leaf.name)

    def print_names(self):
		"""function for printing leafs names"""
		for name in self._names:
			print(name)

    def create_array(self):
		"""creates array of weights and fills it with value according to its node"""
		self.weightsArray = dict()
		for name in self._names:
			self.weightsArray[name] = 0
		if self.name != '':
			self.weightsArray[self.name] = 1

    def add_node_array(self):
		"""adds weights array to every node in the tree"""
		for node in self.traverse('postorder'):
			node.create_array()

    def calculate_weights(self,probability_matrix):
		"""calculates the values in weights array in each node"""

		#fudge factor constant to prevent 0 in the weights array
		fugde_factor = 0.1

		#traverse the tree and compute values in each node
		for node in self.traverse('postorder'):

			#get children nodes of actual node
			children = node.get_children()

			#if no children found, continue with next node
			if not children:
				continue
			else:

				i = 0
				#array where value of multiplication for each item in array is stored
				vals = [1]*250

				#calculate value for each child
				for child in children:

					for parentItem in node._names:
						result = 0
						seq2 = node._idArray[parentItem]

						for childItem in child._names:

							#calculate probability of changing child sequence to parent sequence
							seq1 = child._idArray[childItem]
							probability = probability_matrix.find_pair(seq1,seq2)

							#vzorec Pi*Li*t
							result += probability * child.weightsArray[childItem] * (child.dist + fugde_factor)

						#value from each child needs to be multiplicated
						vals[i] *= result
						#store actual value to weightsArray item in parent node
						node.weightsArray[parentItem] = vals[i]

						i+=1
					i = 0

			#print(node.weightsArray.values())

		return self.get_tree_root().weightsArray


class MutationRecord():
    """class for computing mutation record features"""

    """encoding secondary structure"""
    __encode_sec_struc = {
    	'C':'1','H':'2','S':'3'
    }

    """encoding conservation"""
    __encode_conservation = {
    	'+':'1','-':'-1'
    }

    """encoding charge"""
    __encode_charge = {
    	'+':'1','-':'-1','0':'0'
    }

    """encoding polarity"""
    __encode_polarity = {
    	'+':'1','-':'-1','0':'0'
    }

    """encoding accessible surface area"""
    __encode_asa = {
    	'burried':'10','partially_burried':'11','exposed':'12'
    }

    """encoding size"""
    __encode_size = {
    	'+':'1','-':'-1','0':'0'
    }


    """dictionary for polarity of amino acids
     +  polar
     -  nonpolar
     """
    __polarity = {'A':'-','R':'+','N':'+','D':'+','C':'-','E':'+',
                'Q':'+','G':'-','H':'+','I':'-','L':'-','K':'+',
                'M':'-','F':'-','P':'-','S':'+','T':'+','W':'-',
                'Y':'+','V':'-'
    }
    """dictionary for charge of amino acids
     0 neutral
     + positive
     - negative
    """
    __charge = {'A':'0','R':'+','N':'0','D':'-','C':'0','E':'-',
                'Q':'0','G':'0','H':'0','I':'0','L':'0','K':'+',
                'M':'0','F':'0','P':'0','S':'0','T':'0','W':'0',
                'Y':'0','V':'0'
    }

    """dictionary for hydropathy index of amino acid"""
    __hydro_index = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,
                   'E':-3.5,'Q':-3.5,'G':-0.4,'H':-3.2,'I':4.5,
                   'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,
                   'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2
    }

    """dictionary for size change of amino acids"""
    __size_change1 = {'G':75.1,'A':89.1,'S':105.1,'P':115.1, #first interval
    				'V':117.1,'T':119.1,'C':121.2,'I':131.2,'L':131.2,'N':132.1,'D':133.1,
    				'Q':146.1,'K':146.2,'E':147.1,'M':149.2,'H':155.2, #second interval
    				'F':165.2,'R':174.2,'Y':181.2, #3rd interval
    				'W':204.2 #last interval
    }

    def __init__(self,pdb_id,chain,position,wild_type,mutant):
        self.pdb_id = pdb_id
        self.chain = chain
        self.position = position
        self.wild_type = wild_type
        self.mutant = mutant

    #get pdb file for entered pdb id
    def get_pdb(self):
		protein_file = self.pdb_id + ".pdb"
		pdb_output  = open(protein_file, "w")
		url_pdb = "https://files.rcsb.org/download/%s.pdb" %self.pdb_id
		try:
			handle = urllib2.urlopen(url_pdb)
		except URLError as error:
			print(error.reason)
			sys.exit(1)
		pdb_output.write(handle.read())
		pdb_output.close()

    def PDB_parse(self,name):
		p = PDBParser()
		structure = p.get_structure(name,name+".pdb")
		model = structure[0]
		#pridat try na jednotlive chainy
		try:
			chain = model['A']
		except KeyError as error:
			try:
				chain = model['B']
			except KeyError as error:
				try:
					chain = model['C']
				except KeyError as error:
					try:
						chain = model['I']
					except KeyError as error:
						try:
							chain = model['X']
						except KeyError as error:
							print("Cannot find this type of chain.")
							sys.exit(1)
						else:
							pass
					else:
						pass
				else:
					pass
			else:
				pass
		else:
			pass
		#always returns position of first chain which could no be correct
		residue_list = Selection.unfold_entities(chain,'A')
		#print(residue_list[0].get_full_id()[3][1])
		residue_start = residue_list[0].get_full_id()[3][1]
		return residue_start

    def get_fasta(self):
    	fasta_file = self.pdb_id + ".fasta"
    	fasta_output = open(fasta_file,"w")
        url_f = "https://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=fastachain&compression=NO&structureId=%s&chainId=%s"%(self.pdb_id,self.chain)
    	#url_f = "http://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=%s&compressionType=uncompressed"%self.name
    	try:
    		h = urllib2.urlopen(url_f)
    	except urllib2.URLError as error:
    		print(error.reason)
    		sys.exit(1)
    	fasta_output.write(h.read())
    	fasta_output.close()

    #parse XML file from BLASTP output and create file with sequences
    def parse_XML(self):
    	i=0
    	output = open(self.pdb_id+".txt","w")
    	f = open(self.pdb_id+".xml","r")
    	blast = NCBIXML.parse(f)
    	names = []
    	protein = ''
    	for record in blast:
    		for align in record.alignments:
    			for hsp in align.hsps:
    				i+= 1
    				protein = '>'+align.hit_id+align.hit_def
    				if(protein in names):
    					break
    				else:
    					names.append(protein)
    					output.write('>'+align.hit_id+align.hit_def+'\n')
    				output.write(hsp.sbjct+'\n') #find out
    	f.close()
    	output.close()

    #run standalone BLASTP for searching homologue sequences
    def run_blast(self):
        FNULL = open(os.devnull, 'w')
    	subprocess.call(['./tools/blastp','-query','%s.fasta'%self.pdb_id,'-db','nr','-outfmt','5','-out','%s.xml'%self.pdb_id,'-max_target_seqs','250','-remote'],stdout=FNULL, stderr=subprocess.STDOUT)
        FNULL.close()
    #run CLUSTAL OMEGA to create multiple sequence alignement
    def run_clustal(self):
    	in_file = self.pdb_id + ".txt"
    	out_file = self.pdb_id + "clustal.fasta"
        FNULL = open(os.devnull, 'w')
    	clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
    	subprocess.call(['./tools/clustalo','-i','%s'%in_file,'--outfmt=vie','-o','%s'%out_file,'--auto','-v','--force'],stdout=FNULL, stderr=subprocess.STDOUT)
        FNULL.close()

    #create newick file with phylogenetic tree for conservation
    def create_newick_file(self):
        self.run_blast()
    	self.parse_XML()
    	self.run_clustal()
    	protein = self.pdb_id+'clustal.fasta'
    	output = self.pdb_id+'_newick'
        FNULL = open(os.devnull, 'w')
        subprocess.call(['./tools/FastTree','-out','%s'%output,'%s'%protein],stdout=FNULL, stderr=subprocess.STDOUT)
        FNULL.close()
    """compute conservation of mutated position"""
    def conservation(self,newick_file,clustal_file,residue_start):
        #create blosum matrix and probability_matrix for conservation computation
        matrix = BlosumMatrix('./tools/blosum62.txt')
        probability_matrix = ProbabilityMatrix(clustal_file,matrix)

        t = Etree(newick_file)
        t.create_alignement_table(clustal_file)

        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)

        t.get_names()

        t.add_node_array()
        t.create_names_table(clustal_file)

        t.create_ID_table()
        rootWeightsArray = t.calculate_weights(probability_matrix)

        #get starting position of residue
        #start_pos = t.PDB_parse(self.pdb_id)

        #just for testing purpose
        """f = open(self.pdb_id+'_NEW.txt','r')
        out = open(+'_conservation_results1.txt','w')
        for line in f.readlines():
        	original_acid = line[0]
        	out.write(original_acid+ " ")
        	position = int(line[1:])
        	out.write(str(position)+ ' ')
            conservation_score = t.compute_conservation(clustal_file,start_pos,self.position,rootWeightsArray,self.wild_type)
         	out.write(str(conservation_score)+ '\n')
            """
        #compute conservation score
        conservation_score = t.compute_conservation(clustal_file,residue_start,self.position,rootWeightsArray,self.wild_type)
        print(conservation_score)
        return conservation_score

    """encode conservation score to int value"""
    def conservation_score_encoding(self,conservation_score):
        if(conservation_score >=0 and conservation_score <= 0.2):
            return 0
        elif(conservation_score > 0.2 and conservation_score <= 0.4):
            return 1
        elif(conservation_score > 0.4 and conservation_score <= 0.6):
            return 2
        elif(conservation_score > 0.6 and conservation_score <= 0.8):
            return 3
        else:
            return 4

    """compute index of mutated position in aligned sequence"""
    def compute_conservation(self,file,residue_start,index):
        count_mezera = 0
        count_basic_acid = 0
        count_mutated_acid = 0
        count_else=0
        all_count =0
        count_pos=0
        start_position = 1 #meni sa
        pos = 0
        handle = open(file,"r")

        lines = iter(handle.readlines())
        for line in lines:
            if(line.startswith('>')):
                continue
            else:
                for word in line.split():
                    #if(word[0] == '-'):
	 			    #break
		 			#if(word[0] == 'M'):
		 			#		count_pos -=1#-residue_start+1
                    if(residue_start > len(word)):
                        print(residue_start)
                        print(index)
                        count_pos = residue_start
                        print(count_pos)
                        for i in range(0,len(word),1):
                            if(word[i] != '-'):
                                count_pos +=1
                                if(count_pos == residue_start+index):
                                    pos = i
                                    print(word[i])
                                    break
                    else:
                        print(residue_start)
                        print(index)
                        count_pos = residue_start
                        if(residue_start < 0):
                            chain_res = index#+residue_start #+ abs(residue_start) + abs(residue_start) -1
                        elif (residue_start == 1):
                            chain_res= index+residue_start
                        else:
                            chain_res= index+residue_start+2

                        for i in range(0,len(word),1):
                            if(word[i] != '-'):
                                count_pos +=1
                                if(count_pos == chain_res):
                                    pos = i
		 							#print(pos)
                                    print("ACID:" + word[i])
                                    break
            break
        print(pos)
        return pos

    """create file sequences only"""
    def create_sequence_file(self,correlation_file,clustal_file):
        with open(clustal_file,'r') as f:
            lines = iter(f.readlines())
            with open(correlation_file,'w') as corr:
                for line in lines:
                    if(line.startswith('>')):
                        name = line.strip('>').strip('\n')
                        sequence = lines.next()
                        corr.write(sequence)
            corr.close()
        f.close()

    def correlation_scores(self,correlation_file,position,column_index,residue_start):
        #residue_start = self.PDB_parse(self.pdb_id)
        arr = []
        resArr = dict()
        nucleoColumn1 = dict()
        nucleoColumn2 = dict()

        with open(correlation_file,'r') as f:
            lines = list(iter(f.readlines()))
            correlation_array = []
            #iterate through all colums in sequence
            new_arr = dict()
            #compute correlaton scores
            for index in range(len(lines[0])-1):
                for line in lines:
                    #create arrays for mutated columns and other columns
                    nucleoColumn1[line[index]] = 0
                    nucleoColumn2[line[column_index]] = 0
                    #array for storing pairs from 2 columns
                    resArr[line[index]+line[column_index]] =0

                for line in lines:
                    #count number of acids/pairs in columns
                    nucleoColumn1[line[index]] += 1
                    nucleoColumn2[line[column_index]] +=1
                    resArr[line[index]+line[column_index]] +=1
                #count occurance of acid/pair in columns(in )
                for i in nucleoColumn1:
                    nucleoColumn1[i] /= 250
                for i in nucleoColumn2:
                    nucleoColumn2[i] /= 250
                for i in resArr:
                    resArr[i] /= 250

                correlation_score = 0
                for i in nucleoColumn1:
                    for j in nucleoColumn2:
                        key = i+j
                        if(key in resArr.keys()):
                            value = resArr[key] / (nucleoColumn1[i] * nucleoColumn2[j])
                            correlation_score += resArr[key] * math.log(value,2)
                correlation_array.append(correlation_score)
                #fl.write(str(correlation_score)+'\n')
                #print(correlation_array)
                #new_arr = resArr
                correlation_score = 0
                resArr = dict()
                nucleoColumn1 = dict()
                nucleoColumn2 = dict()
        f.close()
        return correlation_array

    #compute correlation of mutated position
    def correlation(self,position,column_index,residue_start):
        corr_scores = self.correlation_scores(self.pdb_id+"correlation.txt",position,column_index,residue_start)
        threshold = 1
        is_correlated = 0
        above_threshold = 0
        for item in corr_scores:
            if(item >= threshold):
                above_threshold +=1
            if(above_threshold > 0):
                is_correlated = 1
            else:
                is_correlated = 0
        return is_correlated

    """computes size change of original amino acid and mutated amino acid
     0 - acids are in the same size group
     + change from small to big one
     - change from big to small one
    """
    def size_change(self,x,y):
    	if(x in ['G','A','S','P']):
    		if(y in ['V','T','C','I','L','N','D','Q','K','E','M','H']):
    			return '+'
    		elif(y in ['G','A','S','P']):
    			return '0'
    		elif(y in ['F','R','Y']):
    			return '+'
    		elif(y == 'W'):
    			return '+'
    	elif(x in ['V','T','C','I','L','N','D','Q','K','E','M','H']):
    		if(y in['G','A','S','P']):
    			return '-'
    		elif(y in ['F','R','Y','W']):
    			return '+'
    		elif(y in ['V','T','C','I','L','N','D','Q','K','E','M','H']):
    			return '0'
    	elif(x in ['F','R','Y']):
    		if(y in ['G','A','S','P']):
    			return '-'
    		elif(y in ['V','T','C','I','L','N','D','Q','K','E','M','H']):
    			return '-'
    		elif(y in ['F','R','Y']):
    			return '0'
    		elif(y == 'W'):
    			return '+'
    	elif(x == 'W'):
    		if(y in ['G','A','S','P']):
    			return '-'
    		elif(y in ['V','T','C','I','L','N','D','Q','K','E','M','H']):
    			return '-'
    		elif(y in ['F','R','Y']):
    			return '-'
    		elif(y == 'W'):
    			return '0'


    """compute polarity change for wild-type amino acid and mutated acid"""
    def compute_polarity(self,x,y):
    	if self.__polarity[x] == '-':
    		if self.__polarity[y] == '+':
    			return '+'
    		elif self.__polarity[y] == '-':
    			return '0'
    	elif self.__polarity[x] == '+':
    		if self.__polarity[y] == '+':
    			return '0'
    		elif self.__polarity[y] == '-':
    			return '-'

    """computing difference in charge of wild-type acid and mutated acid"""
    def compute_charge(self,x,y):
    	if self.__charge[x] == '+':
    		if self.__charge[y] == '-':
    			return '-'
    		elif self.__charge[y] == '0':
    			return '-'
    		elif self.__charge[y] == '+':
    			return '0'
    	elif self.__charge[x] == '0':
    		if self.__charge[y] == '-':
    			return '-'
    		elif self.__charge[y] == '0':
    			return '0'
    		elif self.__charge[y] == '+':
    			return '+'
    	elif self.__charge[x] == '-':
    		if self.__charge[y] == '-':
    			return '0'
    		elif self.__charge[y] == '0':
    			return '+'
    		elif self.__charge[y] == '+':
    			return '+'

    """compute hydrophobicity index from difference of wild-type and mutant"""
    def compute_hydro_index(self,x,y):
    	if self.__hydro_index[x] > self.__hydro_index[y]:
    		return self.__hydro_index[y] - self.__hydro_index[x]
    	elif self.__hydro_index[x] < self.__hydro_index[y]:
    		return self.__hydro_index[y] - self.__hydro_index[x]
    	elif self.__hydro_index[x] == self.__hydro_index[y]:
    		return 0

    """encode secondary structure to one of 3 options
    H,G,I   H(helix)
    E,B     S(sheet)
    T,S,-   C(coil)
    """
    def sec_struc_code(self,id):
    	if(id == 'H' or id == 'G' or id == 'I'):
    		sec_struc_res = 'H'
    	elif(id == 'E' or id == 'B'):
    		sec_struc_res = 'S'
    	elif(id == 'T' or id == 'S' or id == '-'):
    		sec_struc_res = 'C'
    	return sec_struc_res

    """compute ASA value and secondary structure"""
    def compute_ASA(self):
    	p = PDBParser()
    	structure = p.get_structure("1a2p","1a2p.pdb")
    	model = structure[0]
        #treba nakoniec urobit strukturu s potrebnymi skriptami a tu dat cestu k dssp
    	dssp = DSSP(model,"1a2p.pdb",dssp='/home/juraj/dssp-2.0.4-linux-amd64')
    	#index is position of mutation
    	#zatial je asi potrebne mat residue_start+1 aby sedela hodnota ASA
        residue_start = self.PDB_parse(self.pdb_id)
        index = self.position
    	asa_key = list(dssp.keys())[index -residue_start+1-1]
    	asa = str(dssp[asa_key][3])
    	sec_structure = str(dssp[asa_key][2])
    	sec_struc_id = self.sec_struc_code(sec_structure)
    	#print("Skratka pre sek strukturu:%s" %sec_struc_id)
    	#print("Secondary structure %s" %dssp[a_key][2])
    	#print("ASA %.8f" %dssp[a_key][3])
    	return sec_struc_id,asa

    #divide asa values into 3 groups for predictor
    def divide_asa(self,asa):
    	encoded_asa = float(asa)
    	if(encoded_asa < 0.25 ):
    		return "burried"
    	elif (encoded_asa >= 0.25 and encoded_asa <= 0.5):
    		return "partially_burried"
    	elif (encoded_asa > 0.5):
    		return "exposed"

    #prepare mutation record
    def create_record(self):
        x = self.wild_type
        y = self.mutant

        self.get_fasta()
        self.get_pdb()
        self.create_newick_file()
        residue_start = self.PDB_parse(self.pdb_id)
        self.create_sequence_file(self.pdb_id+"correlation.txt",self.pdb_id+"clustal.fasta")
        mutated_column_index = self.compute_conservation(self.pdb_id+"correlation.txt",residue_start,self.position)

        #compute correlation of mutated position
        correlation = self.correlation(self.position,mutated_column_index,residue_start)
        #compute conservation of mutated position
        conservation = self.conservation(self.pdb_id+"_newick",self.pdb_id+"clustal.fasta",residue_start)
        conservation_code = self.conservation_score_encoding(conservation)
        #compute and encode polarity
        f_polarity = self.compute_polarity(x,y)
        f_polarity = self.__encode_polarity[f_polarity]
        #compute and encode charge
        f_charge = self.compute_charge(x,y)
        f_charge = self.__encode_charge[f_charge]
        #compute and encode hydrophobicity index
        f_hydro_index = str(self.compute_hydro_index(x,y))
        size = self.size_change(x,y)
        size = self.__encode_size[size]
        #compute asa and seconda structure
        struc_id,asa = self.compute_ASA()
        asa_val = self.divide_asa(asa)
        asa_val = self.__encode_asa[asa_val]
        struc_id = self.__encode_sec_struc[struc_id]
        #save inforrmations as a record for predictor
        #record is stored in file pred_record.txt
        record = open("mutation_record.csv","w")
        record.write('%s,%s,%s,%s,%s,%s,%s,%s\n' %("correlation","conservation","polaritychange","chargechange","hydroindexchange","secondarystruc","asa","sizechange"))
        record.write('%s,%s,%s,%s,%s,%s,%s,%s' %(correlation,conservation_code,f_polarity,f_charge,f_hydro_index,struc_id,asa_val,size))
        record.close()
