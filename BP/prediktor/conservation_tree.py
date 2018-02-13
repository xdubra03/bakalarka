from __future__ import division
import pandas as pd
import Bio
from Bio.PDB import *
import urllib2
import os
import shutil
import sys
import subprocess
from ete3 import Tree
import copy
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

    #get fasta file for entered pdb id and chain
    def get_fasta(self):
    	fasta_file = self.name + ".fasta"
    	fasta_output = open(fasta_file,"w")
        url_f = "https://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=fastachain&compression=NO&structureId=%s&chainId=%s"%(self.name,self.chain)
    	#url_f = "http://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=%s&compressionType=uncompressed"%self.name
    	try:
    		h = urllib2.urlopen(url_f)
    	except URLError as error:
    		print(error.reason)
    		sys.exit(1)
    	fasta_output.write(h.read())
    	fasta_output.close()
        i = 0
        fasta_cleaned = open(self.name+"FASTA.fasta","w")
        handle = open(fasta_file,"r")
        lines = iter(handle.readlines())

        for line in lines:
        	if(line.startswith('>')):
        		i+=1
        	if(i < 2):
        		fasta_cleaned.write(line)
        fasta_cleaned.close()

    def parse_XML(self,name):
    	i=0
    	output = open(name+".txt","w")
    	f = open(name+".xml","r")
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

    def run_blast(self,name):
    	subprocess.call(['./blastp','-query','%sFASTA.fasta'%name,'-db','nr','-outfmt','5','-out','%s.xml'%name,'-max_target_seqs','250','-remote'])

    def run_clustal(self,name):
    	in_file = name + ".txt"
    	out_file = name +"clustal.fasta"
    	clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
    	subprocess.call(['./clustalo','-i','%s'%in_file,'--outfmt=vie','-o','%s'%out_file,'--auto','-v','--force'])

    def create_newick_file(self,name):
        self.get_fasta()
        self.run_blast(name)
    	self.parse_XML(name)
    	self.run_clustal(name)
    	protein = name+'clustal.fasta'
    	output = name+'output'
    	subprocess.call(['./FastTree','-out','%s'%output,'%s'%protein])

    #get pdb file for entered pdb id
    def get_pdb(self,name):
		protein_file = self.name + ".pdb"
		pdb_output  = open(protein_file, "w")
		url_pdb = "https://files.rcsb.org/download/%s.pdb" %self.name
		try:
			handle = urllib2.urlopen(url_pdb)
		except URLError as error:
			print(error.reason)
			sys.exit(1)
		pdb_output.write(handle.read())
		pdb_output.close()

     #parse PDB file to find stating position
    def PDB_parse(self,name):
		p = PDBParser()
		structure = p.get_structure(self.name,self.name+".pdb")
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


    def calculate_weights(self):
		"""calculates the values in weights array in each node"""

		#fudge factor constant to prevent 0 in the weights array
		fugde_factor = 0.1

		#traverse the tree and compute values in each node
		for node in t.traverse('postorder'):

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



		return t.get_tree_root().weightsArray
