#!/usr/bin/python2
from __future__ import division
import sys
import Bio
from Bio import SeqIO


file = sys.argv[1]
original_acid = sys.argv[2]
mutated_acid = sys.argv[3]
position = sys.argv[4]
index = int(position)
count_mezera = 0
count_basic_acid = 0
count_mutated_acid = 0
count_else=0
all_count =0
count_pos=0
start_position = 1 #meni sa
pos = 0
#o = open(file+"output.txt","w")
#for seq_record in SeqIO.parse(file+"clustal.fasta","fasta"):
#	o.write(str(seq_record.seq)+'\n')
	#print(str(seq_record.seq))

# 0 index = zaciatocna pozicia z pdb 

#o.close()
handle = open(file+"output.txt","r")
for line in handle.readlines():
	
 	for word in line.split():
 		#if(word[0] == '-'):
 		#	break
 		for i in range(len(word)):
 			if(word[i] != '-'):
 				count_pos +=1
 				if(count_pos == index):
 					pos = i
 					print(pos)
 					print(word[i])
 					break
 	break				
handle.close()
handle = open(file+"output.txt","r")
for line in handle.readlines():
	all_count += 1
	for word in line.split():
		if word[pos] == "-":
			count_mezera += 1
		elif word[pos] == original_acid:
			count_basic_acid +=1
		elif word[pos] == mutated_acid:
			count_mutated_acid += 1
		else:
			count_else +=1				

handle.close()

print ("---------------\n")
print ("original acid: %s" %original_acid)
print ("mutated acid: %s" %mutated_acid)
print ("at position: %s\n" %position)
print ("Statistics:\n")
print ("No all proteins: %d" % all_count)
print ("Is same at position: count = %d" %count_mezera)
print ("No. basic acids at position: %d" %count_basic_acid)
print ("No. mutated acids at position: %d" %count_mutated_acid)
print ("No. else acids at position: %d" %count_else)
percentage = (float((count_basic_acid) / (all_count))) *100
percentage_mutated = (float(count_mutated_acid) / (all_count))*100
percentage_else = (float(count_else) / (all_count))*100
print ("Originals : %f" %percentage)
print ("Mutated : %f" %percentage_mutated)
print ("Else amino acid: %f"%percentage_else)