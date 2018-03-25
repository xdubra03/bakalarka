#!/usr/bin/python

import sys

#compute index of mutated position in aligned sequence
def compute_conservation(file,residue_start,index,wild_type):
    count_pos = 0
    pos = 0
    handle = open(file,"r")

    lines = iter(handle.readlines())
    for line in lines:
        if(line.startswith('>')):
            continue
        else:
            for word in line.split():
                if(residue_start > len(word)):
                    count_pos = 0
                    print(count_pos)
                    for i in range(0,len(word),1):
                        if(word[i] != '-'):
                            count_pos +=1
                        if(count_pos == residue_start+index-residue_start):
                            pos = i
                            print(word[i])
                            break
                else:
                    print(residue_start)
                    print(index)
                    count_pos = 0
                    if(residue_start < 0):
                        chain_res = index + abs(residue_start)#+residue_start #+ abs(residue_start) + abs(residue_start) -1
                    elif (residue_start == 1):
                        chain_res= index+residue_start-1
                    else:
                        chain_res= index+residue_start-1

                    for i in range(0,len(word),1):
                        if(word[i] != '-'):
                            print(word[i])
                            count_pos +=1
                            if(count_pos == chain_res):
                                pos = i
                                print(pos)
                                print("ACID:" + word[i])
                                if(word[i] == wild_type):
                                    print('Prvy pokus OK')
                                    break
                                else:
                                    count_pos = 0
                                    chain_res = index + residue_start - residue_start
                                    for i in range(0, len(word), 1):
                                        if (word[i] != '-'):
                                            print(word[i])
                                            count_pos += 1
                                            if (count_pos == chain_res):
                                                pos = i
                                                print(pos)
                                                print("ACID:" + word[i])
                                                if (word[i] == wild_type):
                                                    print('Druhy pokus OK')
                                                    return pos
                                                else:
                                                    print('Neuspech')

            break
        print(pos)
        return pos

compute_conservation(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4])