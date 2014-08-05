#! /usr/bin/env python

# For parsing hmmalign output. Makes a dictionary where keys are sequence names,
# and values are a list of the following four elements: list of aligned residues,
# list of posterior values for each sequence, list of posterior values for all
# columns, reference assignments for all columns.  Then takes only the columns 
# that have been assigned as consensus homologous columns (those marked with an
# 'x' in the RF rows), and outputs these in fasta format.

# Output is not in 'aligned' order.  Could fix this by collecting seqs as
# 'Seq' objects with Biopython in a list and writing out.

import random,sys

def add_sequence(line,D):
    '''Adds sequence to dictionary as first element in value list.
    key is sequence name'''
    tup = (line.split(' ')[0],line.split(' ')[-1])
    if tup[0] not in D: #is sequence name in dictionary already?
    	seq = [aa for aa in tup[1]]
    	D[tup[0]] = [seq] #add name and list of AAs
    else:
    	D[tup[0]][0] = D[tup[0]][0] + [aa for aa in tup[1]]
    return D

def add_GR(line,D):
    '''Adds posterior value for each site of a sequence to dictionary as second
    element of value list.'''
    parts = line.split(' ')
    if len(D[parts[1]]) < 2: #does entry have a GR dimension yet?                  
    	D[parts[1]] = D[parts[1]] + [[gr for gr in parts[-1]]]
    else:
    	D[parts[1]][1] = D[parts[1]][1] + [gr for gr in parts[-1]]    
    return D
    
def add_PP(line,D):
    '''Adds posterior value for each column to all sequences in dictionary.
    Third element in value list'''
    parts = line.split(' ')
    if len(D[random.choice(D.keys())]) < 3: #does entry have a PP dimension yet?  
    	for k in D.iterkeys():
    	    D[k] = D[k] + [[pp for pp in parts[-1]]]
    else:
    	for k in D.iterkeys():
    	    D[k][2] = D[k][2] + [pp for pp in parts[-1]]    
    return D
    
def add_RF(line,D):
    '''Adds reference assignments for each column to all sequences in dictionary.
    Fourth element in value list'''
    parts = line.split(' ')
    if len(D[random.choice(D.keys())]) < 4: #does entry have a RF dimension yet?  
    	for k in D.iterkeys():
    	    D[k] = D[k] + [[rf for rf in parts[-1]]]
    else:
    	for k in D.iterkeys():
    	    D[k][3] = D[k][3] + [rf for rf in parts[-1]]    
    return D

def parse_hmmalign(infile):
    '''Step through infile and make dictionary where keys are sequence names,
    and values are a list of four elements made by add_sequence, add_GR,
    add_PP, and add_RF.'''
    D = {}
    with open(infile) as f:
        for line in f:                                 
            line = line.strip('\n').strip('//').strip()
            if not line.startswith('#') and line !='':     
        	D = add_sequence(line,D)
            elif line.startswith('#=GR'):    
        	D = add_GR(line,D)
            elif line.startswith('#=GC PP_cons'):
        	D = add_PP(line,D)
            elif line.startswith('#=GC RF'):                             
        	D = add_RF(line,D)
            elif line.startswith('# STOCKHOLM'):         
            	pass
            else:
            	pass                                           
    return D 
    
def print_fasta(infile): 
    '''Prints only those columns marked by 'x' in RF row.  Output is in fasta.'''
    D = parse_hmmalign(infile)
    for k in D.iterkeys():
    	count = 0
    	seq = ''                                      
    	if D[k] != '':
    	    #sys.stderr.write(k+'\n')
    	    for item in D[k][3]:
    	    	if item == 'x':
    	    	    seq += D[k][0][count]
    	        count +=1                           
    	    print '>'+k                                     
	    print seq	    	                         

if __name__ == '__main__':
    try:
    	infile = sys.argv[1]                     
    	print_fasta(infile)  
    except IndexError:                 
    	sys.stderr.write("Takes one argument = infile\n")
