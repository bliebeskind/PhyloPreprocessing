#! /usr/bin/env python

## For parsing tab delimited files downloaded from Ensembl after an ortholog
## search. Script will take a list of files (different orthologous groups) and
## print out a file with headers for each species, and a list of collected 
## genes for that species.  These can then be used to download the sequences
## from Biomart.

import sys,re

def make_ids_D(src, D, infile):
	'''Skip header (first line) and split subsequent rows of a tab-delimited
	file. Extract Ensembl ids or print error message if id doesn't match 
	Ensembl format. Make a dictionary where species names are keys and a 
	list of gene ids are the values.'''
	line = src.readline() #skip header
	rowcount = 1
	for row in src:
		splitrow = row.split('\t')
		species = splitrow[0]
		geneid = re.search("ENS\w+\d{11}", splitrow[3])
		if geneid and species not in D: #first entry for species?
			D[species] = [geneid.group()]
		elif geneid:
			D[species].append(geneid.group())
		else:
			sys.stderr.write("no match row %i in %s\n" % (rowcount,infile))
		rowcount +=1
	return D

def concat_Ds(src):
	'''Loop over list of infiles and update dictionary for each by calling
	make_ids_D'''
	D = {}
	for infile in src:
		with open(infile) as f:
			D = make_ids_D(f, D, infile)
	return D

def extract_ids_from_D(src):
	'''Call concat_Ds to get dictionary. Loop over keys, which are species
	names, and print all the genes associated with that species under the 
	species name as a header.  Also print number of genes for each species.'''
	D = concat_Ds(src)
	for species in D.iterkeys():
		gene_list = D[species]
		print species + " %i gene(s)" % len(gene_list)
		for gene in gene_list:
			print gene
		print '' #print blank line

if __name__ == '__main__':
	file_list = sys.argv[1:]
	extract_ids_from_D(file_list)