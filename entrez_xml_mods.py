#! /usr/bin/env python

def entrez_read(infile):
	from Bio import Entrez
	with open(infile) as f:
		record = Entrez.read(f)
	return record
	
def entrez_parse(infile):
	from Bio import Entrez
	with open(infile) as f:
		record = Entrez.parse(f)
	return record	

def get_taxonomy(record):
	rec = entrez_read(record)
	return rec[0]['GBSeq_taxonomy']

def get_accession_version(record):
	rec = entrez_read(record)
	return rec[0]['GBSeq_accession-version']	

def xml_dict_fields(field):
	return field.keys()

def xml_list_fields(field):
	return [item for item in field]

def get_all_taxonomies():
	"""returns dictionary with accession = key, and taxonomy = value
	 for all genbank xml files in directory"""
	import glob
	D = {}
	file_list = glob.glob('*.xml')
	for file in file_list:
		record = entrez_read(file)
		key = record[0]['GBSeq_accession-version']
		val = record[0]['GBSeq_taxonomy']
		D[key] = val
	return D
		
def get_accession_and_taxonomy(accession):
	"""Takes accession number of entrez xml file as input.
	Opens saved file in working directory.
	Returns list accession=[0] and taxonomy = [1]
	 """
	import glob
	key_value = []
	f = glob.glob(accession + '.xml')
	for filename in f:
		record = entrez_read(filename)
		key = record[0]['GBSeq_accession-version']
		val = record[0]['GBSeq_taxonomy']
		key_value.append(key)
		key_value.append(value)
	return key_value
			
def get_taxon_taxonomy_accession(accession):
	"""Takes accession number of entrez xml file as input.
	Opens saved file in working directory.
	Returns list [taxon, taxonomy, accession]
	 """
	import glob
	org_tax_acc = []
	f = glob.glob(accession + '.xml')
	for filename in f:
		record = entrez_read(filename)
		org_tax_acc.append(record[0]['GBSeq_organism'])
		org_tax_acc.append(record[0]['GBSeq_taxonomy'])
		org_tax_acc.append(record[0]['GBSeq_accession-version'])
	return org_tax_acc

def count_otu(dictionary, otu):
        """returns number of a certain taxonomic unit present in 
        dictionary like that output by get_all_taxonomies"""
        count = 0
        for val in dictionary.itervalues():
                if otu in val:
                        count +=1
        return count          
                                
                        
                
