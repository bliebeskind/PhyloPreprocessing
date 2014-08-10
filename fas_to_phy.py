#! /usr/bin/env python

import sys
from Bio import SeqIO

def fas_to_phy(infile):
	'''Convert a fasta alignment into phylip format'''
	records = list(SeqIO.parse(infile,'fasta'))
	first_line = True
	for rec in records:
		if first_line:
			print " " + str(len(records)) + " " + str(len(rec.seq))
			first_line = False
		print rec.id + ' ' + str(rec.seq)
		
if __name__ == '__main__':
	Usage = """
	####### Usage #######
	fas_to_phy.py <myFasta> >> my_phylip.phy
	"""
	try:
		fas_to_phy(sys.argv[1])
	except IndexError:
		print Usage
