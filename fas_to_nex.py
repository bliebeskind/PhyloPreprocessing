#! /usr/bin/env python

import sys
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

def fas_to_nex(infile,outfile,protein=True):
	'''Convert fasta infile to nexus and write to outfile. Uses BioPython'''
	if protein:
		aln = AlignIO.read(infile,'fasta',alphabet=Gapped(IUPAC.extended_protein))
	else:
		aln = AlignIO.read(infile,'fasta',alphabet=Gapped(IUPAC.unambiguous_dna))
	AlignIO.write(aln,outfile,'nexus')
	
if __name__ == '__main__':
	Usage = """
	#### Usage ####
	fas_to_nex.py <infile> <outfile>
	"""
	try:
		infile,outfile = sys.argv[1], sys.argv[2]
		fas_to_nex(infile,outfile)
	except IndexError:
		print Usage
