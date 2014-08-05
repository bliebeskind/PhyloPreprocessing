#! /usr/bin/env python

from Bio import SeqIO
import sys

def remove_repeats(infile,format='fasta'):
	'''Return generator over non-redundant sequences. Only one sequence of
	a pair whose sequences match will be returned.'''
	nonred_seqs = []
	with open(infile) as f:
		records = SeqIO.parse(infile,format)
		for rec in records:
			if str(rec.seq) not in nonred_seqs:
				yield rec
				nonred_seqs.append(str(rec.seq))
			else:
				print rec.id    

def print_nonred_seqs(infile,outfile,informat='fasta',outformat='fasta'):
	'''Write non-redundant sequences from infile to outfile.'''
	SeqIO.write(remove_repeats(infile,informat),outfile,outformat)
    	
if __name__ == '__main__':
	Usage = """
###     Usage     ###
aln_cleaner.py <fasta_infile> <outfile> >> removed_seqs.txt
	"""
	try:
		infile = sys.argv[1]
		outfile = sys.argv[2]
		print_nonred_seqs(infile,outfile)
	except IndexError:
		if sys.argv[1] == 'h':
			print Usage
		else:
			sys.stderr.write('Takes two arguments:\n' + Usage)
