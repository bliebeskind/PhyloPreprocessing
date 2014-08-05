#! /usr/bin/env python

## For getting longest open reading frames from a file of nucleotide sequences
## *** Only forward frames are evaluated ***
##
## Function frames() is adapted from:
## 		http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec352


from Bio import SeqIO
import sys

def frames(rec, min_length=100):
	'''Given a Bio.Seq object, find forward reading frames and return Dictionary
	of structure: {orf_length: (start, end)}. Default minimum protein
	length if 100 AAs. Removes stop codons.'''
	reading_frames = {}
	seq_len = len(rec.seq)
	for frame in range(3):
		trans = str(rec.seq[frame:].translate())
		trans_len = len(trans)
		aa_start, aa_end = 0, 0
		while aa_start < trans_len:
			aa_end = trans.find("*", aa_start)
			if aa_end == -1: # if "*" not found
				aa_end = trans_len
			if aa_end - aa_start >= min_length:
				start = frame+aa_start*3
				end = min(seq_len, frame+aa_end*3+3) - 3 #remove stop codon
				reading_frames[(end - start)] = (start, end)
			aa_start = aa_end + 1 #start again after "*"
	return reading_frames
	
def get_orf(rec, min_length=100):
	'''Calls frames() and returns the nucleotide sequence of the longest orf.'''
	D = frames(rec, min_length)
	top_orf = sorted(D.keys())[-1]
	sequence = rec.seq[D[top_orf][0]:D[top_orf][1]]
	return sequence
	
def get_orfs_from_file(infile, min_length=100):
	'''Open fasta infile and return iterator of SeqRecords with ORF sequences.'''
	records = SeqIO.parse(infile, 'fasta')
	for rec in records:
		try:
			rec.seq = get_orf(rec, min_length)
		except IndexError: # frames() found nothing above min length
			continue
		yield rec
		
def get_prots_from_file(infile, min_length=100):
	'''Open fasta infile and return iterator of SeqRecords with protein sequences.'''
	records = SeqIO.parse(infile, 'fasta')
	for rec in records:
		try:
			rec.seq = get_orf(rec, min_length).translate()
		except IndexError: # frames() found nothing above min length
			continue
		yield rec

	
if __name__ == '__main__':
	infile = sys.argv[1]
	outfile = sys.argv[2]
	prots = bool(sys.argv[3])
	try:
		min_length = int(sys.argv[4])
		if prots:
			outrecs_iter = get_prots_from_file(infile, min_length)
		else:
			outrecs_iter = get_orfs_from_file(infile, min_length)
	except IndexError: # no min length specified
		if prots:
			outrecs_iter = get_orfs_from_file(infile, min_length=100)
		else:
			outrecs_iter = get_orfs_from_file(infile, min_length=100)
	with open(outfile,'w') as f:
		SeqIO.write(outrecs_iter, f, 'fasta')
