#! /usr/bin/env python

from Bio import AlignIO
import sys

input_alignment = sys.argv[1]

def get_pps(infile):
	D = {}
	aln = AlignIO.read(infile,'stockholm')
	for rec in aln:
		D[rec.id] = rec.letter_annotations['posterior_probability']
	return D
	
def print_data(infile):
	pp_bins = ['.','1','2','3','4','5','6','7','8','9','*']
	D = get_pps(infile)
	print "Ins,.1,.2,.3,.4,.5,.6,.7,.8,.9,1"
	for rec in D.iterkeys():
		ppstring = ''
		for b in pp_bins:
			ppstring += str(D[rec].count(b)) + ','
		print rec + ',' + ppstring.rstrip(',')

if __name__ == '__main__':
	print_data(input_alignment)
