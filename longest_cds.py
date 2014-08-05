#! /usr/bin/env python

###	For handling fasta files of Ensembl genes with multiple transcripts.  Outputs
### the longests transcripts associated with each gene and a file called Transcript
### Counts which lists the number of transcripts associated with each gene.


from Bio import SeqIO
import sys

def geneid_dict(infile):
	'''Opens a fasta file of Ensembl CDS sequences, some of which are transcripts
	of the same gene.  Returns a dictionary with gene ids as keys and SeqRecords
	of the CDSs as values.'''
	D = {}
	with open(infile) as f:
		records = SeqIO.parse(f, "fasta")
		for rec in records:
			gene_id = rec.id.split('|')[0]
			if gene_id in D:
				D[gene_id].append(rec)
			else:
				D[gene_id] = [rec]
	return D
	
def longest_sequence(lst):
	'''Takes list of SeqRecords objects and returns the SeqRecord with the longest
	sequence, not counting "N" characters.'''
	longestseq = ''
	for i in lst:
		seq = str(i.seq).replace("N", "")
		try:
			assert len(seq) != len(longestseq) # Is it a tie?
		except AssertionError:
			sys.stderr.write("Tie between %s and %s\n" % (longestrec.id, i.id))
		if len(seq) > len(longestseq):
			longestrec = i
			longestseq = str(i.seq)
	return longestrec
	
def longest_list(infile):
	'''Calls geneid_dict.  Loops over entries, writes out trancript counts for 
	each gene id, and returns a list of the longest transcript for each gene id.'''
	D = geneid_dict(infile)
	longest_cdss = []
	transcript_counts = ''
	for k, v in D.iteritems():
		#Can be validated with grep ">" file | cut -d "|" -f 1 | sort | uniq -c
		transcript_counts += "Gene: %s, Transcript Count: %i\n" % (k, len(v))
		longest_cdss.append(longest_sequence(v))
	return longest_cdss, transcript_counts
	
	
if __name__ == '__main__':
	Usage = "## Usage: longest_cds.py <infile> <outfile> ##"
	try:
		infile, outfile = sys.argv[1], sys.argv[2]
		CDSs, transcript_counts = longest_list(infile)
		SeqIO.write(CDSs, outfile, "fasta")
		with open("Transcript_Counts", "a") as f:
			f.write(transcript_counts)
	except IndexError:
		sys.stderr.write(Usage)
