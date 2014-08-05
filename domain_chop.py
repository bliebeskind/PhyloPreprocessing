#! /usr/bin/env python

## For cutting four-domain proteins into constituent domains
## Used in Liebeskind et al. 2013, Current Biology

from Bio import AlignIO
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def create_domain(seq, name):
        '''Return SeqRecord object from sliced domain, removing gap characters.'''                                             
        seq_record = SeqRecord(Seq(seq.replace('-', ''), generic_protein), id=name,\
	description = '')
        return seq_record

def get_domains(alignment,cut_list):
	'''
	Input is a BioPython Alignment object and a list of cutting sites that
	divide the domains of the constituent sequences. Returns a list of the
	domains divided by the cut sites as ungapped SeqRecord objects.
	'''
	domains = []                                 
	for record in alignment:
		first_cut = True
		for cut in range(len(cut_list)):
			assert cut_list[cut] <= len(record), "Cut out of alignment range"
			assert cut_list[cut] is int, "Cuts should be integers"
			if first_cut:
				domain = str(record.seq[:cut_list[cut]])
				rec = create_domain(domain, record.id + "_D%i" % (cut+1))
				domains.append(rec)
				first_cut = False
			else:
				domain = str(record.seq[cut_list[cut-1]:cut_list[cut]])
				rec = create_domain(domain, record.id + "_D%i" % (cut+1))
				domains.append(rec)
	return domains



if __name__ == '__main__':
	
	Usage = """
	## Usage ##
  domain_chop.py <infile> <outfile> <cut1> <cut2> <cut3>                                                       
  Outputs contituent domains with gap characters removed.
  Alignment must be in fasta format
"""
	
	import sys
	from Bio import SeqIO
	
	try:
		alignment = AlignIO.read(sys.argv[1],'fasta')
		outfile = sys.argv[2]
		cut_list = [int(cut) for cut in sys.argv[3:]]
		domains = get_domains(alignment,cut_list)
		SeqIO.write(domains,outfile,'fasta')
	except:
		print Usage
