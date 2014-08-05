#! /usr/bin/env python

## For dividing two pore channels into their constituent domains.  Takes an 
# alignment and divides each seq in two.  Outputs raw sequences with gap 
# characters removed.

Usage = """
Usage:
   domain_chop.py <infile> <outfile> <cut1 (location in alignment)> 

Outputs contituent domains with gap characters removed.         
"""

import sys
from Bio import AlignIO, SeqIO
from Bio.Alphabet import generic_protein                             
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment                    

records = []
def create_domain(seq, name):
        """Makes seq record objects, taking seq, name and description
        arguments.  Appends to pre-defined list"""                                             
        seq_record = SeqRecord(Seq(seq.replace('-', ''), generic_protein), id=name,\
	description = '')
        records.append(seq_record)

try:
        infile = sys.argv[1]
        outfile = sys.argv[2]                                        
        cut1 = int(sys.argv[3])
        alignment = AlignIO.read(infile, 'fasta')
        for record in alignment:
                domain1 = str(record.seq[:cut1]) #makes cuts
                domain2 = str(record.seq[cut1:])
                create_domain(domain1, record.id + '_D1')       
                create_domain(domain2, record.id + '_D2')
        SeqIO.write(records, outfile, "fasta")                           
         
except IndexError:                                                   
        print 'Takes 3 variables: infile, outfile, cut1'
        print Usage                                               
                                                            
                 
 
