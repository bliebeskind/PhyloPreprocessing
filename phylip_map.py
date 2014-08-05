#! /usr/bin/env python

## Adapted from BioPython Cookbook 6.2.1
# Takes a newick tree as input and replaces names with seq0, seq1...

Usage = """
phylip_map.py <fasta_alignment> 

Prints name mapping to stdout; pipe to file to save mapping
Saves a file called "your_fasta_file" + .phy
"""

from Bio import AlignIO
from pprint import pprint as pp
import sys

infile = sys.argv[1]
infilecore = infile.split('.')[0]

alignment = AlignIO.read(infile, "fasta")
name_mapping = {}
for i, record in enumerate(alignment):
    i = "seq%i" % i #convert enumerated digit to string 'seq0','seq1'...
    name_mapping[i] = record.id
    record.id = i
pp(name_mapping)

AlignIO.write([alignment], infilecore + ".phy", "phylip")

