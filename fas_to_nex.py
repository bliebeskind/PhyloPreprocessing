#! /usr/bin/env python

from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

def fas_to_nex(infile,outfile,protein=True):
	'''Convert fasta infile to nexus and write to outfile. Uses BioPython'''
	if protein:
		aln = AlignIO.read(infile,'fasta',alphabet=Gapped(IUPAC.extended_protein))
	else:
		aln = AlignIO.read(infile,'fasta',alphabet=Gapped(IUPAC.unambiguous_dna))
	AlignIO.write(aln,outfile,'nexus')
