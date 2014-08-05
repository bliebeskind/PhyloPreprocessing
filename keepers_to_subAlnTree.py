#! /usr/bin/env python

### For getting subalignment and tree with an text file of taxa you wish to 
### retain.  Text file must be in format of a Python list
### Default input/output of alignment is fasta/fasta
### Default input/output of tree is nexus/newick - suppresses all internal labels

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
import dendropy, sys

def get_keepers(infile):
	'''Infile should be text of a python list of sequences to keep'''
	with open(infile) as f:
		return eval(f.read())

def aln_keepers(alignment, keepers, format='fasta'):
	'''Given list of keepers and an alignment return alignment that's been 
	stripped of all taxa except keepers'''
	records = SeqIO.parse(alignment, format)
	outrecs = []
	outnames = []
	for rec in records:
		species = rec.id[:rec.id.find("SCN")]
		assert len(species) != len(rec.id), "Didn't find 'SCN' in %s" % rec.id
		if species in keepers:
			outrecs.append(rec)
			outnames.append(rec.id)
	return outrecs, outnames
		
def trim_tree(input_tree, keepers, format='nexus'):
	'''Trim dendropy tree to tips in keepers'''
	tree = dendropy.Tree.get_from_path(input_tree,format)
	outnodes = []
	taxa = tree.infer_taxa()
	for t in taxa:
		taxon = str(t).replace(' ', '_')
		species = taxon[:taxon.find("SCN")]
		if species in keepers:
			outnodes.append(str(t))
	tree.retain_taxa_with_labels(outnodes)
	return tree, outnodes
	
if __name__ == '__main__':
	Usage = '''
############# Usage ###############
keepers_to_subAlnTree.py <keeper file> <alignment> <tree>
keeper file must be in format of python list
	'''
	try:
		keepers = sys.argv[1]
		alignment = sys.argv[2]
		input_tree = sys.argv[3]
	except IndexError:
		print Usage
		
	keeper_list = get_keepers(keepers)
	outrecs, outnames = aln_keepers(alignment, keeper_list)
	tree, outnodes = trim_tree(input_tree, keeper_list)
	assert set(outnames) == set(outnodes), "Alignment and Tree sets don't match"
	tree.write_to_path("trimmed.nhx",'newick',suppress_edge_lengths=True,\
suppress_internal_node_labels=True)
	AlignIO.write(MultipleSeqAlignment(outrecs),'trimmed.fas','fasta')
