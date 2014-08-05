#! /usr/bin/env python

# For constraining a tree to only have the taxa in an input alignment.
# Assumes the alignment and desired tree are sub-sets of the originals.

import dendropy, sys
from Bio import SeqIO

def get_keeper_tips(inalignment, format):
	'''Take input alignment and return list of taxa'''
	records = SeqIO.parse(inalignment, format)
	tips = []
	for rec in records:
		if rec.id not in tips:
			tips.append(rec.id)
	return tips
	
def trim_tree(intree, tips, outfile, informat='nexus', outformat='newick'):
	'''Takes a tree and a list of taxa and trims tree to only this list'''
	tree = dendropy.Tree.get_from_path(intree, informat)
	tree.retain_taxa_with_labels(tips)
	tree.write_to_path(outfile, outformat, suppress_internal_node_labels=True,\
suppress_edge_lengths=True)
	
	
if __name__ == '__main__':
	Usage = """
######################## Usage ############################
aln_to_subtree.py <input nexus tree> <input nexus alignment>
	"""
	try:
		intree = sys.argv[1]
		inalignment = sys.argv[2]
		core = inalignment.split('.')[0]
		tips = get_keeper_tips(inalignment, 'fasta')
		trim_tree(intree, tips, core+'.nhx')
	except IndexError:
		print Usage
