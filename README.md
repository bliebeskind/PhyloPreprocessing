PhyloPreprocessing
========

Collection of little scripts for wrangling files before a phylogenetic analysis.

Script | Useful for you? | Description |
------ | ----------------- | ----------- |
aln_to_subtree | yes | remove taxa from a tree to match those in an alignment.
domain_chop | no | divide 4 domain protein into constituent domains. Used in [Liebeskind et al. 2013](http://www.cell.com/current-biology/abstract/S0960-9822%2813%2901141-X)
ensembl_tsv_ids | no | for parsing tab delimited files downloaded from Ensembl after an ortholog search
entrez_xml_mods | no | functions for dealing with entrez xml. Most are broken
fas_to_nex | yes | convert fasta file to nexus
fas_to_phy | yes | convert fasta file to phylip (good for PAML)
hmmalign_parser | no | parse a stockholm alignment output by hmmer
keepers_subAlnTree | yes | get a subalignment and tree with an text file of taxa you wish to retain
longest_cds | yes | for getting longest transcripts from files of Ensembl genes with multiple transcripts
open_reading_fram | yes | get longest open reading frames for each sequence in a fasta file (good for NGS data but only looks at forward frame, so it works for RNAseq)
phylip_map | maybe | replace long description lines with short names, and keep track of changes
remove_redundant | yes | remove redundant sequences from a sequence file
stockholm_pps | no | explore posterior probabilities on stockholm alignment such as those put out by hmmer
TPC_chop | no | divide 2 domain protein into constituent domains. Used in [Liebeskind et al. 2013](http://www.cell.com/current-biology/abstract/S0960-9822%2813%2901141-X) 