#! /usr/bin/env python
                             
def get_accessions(infile):
        """Open fasta file of sequences with gi lines and
        return list of accession numbers"""    
        with open(infile) as f:
                return [line.split('|')[3] for line in f if line.startswith('>gi')] 
                                     
def make_info_list(accessions):
        """Takes list of accession numbers
        Calls entrez_xml_mods
        Returns list of 3 part lists = [taxon, taxonomy, accession]
        """                                 
        return [entrez_xml_mods.get_taxon_taxonomy_accession(a) for a in accessions] 
        
def check_non_gbs(infile):
        #import sys  
        with open(infile) as f:
                for line in f: 
                        if line.startswith('>') and line.startswith('>gi') == False:
                                sys.stderr.write("%s not in GB format \n" % line.strip('>')) 
                                         
# **Driver** 
# Print information in tab-delimited file  
if __name__ == '__main__':
	try:
	        import entrez_xml_mods
		import sys
		infile = sys.argv[1]
		for info in make_info_list(get_accessions(infile)):
		        #prints only the 4th value of taxonomy          
		        print info[0], '\t', info[1].split(';')[3], '\t', info[2]  
                check_non_gbs(infile)
	except IndexError:
		print Usage
                                                   

