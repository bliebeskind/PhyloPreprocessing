#! /usr/bin/env python

from Bio import SearchIO
import numpy as np
import sys

np.seterr(all='raise') # set numpy warnings to actually be warnings

def parse_hmmscan_tab(infile, print_header=True):
    '''Parse hmmscan output in --tblout format'''
    if print_header:
        yield "query","top hit","evalue","certainty","num sig hits"
    records = SearchIO.parse(infile,'hmmer3-tab')
    for rec in records:
        query = rec.id
        if len(rec) > 1:
            hit1,hit2 = rec.hits[0],rec.hits[1]
            eval1,eval2 = hit1.evalue,hit2.evalue
            if eval1 != 0: # convert to -log10 evalue
                eval1 = -np.log10(eval1)
            if eval2 != 0:
                eval2 = -np.log10(eval2)
            if eval1 == 0 and eval2 != 0: # this may be a hack, I don't care
                certainty = 1
            elif eval1 == 0 and eval2 == 0:
                certainty = 0
            else: # calculate certainty with info theoretic calc.
                total = eval1 + eval2
                p1,p2 = eval1/total, eval2/total
                certainty = 1 + (p1 * np.log2(p1)) + (p2 * np.log2(p2))
        else:
            certainty = 1
        yield query, rec.hits[0].id, rec.hits[0].evalue, certainty, len(rec)

def write_hmmscan_tab(infile,outfile,sep='\t',print_header=True):
    '''Write a file with four fields:
        Query name
        Top hit
        E-value
        Certainty (1 + p1*log2(p1) + p2*log2(p2), where p1 and p2 are -log10 pvalues of first 2 hits)
        Number of hits'''
    count = 0
    with open(outfile,'w') as out:
        for line in parse_hmmscan_tab(infile,print_header):
            out.write(sep.join(map(str,line)) + '\n')
            count += 1
            if count % 10 == 0:
                print count

if __name__ == '__main__':
    try:
        infile,outfile = sys.argv[1], sys.argv[2]
    except IndexError:
        raise Exception("Usage: hmmscan_results.py <infile> <outfile>")
    write_hmmscan_tab(infile,outfile)
