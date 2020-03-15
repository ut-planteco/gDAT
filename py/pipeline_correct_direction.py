#!/usr/bin/env python
from __future__ import division

import os
import argparse
import sys

"""
	Reverse complements FASTA sequences based on BLAST+ by searching for +/- strand hits.
"""
parser = argparse.ArgumentParser(description = """ 
	Reverse complements FASTA sequences based on BLAST+ by searching for +/- strand hits.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', required = True, type = open, help = """
	specify a FASTA file to be corrected
	""")
parser.add_argument(
	'-b', metavar = 'BLAST_FILE', required = True, type = open, help = """
	specify a BLAST+ file containing hit description with query strands
	""")

args = parser.parse_args()

def reverse_complement(seq):   
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	bases = list(seq) 
	bases = reversed([complement.get(base,base) for base in bases])
	return ''.join(bases)

"""
	Read strand info into memory
"""
lookup = {}
for r in args.b:
	col = r.strip().split("\t")
	if col[7] is "1/-1":
		lookup[col[0]] = True
"""
	Read FASTA file and reverse complement sequences with wrong strand
	Read FASTA and QUALITY file simultaneously and write output
"""
header = ""
sequence = ""
name = args.f.name.replace(".fasta", "").replace(".fna", "")
out = open("%s.cs.fasta" % name, 'w+')

for r in args.f:
	r = r.strip()
	if r.startswith(">"):
		if len(sequence) > 0:
			if header in lookup:
				sequence = reverse_complement(sequence)
			out.write("%s\n%s\n" % (header, sequence))
		header = r
		sequence = ""
	else:
		sequence += r

if len(sequence) > 0:
	if header in lookup:
		sequence = reverse_complement(sequence)
	out.write("%s\n%s\n" % (header, sequence))

if out is not None:
	out.close()