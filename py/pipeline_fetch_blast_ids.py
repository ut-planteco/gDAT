#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import fileinput
import sys

"""
	This script helps to retrieve reference IDs from the second column of a BLAST+ output file. These values are used to fetch sequences from the reference database, and the output is written to a FASTA file. Input is read from STDIN or by defining files.
"""

parser = argparse.ArgumentParser(description = """
	This script helps to retrieve reference IDs from the second column of a BLAST+ output file. These values are used to fetch sequences from the reference database, and the output is written to a FASTA file. Input is read from STDIN or by defining files.
	""")
parser.add_argument(
	'-i', metavar = 'BLAST_FILES', required = False, type = str, help = """
	specify multiple BLAST+ files, separating them by commas
	""")
	
args = parser.parse_args()

hits = {}
lookup = {}

if args.i:
	files = args.i.split(",")
	for f in files:
		with open(f) as infile:
			for line in infile:
				col = line.strip().split("\t")
				if len(col) < 15:
					continue
				if col[1] not in lookup:
					lookup[col[1]] = True
else:
	for line in fileinput.input():
		col = line.strip().split("\t")
		if len(col) < 15:
			continue
		if col[1] not in lookup:
			lookup[col[1]] = True

for key in lookup:
	col = key.split("|")
	sys.stdout.write("%s\n" % key)