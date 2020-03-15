#!/usr/bin/env python
from __future__ import division

import os
import sys
import argparse
import fileinput

"""
	Maps back chimeric occurences for indvidual reads in the clustered file.
"""

parser = argparse.ArgumentParser(description = """
	Maps back chimeric occurences for indvidual reads in the clustered file.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', type = open, required = True, help = """
	specify a chimera free FASTA file that was written by vsearch
	""")
parser.add_argument(
	'-uc', metavar = 'CLUSTER_FILE', type = open, help = """
	specify the clustered UC file to map back chmeric counts to individual reads
	""")
	
args = parser.parse_args()

headers = {}
count = 0

sys.stderr.write("Recalculating non-chimera count based on clustered information\n")

for f in args.f:
	f = f.strip()
	if len(f) > 0 and f[0] == ">":
		headers[f[1:].strip()] = True
		count += 1

for f in args.uc:
	f = f.strip()
	col = f.split("\t")
	if len(col) > 9 and col[9] in headers:
		count += 1

sys.stderr.write("Found %s sequences with non-chimeras\n" % count)