#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import sys
import operator
import re
import gzip
import time

"""
	Filters BLAST output to conserv space based on identity threshold and keeps track how many sequences have been parsed.
"""
parser = argparse.ArgumentParser(description = """
	Filters BLAST output to conserv space based on identity threshold and keeps track how many sequences have been parsed.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA FILE', required = True, type = open, help = """
	FASTA input file that is used for the BLAST to count sequence positions in the file
	""")
parser.add_argument(
	'-i', metavar = 'IDENTITY', required = False, type = int, default = 95, help = """
	specify minimum allowed BLAST identity to output hits
	""")
	
args = parser.parse_args()
headers = {}
l = 0

for f in args.f:
	if f[0] == ">":
		l += 1
		if l % 25000 == 0:
			headers[f[1:].strip()] = l
sys.stderr.write("Found %d sequences\n" % l)
if l == 0:
	l = 1

last_l = 0
start = time.time()

for line in sys.stdin:
	col = line.split("\t")
	try:
		if col[0] in headers and headers[col[0]] > last_l:
			done = headers[col[0]] / l * 100.0
			last_l = headers[col[0]]
			diff = time.time() - start
			remaining = round(diff / done * 100.0 - diff)
			sys.stderr.write("Processed %d/%d BLAST hits (%.2f%%), %ds remaining\n" % (headers[col[0]], l, done, remaining))
		if float(col[4]) >= args.i:
			sys.stdout.write(line)
	except:
		continue
