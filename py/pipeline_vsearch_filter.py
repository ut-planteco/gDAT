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
	Converts vsearch global alignment search to match the BLAST output, adds 100% alignment columns at the end.
"""
parser = argparse.ArgumentParser(description = """
	Converts vsearch global alignment search to match the BLAST output, adds 100% alignment columns at the end.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA FILE', required = True, type = open, help = """
	FASTA input file that is used for the vsearch to count sequence positions in the file
	""")

args = parser.parse_args()
headers = {}
finished = {}
l = 0

for f in args.f:
	if f[0] == ">":
		l += 1
		if l % 25000 == 0:
			headers[f[1:].strip()] = l
sys.stderr.write("Found %d sequences\n" % l)
if l == 0:
	l = 1

start = time.time()

for line in sys.stdin:
	col = line.strip().split("\t")
		if col[0] in headers and col[0] not in finished:
			done = headers[col[0]] / l * 100.0
			finished[col[0]] = True
			diff = time.time() - start
			remaining = round(diff / done * 100.0)
			sys.stderr.write("Processed %d/%d vsearch hits (%.2f%%), %ds remaining\n" % (headers[col[0]], l, done, remaining))
		newline = [col[0], col[1].split(" ")[0], col[1], "0", col[2], "%s" % (int(col[3]) - int(col[4])), col[3], "1/1", col[6], col[7], col[8], col[9], col[7], col[9], col[3]]
		sys.stdout.write("%s\n" % ("\t".join(newline)))	
	except:
		continue
