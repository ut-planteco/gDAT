#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import sys
import operator
import re
import gzip

"""
	Splits interleaved reads into two separate files to be handled with vsearch paired-end merger.
"""
parser = argparse.ArgumentParser(description = """
	Splits interleaved reads into two separate files to be handled with vsearch paired-end merger.
	""")
parser.add_argument(
	'-i', metavar = 'FASTQ FILE', required = True, type = open, help = """
	FASTQ input file where forward reads and reverse reads are at odd and even positions respectively in the file
	""")

args = parser.parse_args()
headers = {}
finished = {}
l = 0

for f in args.i:
	l += 1
sys.stderr.write("Found 2x%d paired-end reads\n" % (l / 8))

f1 = open("%s.r1.fastq" % args.i.name, "w")
f2 = open("%s.r2.fastq" % args.i.name, "w")

l = 0
with open(args.i.name) as f:
	for r in f:
		if l % 8 < 4:
			f1.write(r)
		else:
			f2.write(r)
		l += 1
f1.close()
f2.close()