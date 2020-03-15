#!/usr/bin/env python
from __future__ import division

import os
import argparse
import gzip
import console
import fileinput

"""
	Combines multiple BLAST+ outputs into a single BLAST+ file by selecting BLAST+ hits based on BLAST+ score. If input and output is not defined then they are taken from STDIN and STDOUT respectively. For input parameters, separate files with commas.
"""

parser = argparse.ArgumentParser(description = """
	Combines multiple BLAST+ outputs into a single BLAST+ file by selecting BLAST+ hits based on BLAST+ score. If input and output is not defined then they are taken from STDIN and STDOUT respectively. For input parameters, separate files with commas.
	""")
parser.add_argument(
	'-f', metavar = 'BLAST_FOLDER', required = False, type = str, help = """
	specify a folder location for multiple BLAST+ files, selects files with *.0.blast extension
	""")
parser.add_argument(
	'-i', metavar = 'BLAST_FILES', required = False, type = str, help = """
	specify multiple BLAST+ files, separating them by commas
	""")
parser.add_argument(
	'-allhits', metavar = '[0|1]', required = False, type = int, default = 0, help = """
	specify whether to select all hits, otherwise only the best hit is selected for each query
	""")
	
args = parser.parse_args()

hits = {}

def packedFile(filename):
	ext = filename.split(".")[-1]
	if ext.lower() == "gz":
		return True
	else:
		return False

if args.f:
	files = []
	for filename in os.listdir(args.f):
		# find BLAST output files and do not use original BLAST file in the search
		if filename.endswith(".blast") and not filename.endswith(".fasta.blast"):
			files.append("%s/%s" % (args.f, filename))
	args.i = ",".join(files)
if args.i:
	files = args.i.split(",")
	for f in files:
		# ignore combined BLAST file
		if "all.blast" in f:
			continue
		if packedFile(f):
			fh = gzip.open(f, "r")
		else:
			fh = open(f, "r")
		for line in fh:
			col = line.strip().split("\t")
			if len(col) < 15:
				continue
			if args.allhits:
				print(line.strip())
				continue
			seq = col[0]
			score = col[14]
			if seq in hits:
				comp_score = hits[seq][14]
				if score > comp_score:
					hits[seq] = col
			else:
				hits[seq] = col
		if fh:
			fh.close()

	if not args.allhits:
		for key in hits:
			print("\t".join(hits[key]))

else:

	for line in fileinput.input():
		col = line.strip().split("\t")
		if len(col) < 15:
			continue
		if args.allhits:
			print(line.strip())
			continue
		seq = col[0]
		score = col[14]
		if seq in hits:
			comp_score = hits[seq][14]
			if score > comp_score:
				hits[seq] = col
		else:
			hits[seq] = col
	
	if not args.allhits:
		for key in hits:
			print("\t".join(hits[key]))
