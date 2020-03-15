#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import sys
import operator
import random
import math

"""
	Picks sequences from the BLAST+ results based on a keyword(s) search.
"""
parser = argparse.ArgumentParser(description = """ 
	Picks sequences from the BLAST+ results based on a keyword(s) search.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', required = True, type = open, help = """
	specify the FASTA file that was used for the BLAST+
	""")
parser.add_argument(
	'-uc', metavar = 'UC_FILE', type = open, help = """
	specify a clustered UC file
	""")
parser.add_argument(
	'-b', metavar = 'BLAST_FILE', required = True, type = open, help = """
	specify a BLAST+ output
	""")
parser.add_argument(
	'-l', metavar = 'LOOKUP_FILE', required = False, type = open, help = """
	specify a lookup file to convert GB accession codes to taxa, if not defined, third column of BLAST+ hit description from BLAST+ output is used to search entries
	""")
parser.add_argument(
	'-i', metavar = 'IDENTITY', required = False, default = 90, type = int, help = """
	specify the minimum identity allowed for a BLAST+ hit
	""")
parser.add_argument(
	'-a', metavar = 'ALIGNMENT%', required = False, default = 90, type = int, help = """
	specify the minimum alignment length allowed for a BLAST+ hit
	""")
parser.add_argument(
	'-c', metavar = 'COUNT', required = False, default = 0, type = int, help = """
	specify the number of sequences to pick for each hit found using the keyword list
	""")
parser.add_argument(
	'-d', metavar = 'DESCRIPTION', required = False, default = False, type = bool, help = """
	include BLAST+ description/taxonomy information in the FASTA headers
	""")
parser.add_argument(
	'-s', metavar = 'KEYWORD', required = False, default = 'mucoromycotina', type = str, help = """
	specify keywords to be searched among the BLAST+ results, use commas to separate keywords
	""")
	
args = parser.parse_args()

keywords = args.s.lower().split(",")
lookup = {}

if args.l:
	console.log("Building lookup hash\n")
	with open(args.l.name) as infile:
		for f in infile:
			col = f.strip().split("\t")
			lookup[int(col[0])] = col[1]
	console.log("Found %s lookup\n" % len(lookup))

found = {}
count = {}
ucs = {}

if args.uc is not None:
	console.log("Looking up clusters to map sequences back to individuals\n")
	for f in args.uc:
		f = f.strip()
		col = f.split("\t")
		if col[0] == "S":
			_cluster = col[8]
			if _cluster not in ucs:
				ucs[_cluster] = [_cluster]
		elif col[0] == "H":
			_member = col[8]
			_cluster = col[9]
			# if for some reason centroids and members are outputed in random order, add centroid
			if _cluster not in ucs:
				ucs[_cluster] = [_cluster]
			ucs[_cluster].append(_member)
	console.log("Found %s clusters\n" % len(ucs))

console.log("Parsing BLAST output\n")
with open(args.b.name) as infile:
	for f in infile:
		col = f.strip().split("\t")
		if len(col) > 10:
			mlen = min(int(col[12]), int(col[13]))
			alen = int(col[6])
			identity = float(col[4])
			if identity >= args.i and alen >= mlen * (float(args.a) / 100):
				try:
					col[1] = int(col[1])
				except ValueError:
					try:
						col[1] = int(col[1].split("|")[1])
					except ValueError:
						col[1] = col[1]
					except IndexError:
						col[1] = col[1]
				if col[1] in lookup:
					comp = lookup[col[1]].lower()
				else:
					comp = col[2].lower()
				for key in keywords:
					if key in comp:
						if args.c and args.c > 0:
							if key not in count:
								count[key] = 1
							else:
								count[key] += 1
							if count[key] > args.c:
								continue
						if args.d:
							if args.l:
								found[col[0]] = comp	
							else:
								found[col[0]] = col[2]
						else:
							found[col[0]] = key
						if args.uc and col[0] in found and col[0] in ucs:
							# also count cluster members
							for t in ucs[col[0]]:
								found[t] = found[col[0]]
	console.log("Found %d hits\n" % len(found))

selected = 0
write = False
with open(args.f.name) as infile:
	for f in infile:
		f = f.strip()
		if len(f) > 0 and f[0] == ">":
			header = f[1:].split(" ")[0].split("\t")[0]
			if header in found:
				write = True
				sys.stdout.write("%s\t%s\n" % (f, found[header]))
				selected += 1
			else:
				write = False
		elif write:
			sys.stdout.write("%s\n" % f)

console.log("Found %d sequences\n" % selected)
