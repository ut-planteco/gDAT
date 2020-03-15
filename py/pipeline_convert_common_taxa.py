#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import fileinput
import sys

"""
	Converts GB/NCBI BLAST+ output to include common taxonomy hits.
"""

parser = argparse.ArgumentParser(description = """
	Converts GB/NCBI BLAST+ output to include common taxonomy hits.
	""")
parser.add_argument(
	'-b', metavar = 'BLAST_FILE', required = True, type = open, help = """
	specify a BLAST+ file containing multiple best hits
	""")
parser.add_argument(
	'-lookup', metavar = 'LOOKUP_FILE', required = True, type = open, help = """
	specify a LOOKUP file to convert GB/NCBI accession codes to taxonomy
	""")
parser.add_argument(
	'-i', metavar = 'IDENTITY', required = True, type = int, help = """
	specify the minimum hit identity in percentage
	""")
parser.add_argument(
	'-l', metavar = 'ALIGNMENT', required = True, type = int, help = """
	specify the minimum allowed hit aliginment length in percentage
	""")
	
args = parser.parse_args()

lookup = {}
for f in args.lookup:
	f = f.strip()
	col = f.split("\t")
	if len(col) > 1:
		lookup[col[0]] = col[1]

i = 0
multiple_queries = {}
multiple_hits = {}
last_consensus_counts = {}
besthits = {}
hits_commontaxa = {}
hits_commontaxa_counted = {}
console.log("Parsing BLAST+ to calculate common taxa\n")
fh = open(args.b.name, "r")
for f in fh:
	f = f.strip()
	col = f.split("\t")
	if len(col) < 14:
		continue
	if int(col[13]) - int(col[11]) < 10 and int(col[12]) - int(col[9]) >= 10:
		# if reference sequence ends, substract the missing part (query sequence runs over refernce sequence)
		mlen = min(int(col[12]), int(col[13])) - (int(col[12]) - int(col[9]))
	else:
		mlen = min(int(col[12]), int(col[13]))
	alen = int(col[6])
	if float(col[4]) >= args.i and alen >= mlen * (args.l / 100.0):
		try:
			col[1] = int(col[1])
		except ValueError:
			try:
				col[1] = int(col[1].split("|")[1])
			except:
				col[1] = col[1]
		hit = str(col[1])
		if col[0] in besthits:
			if int(col[14]) > int(besthits[col[0]][14]):
				besthits[col[0]] = col
		else:
			besthits[col[0]] = col
		if hit in lookup:
			hit = lookup[hit]
			if col[0] not in multiple_hits:
				multiple_hits[col[0]] = [hit]
			else:
				multiple_hits[col[0]].append(hit)

console.log("Building common taxa\n")
for k in multiple_hits:
	res = {}
	consensus = []
	tmp = multiple_hits[k]
	i = 0
	longest = 0
	# use nested dictionary
	for taxa in tmp:
		t = taxa.split("; ")
		i = 0
		for name in t:
			if i not in res:
				res[i] = {}
			if name in res[i]:
				res[i][name] += 1
			else:
				res[i][name] = 1
			i += 1
			if i > longest:
				longest = i
	# calculate longest consensus
	unknowns = []
	last_consensus_count = 0
	for i in range(longest):
		if i == 0:
			for t in res[0]:
				unknowns.append(t)
		maxk = max(res[i], key=res[i].get)
		if res[i][maxk] > len(multiple_hits[k]) * 0.5:
			last_consensus_count = res[i][maxk]
			consensus.append(maxk)
		else:
			break
	if len(consensus) == 0:
		consensus.append("Unknown, multiple hits: %s" % ", ".join(unknowns))
	last_consensus_counts[k] = last_consensus_count
	hits_commontaxa[k] = "; ".join(consensus)

console.log("Converting BLAST+\n")
for f in args.b:
	f = f.strip()
	col = f.split("\t")
	if len(col) < 14:
		continue
	if int(col[13]) - int(col[11]) < 10 and int(col[12]) - int(col[9]) >= 10:
		# if reference sequence ends, substract the missing part (query sequence runs over refernce sequence)
		mlen = min(int(col[12]), int(col[13])) - (int(col[12]) - int(col[9]))
	else:
		mlen = min(int(col[12]), int(col[13]))
	alen = int(col[6])
	if float(col[4]) >= args.i and alen >= mlen * (args.l / 100.0):
		# when using common taxa, get hit based on query sequence
		if col[0] in hits_commontaxa:
			hit = hits_commontaxa[col[0]]
			# if we have counted the hit, no need to recount the same hits
			if col[0] in hits_commontaxa_counted:
				continue
			hits_commontaxa_counted[col[0]] = True
			col = besthits[col[0]]
			col[1] = str(col[1])
			col[2] = hit
			if last_consensus_counts[col[0]] > 0:
				perc = (last_consensus_counts[col[0]] / float(len(multiple_hits[col[0]]))) * 100.0
			else:
				perc = 0.0
			col[3] = "Matches/Total: %s/%s (%.2f%%)" % (last_consensus_counts[col[0]], len(multiple_hits[col[0]]), perc)
			print("\t".join(col))