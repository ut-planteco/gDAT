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
	Picks sequences from BLAST+ results based on keyword(s) search.
"""
parser = argparse.ArgumentParser(description = """ 
	Picks sequences from BLAST+ results based on keyword(s) search.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', required = True, type = open, help = """
	specify a FASTA file that was used for the BLAST+
	""")
parser.add_argument(
	'-b', metavar = 'BLAST_FILE', required = True, type = open, help = """
	specify a BLAST+ output
	""")
parser.add_argument(
	'-uc', metavar = 'CLUSTER_FILE', required = False, type = open, help = """
	specify a cluster information file, which will be used to map BLAST+ hits back
	to individual sequences
	""")
parser.add_argument(
	'-i', metavar = 'IDENTITY', required = False, default = 97, type = int, help = """
	specify the minimum identity allowed for a BLAST+ hit
	""")
parser.add_argument(
	'-a', metavar = 'ALIGNMENT', required = False, default = 95, type = int, help = """
	specify the minimum alignment length allowed for a BLAST+ hit
	""")
parser.add_argument(
	'-t', metavar = 'BLAST_TYPE', required = True, type = int, help = """
	specify BLAST+ database used:
		0 - VTX identifier is used (MaarjAM),
		1 - full description of the hit is used,
		2 - SH identifier is used (UNITE)
	""")
parser.add_argument(
	'-pick', metavar = 'SELECTION_MODE', required = False, default = 3, type = int, help = """
	specify how to selecte sequences:
		1 - sequentialy (selected from the first to last)
		2 - equally divided
		3 - randomly
	""")
parser.add_argument(
	'-ml', metavar = 'MINIMUM_LENGTH', required = False, default = 100, type = int, help = """
	specify the minimum acceptable sequence length (100 bases)
	""")
parser.add_argument(
	'-s', metavar = 'MINIMUM_SEQUENCES', required = False, default = 4, type = int, help = """
	specify number of sequences to pick
	""")
parser.add_argument(
	'-group', metavar = 'GROUP', required = False, default = True, type = bool, help = """
	group clusters together
	""")
parser.add_argument(
	'-sort', metavar = 'SORT', required = False, default = True, type = bool, help = """
	sort clusters in descending order by cluster size (number of sequences)
	""")
parser.add_argument(
	'-persample', metavar = 'PICK_PER_SAMPLE', required = False, default = False, type = bool, help = """
	specify whether to pick representative sequences per sample per taxin
	""")
parser.add_argument(
	'-outputfasta', metavar = 'OUTPUT_FASTA', required = False, default = False, type = bool, help = """
	output in FASTA format, otherwise generates a tabulated file
	""")
	
args = parser.parse_args()

found = {}
hits = {}
seqs = {}
count = {}
hits_i = 0
ucs = {}

if args.uc is not None:
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

console.log("Parsing BLAST output")
with open(args.b.name) as infile:
	for f in infile:
		col = f.strip().split("\t")
		if len(col) > 10 and int(col[12]) > args.ml:
			if int(col[13]) - int(col[11]) < 10 and int(col[12]) - int(col[9]) >= 10:
				# if reference sequence ends, substract the missing part (query sequence runs over refernce sequence)
				mlen = min(int(col[12]), int(col[13])) - (int(col[12]) - int(col[9]))
			else:
				mlen = min(int(col[12]), int(col[13]))
			alen = int(col[6])
			identity = float(col[4])
			alignment = alen / float(mlen) * 100.0
			if alignment >= 100:
				alignment = 100
			if identity >= args.i and alignment >= args.a:
				hits_i += 1
				if args.t == 0:
					hit = col[2].split(" ")[-1]
				elif args.t == 1:
					hit = col[2].split("| ")[-1]
				elif args.t == 2:
					tmp = col[1].split("|")
					if len(tmp) > 2:
						hit = tmp[2]
					else:
						# fallback, not SH entry
						hit = col[1]
				else:
					# fallback to full description
					hit = col[2]
				sample = col[0].split("|")[0]
				index = hit
				if args.persample:
					index = "%s_%s" % (hit, sample)
				if index not in hits:
					hits[index] = []
				hits[index].append([col[0], sample, hit, identity, alignment, col[3], col[12], int(col[14])])
				# populate also clustered hits
				if args.uc and col[0] in ucs:
					for t in ucs[col[0]]:
						hits_i += 1
						if args.persample:
							sample = t.split("|")[0]
							index = "%s_%s" % (hit, sample)
							if index not in hits:
								hits[index] = []
						hits[index].append([t, sample, hit, identity, alignment, col[3], col[12], int(col[14])]) 
	console.log("Found %d hits" % hits_i)

# sort hits by abundance
if args.sort:
	keys = sorted(hits, key=lambda k: len(hits[k]), reverse=True)
else:
	keys = hits

select_sequences = {}

c_i = 0
for k in keys:
	col = hits[k]
	# sort hits by BLAST score so largest scores are always the top ones
	col.sort(key=lambda x: x[7], reverse=True)
	select = args.s
	size = len(col)
	if select > size:
		select = size
	if args.pick == 1:
		for i in range(min(size, select)):
			select_sequences[hits[k][i][0]] = [k, i]
	elif args.pick == 2:
		step = int(math.ceil(size / select))
		if step == 0:
			step = 1
		for i in range(0, args.s, step):
			select_sequences[hits[k][i][0]] = [k, i]
	else:
		while select > 0:
			rnd = random.randint(0, size - 1)
			if hits[k][rnd][0] not in select_sequences:
				select -= 1
				select_sequences[hits[k][rnd][0]] = [k, rnd]

console.log("Found %d sequences to be picked from %d hits\n" % (len(select_sequences), hits_i))
if not args.outputfasta:
	sys.stdout.write("header\tsample\thit\tidentity\talignment\te-value\tlength\tblast_score\tsequence\n")

sequences = {}
write = False
output = {}
with open(args.f.name) as infile:
	for f in infile:
		if len(f) > 0 and f[0] == ">":
			header = f[1:].strip().split("\t")[0].split(" ")[0]
			if header in select_sequences:
				ids = select_sequences[header]
				hit = []
				for v in hits[ids[0]][ids[1]]:
					hit.append(str(v))
				hit = "\t".join(hit)
				if args.outputfasta:
					header_out = ">%s\n" % hit
				else:
					header_out = "%s\t" % hit
				if not args.group:
					sys.stdout.write(header_out)
				else:
					if ids[0] not in output:
						output[ids[0]] = [header_out]
						cluster_i = 0
					else:
						cluster_i = len(output[ids[0]])
						output[ids[0]].append(header_out)
				write = True
			else:
				write = False
		elif write:
			if not args.group:
				sys.stdout.write("%s\n" % f.strip())
			else:
				output[ids[0]][cluster_i] += f.strip()

if args.group:
	for k in output:
		for l in output[k]:
			sys.stdout.write("%s\n" % l)
