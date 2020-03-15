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
	Picks sequences from the cluster.
"""
parser = argparse.ArgumentParser(description = """ 
	Picks sequences from the cluster.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', required = True, type = open, help = """
	specify a FASTA file that was used for the clustering
	""")
parser.add_argument(
	'-c', metavar = 'CLUSTER_FILE', required = True, type = open, help = """
	specify a clustered UC file that was the result of the clustering
	""")
parser.add_argument(
	'-s', metavar = 'SEQUENCES_CLUSTER', required = False, default = 4, type = int, help = """
	specify the number of sequences to pick from each cluster 
	""")
parser.add_argument(
	'-m', metavar = 'MINIMUM_SIZE', required = False, default = 100, type = int, help = """
	specify the minimum cluster size (number of sequences) to be included in the search
	""")
parser.add_argument(
	'-blast', metavar = 'BLAST_FILE', required = False, type = open, help = """
	specify a BLAST+ file to provide hit description for the FASTA headers
	""")
parser.add_argument(
	'-l', metavar = 'LOOKUP_FILE', required = False, type = open, help = """
	specify a lookup file to convert GB accession codes to taxa, if not defined, the third column of the BLAST+ hit description from the BLAST+ output is used to search entries
	""")
parser.add_argument(
	'-centroid', metavar = 'CENTROID', required = False, default = True, type = bool, help = """
	specify whether to pick centroids
	""")
parser.add_argument(
	'-group', metavar = 'GROUP', required = False, default = True, type = bool, help = """
	groups clusters together
	""")
parser.add_argument(
	'-sort', metavar = 'SORT', required = False, default = True, type = bool, help = """
	sorts clusters in descending order by cluster size (number of sequences)
	""")
parser.add_argument(
	'-pick', metavar = 'SELECTION', required = False, default = 3, type = int, help = """
	specify selection method:
		1 - sequentialy (selected from the first to last)
		2 - equally divided
		3 - randomly
	""")
parser.add_argument(
	'-ml', metavar = 'MINIMUM_LENGTH', required = False, default = 100, type = int, help = """
	specify the minimum length for selected sequences
	""")
parser.add_argument(
	'-mi', metavar = 'IDENTITY', required = False, default = 97, type = int, help = """
	specify the minimum identity from the BLAST+ result
	""")
	
args = parser.parse_args()

def sort_by_values_len(dict):
    dict_len= {key: len(value) for key, value in dict.items()}
    import operator
    sorted_key_list = sorted(dict_len.items(), key=operator.itemgetter(1), reverse=True)
    sorted_dict = [{item[0]: dict[item [0]]} for item in sorted_key_list]
    return sorted_dict

# find out clusters and their members
console.log("Building hash table for clusters\n")

clusters = {}
conversion = {}
blasts = {}

if args.blast:
	with open(args.blast.name) as infile:
		for f in infile:
			f = f.strip()
			col = f.split("\t")
			mlen = min(int(col[12]), int(col[13]))
			alen = int(col[6])
			alignment = float(mlen) / alen * 100.0
			if alignment >= 100:
				alignment = 100
			# select only best hit
			if col[0] in blasts:
				if col[14] > blasts[col[0]][4]:
					blasts[col[0]] = [col[2], col[3], col[4], alignment, col[14]]
			else:
				blasts[col[0]] = [col[2], col[3], col[4], alignment, col[14]]
			
with open(args.c.name) as infile:
	for f in infile:
		f = f.strip()
		col = f.split("\t")
		if col[0] == "H" and float(col[2]) >= args.ml and float(col[3]) >= args.mi:
			if len(col) > 9:
				conversion[col[8]] = col[9]
				if col[9] in clusters:
					clusters[col[9]].append(col[8])
				else:
					clusters[col[9]] = [col[9], col[8]]
					conversion[col[9]] = col[9]

console.log("Found %d clusters\n" % len(clusters))
console.log("Remove small clusters\n")

keys = list(clusters.keys())

for k in keys:
	if len(clusters[k]) < args.m:
		clusters.pop(k)

#clusters = sort_by_values_len(clusters)
if args.sort:
	keys = sorted(clusters, key=lambda k: len(clusters[k]), reverse=True)
else:
	keys = clusters

console.log("Found %d clusters\n" % len(clusters))

select_sequences = {}

c_i = 0
for k in keys:
	c_i += 1
	size = len(clusters[k])
	select = args.s
	if select > size:
		select = size
	if args.centroid:
		select -= 1
		select_sequences[clusters[k][0]] = [c_i, size]
	
	if args.pick == 1 or size == select:
		for i in range(select):
			select_sequences[clusters[k][i + 1]] = [c_i, size]
	elif args.pick == 2 or size < select * 2:
		if args.centroid:
			select += 1
		if select > size:
			select = size
		step = int(math.ceil(size / select))
		for i in range(0, size, step):
			select_sequences[clusters[k][i]] = [c_i, size]
	else:
		while select > 0:
			rnd = random.randint(0, size - 1)
			if clusters[k][rnd] not in select_sequences:
				select -= 1
				select_sequences[clusters[k][rnd]] = [c_i, size]

console.log("Found %d sequences to be picked from %d clusters\n" % (len(select_sequences), len(clusters)))

sequences = {}
write = False
output = {}

with open(args.f.name) as infile:
	for f in infile:
		if len(f) > 0 and f[0] == ">":
			header = f[1:].strip().split("\t")[0].split(" ")[0]
			if header in select_sequences:
				cluster = select_sequences[header]
				if args.blast and header in conversion and conversion[header] in blasts:
					blast = blasts[conversion[header]]
					header_out = ">%s_C%s_S%s i%s a%s e%s %s\n" % (header, cluster[0], cluster[1], blast[2], blast[3], blast[1], blast[0])
				else:
					header_out = ">%s_C%s_S%s\n" % (header, cluster[0], cluster[1])
				if not args.group:
					sys.stdout.write(header_out)
				else:
					if cluster[0] not in output:
						output[cluster[0]] = [header_out]
						cluster_i = 0
					else:
						cluster_i = len(output[cluster[0]])
						output[cluster[0]].append(header_out)
				write = True
			else:
				write = False
		elif write:
			if not args.group:
				sys.stdout.write("%s\n" % f.strip())
			else:
				output[cluster[0]][cluster_i] += f.strip()

if args.group:
	for k in output:
		for l in output[k]:
			sys.stdout.write("%s\n" % l)
