#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import fileinput
import sys

"""
	Converts UC cluster file into a pivot table using the clustered file to provide cluster sizes and the BLAST+ file to map back individual sequences to each cluster centroid and to identify the best hit using identity and alignment information.
"""

parser = argparse.ArgumentParser(description = """
	Converts UC cluster file into a pivot table using the clustered file to provide cluster sizes and the BLAST+ file to map back individual sequences to each cluster centroid and to identify the best hit using identity and alignment information.
	""")
parser.add_argument(
	'-uc', metavar = 'CLUSTERED_UC', required = True, type = open, help = """
	specify a clustered UC file
	""")
parser.add_argument(
	'-blast', metavar = 'BLAST_FILE', required = True, type = open, help = """
	specify a BLAST+ file
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', required = False, type = open, help = """
	specify a clustered centroids fasta file to generate a nohits fasta
	""")
parser.add_argument(
	'-lookup', metavar = 'LOOKUP_FILE', required = False, type = open, help = """
	specify a LOOKUP file to convert the BLAST+ results
	""")
parser.add_argument(
	'-sort', metavar = 'SORT', required = False, type = bool, help = """
	sort clusters by number of elements (sequences) from largest to smallest
	""")
parser.add_argument(
	'-i', metavar = 'IDENTITY', required = False, type = int, default = 97, help = """
	specify minimum allowed identity in percentage, this is used to the generate nohits fasta
	""")
parser.add_argument(
	'-l', metavar = 'ALIGNMENT', required = False, type = int, default = 95, help = """
	specify minimum allowed alignment length in percentage, this is used to generate the nohits fasta
	""")
parser.add_argument(
	'-s', metavar = 'CENTROID_SEQUENCE', required = False, type = bool, default = False, help = """
	write out sequences for each centroid in the OTU table
	""")
parser.add_argument(
	'-commontaxa', metavar = 'COMMON_TAXA', type = int, default = 0, help = """
	specify as 1 to generate a common taxonomy, when using GB/NCBI BLAST+ with multiple best hits. For  a given rank at least 51 percent of the hits should have the same name to be accepted as a consensus identity; ranks are checked sequentially from lowest to highest.
	""")
	
args = parser.parse_args()

clusters = {}
cluster_ids = {}
centroids = {}
rev_centroids = {}
rows = {}
cols = {}
hits = {}
blast = {}
lookup = {}
nohits = {}
centroids_i = 0
total = 0

if args.lookup:
	for f in args.lookup:
		col = f.strip().split("\t")
		lookup[col[0]] = col[1]

i = 0
multiple_queries = {}
multiple_hits = {}
hits_commontaxa = {}
if args.commontaxa and args.lookup:
	console.log("Parsing BLAST to calculate common taxa\n")
	fh = open(args.blast.name, "r")
	for f in fh:
		f = f.strip()
		col = f.split("\t")
		if int(col[13]) - int(col[11]) < 10 and int(col[12]) - int(col[9]) >= 10:
			# if reference sequence ends, substract the missing part (query sequence runs over refernce sequence)
			mlen = min(int(col[12]), int(col[13])) - (int(col[12]) - int(col[9]))
		else:
			mlen = min(int(col[12]), int(col[13]))
		alen = int(col[6])
		if float(col[4]) >= args.i and alen >= mlen * (args.l / 100.0):
			if args.commontaxa and args.lookup:
				try:
					col[1] = int(col[1])
				except ValueError:
					try:
						col[1] = int(col[1].split("|")[1])
					except:
						col[1] = col[1]
				hit = str(col[1])
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
		for i in range(longest):
			if i == 0:
				for t in res[0]:
					unknowns.append(t)
			maxk = max(res[i], key=res[i].get)
			if res[i][maxk] > len(multiple_hits[k]) * 0.5:
				consensus.append(maxk)
			else:
				break
		if len(consensus) == 0:
			consensus.append("Unknown, multiple hits: %s" % ", ".join(unknowns))
		hits_commontaxa[k] = "; ".join(consensus)

for f in args.uc:
	col = f.strip().split("\t")
	# ignore short rows
	if len(col) < 9:
		continue
	sample = col[8].split("|")[0]
	found = False
	if col[0] == "S":
		centroids_i += 1
		centroids[col[8]] = centroids_i
		cluster = col[8] 
		if cluster not in clusters:
			clusters[cluster] = [col[8]]
		found = True
	if col[0] == "H":
		cluster = col[9]
		if cluster not in clusters:
			centroids_i += 1
			centroids[col[8]] = centroids_i
			clusters[cluster] = [col[8]]
		clusters[cluster].append(col[8])
		found = True
	if found:
		index = "%s_%s" % (cluster, sample)
		total += 1
		if index in hits:
			hits[index] += 1
		else:
			hits[index] = 1
		if sample in cols:
			cols[sample] += 1
		else:
			cols[sample] = 1
		if cluster in rows:
			rows[cluster] += 1
		else:
			rows[cluster] = 1

if args.sort:
	keys = sorted(clusters, key=lambda k: len(clusters[k]), reverse=True)
	_keys = sorted(cols, key=lambda k: cols[k], reverse=True)
else:
	keys = list(clusters.keys())
	_keys = list(cols.keys())

for f in args.blast:
	col = f.strip().split("\t")
	try:
		_id = int(col[1])
	except ValueError:
		try:
			_id = int(col[1].split("|")[1])
		except ValueError:
			_id = col[1]
		except IndexError:
			_id = col[1]
	alen = int(col[6])
	if int(col[13]) - int(col[11]) < 10 and int(col[12]) - int(col[9]) >= 10:
		# if reference sequence ends, substract the missing part (query sequence runs over refernce sequence)
		mlen = min(int(col[12]), int(col[13])) - (int(col[12]) - int(col[9]))
	else:
		mlen = min(int(col[12]), int(col[13]))
	alignment = float(mlen) / alen * 100
	if alignment > 100:
		alignment = 100
	desc = ""
	found = False
	if args.commontaxa:
		if col[0] in hits_commontaxa:
			desc = hits_commontaxa[col[0]]
		else:
			desc = "poorhit"
		if col[0] in blast:
			if col[14] > blast[col[0]][8]:
				found = True
		else:
			found = True
	else:
		if args.lookup:
			if str(_id) in lookup:
				desc = lookup[str(_id)]
		# select only best hit
		if col[0] in blast:
			if col[14] > blast[col[0]][8]:
				found = True
		else:
			found = True
	if found:
		blast[col[0]] = [str(_id), desc, col[0], col[1], col[2], col[3], col[4], str(alignment), col[14]]

# calculate nohits
for cluster in rows:
	if cluster in blast:
		if float(blast[cluster][6]) < args.i or float(blast[cluster][7]) < args.l:
			nohits[cluster] = True
	else:
		nohits[cluster] = True

sequences = {}

if args.f:
	fh = open("%s.i%d.a%d.nohits.fasta" % (args.f.name, args.i, args.l), "w+")
	write = False
	centroid = False
	for f in args.f:
		if len(f) > 0 and f[0] == ">":
			write = False
			header = f.strip()[1:]
			if header in nohits:
				write = True
				f = f.strip()
				if header in cluster_ids:
					f = "%s|%s" % (f, cluster_ids[header])
				if header in blast:
					col = blast[header]
					f = "%s|%s" % (f, "|".join([col[3], col[4], col[6], col[7]]))
				f = "%s\n" % f
			if args.s and header in centroids:
				centroid = True
			else:
				centroid = False
		elif centroid:
			if header in sequences:
				sequences[header] += f.strip()
			else:
				sequences[header] = f.strip()
		if write:
			fh.write(f)
	fh.close()

sys.stdout.write("#OTU\t")
for _key in _keys:
	sys.stdout.write("%s\t" % _key)
sys.stdout.write("Total\tID\tLookup\tCentroid\tBLAST ID\tBLAST description\tE-value\tIdentity\tAlignment\tScore\tSequence\n")
sys.stdout.write("#Total\t")
for _key in _keys:
	sys.stdout.write("%s\t" % cols[_key])
sys.stdout.write("%s\n" % total)
i = 0
for key in keys:
	i += 1
	if args.sort:
		otu = i
	else:
		otu = centroids[key]
	cluster_ids[key] = "OTU%s" % i
	sys.stdout.write("OTU%s\t" % i)
	for _key in _keys:
		index = "%s_%s" % (key, _key)
		if index in hits:
			sys.stdout.write("%s\t" % hits[index])
		else:
			sys.stdout.write("0\t")
	sys.stdout.write("%s\t" % rows[key])
	if key in blast:
		sys.stdout.write("\t".join(blast[key]))
	else:
		nohits[key] = True
		sys.stdout.write("nohit\tnohit\t%s\t\t\t\t\t\t" % key)
	if args.s and key in sequences:
		sys.stdout.write("\t%s" % sequences[key])
	sys.stdout.write("\n")