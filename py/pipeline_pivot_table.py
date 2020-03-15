#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import sys
import operator

"""
	Summarize BLAST+ results using specified identity and alignment length.
	thresholds.
"""
parser = argparse.ArgumentParser(description = """ 
	Summarize BLAST+ results using specified identity and alignment length.
	thresholds.
	""")
parser.add_argument(
	'-b', metavar = 'BLAST_FILE', required = True, type = open, help = """
	specify a BLAST+ output
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', type = open, help = """
	specify a nonclustered FASTA file to provide sequence counts per sample and to generate a nohits file
	""")
parser.add_argument(
	'-uc', metavar = 'UC_FILE', type = open, help = """
	specify a clustered UC file to map hits back to individual reads
	""")
parser.add_argument(
	'-i', metavar = 'IDENTITY', required = True, type = int, help = """
	specify a minimum identity threshold in percentage for hits
	""")
parser.add_argument(
	'-l', metavar = 'ALIGNMENT', required = True, type = int, help = """
	specify a minimum alignment threshold in percentage for hits
	""")
parser.add_argument(
	'-lbp', metavar = 'ALIGNMENT', required = False, type = int, help = """
	specify minimum hit aliginment length in base pairs to be accepted as an hit
	""")
parser.add_argument(
	'-vqs', metavar = 'VARIABLE_START', required = False, type = int, help = """
	specify the start position of the variable region for the query sequence
	""")
parser.add_argument(
	'-vqe', metavar = 'VARIABLE_END', required = False, type = int, help = """
	specify the end position of variable region for the query sequence
	""")
parser.add_argument(
	'-vs', metavar = 'VARIABLE_START', required = False, type = int, help = """
	specify the start position of the variable region for the reference sequence
	""")
parser.add_argument(
	'-ve', metavar = 'VARIABLE_END', required = False, type = int, help = """
	specify the end position of the variable region for the reference sequence
	""")
parser.add_argument(
	'-t', metavar = 'BLAST_TYPE', required = True, type = int, help = """
	specify which BLAST+ database was used:
	0 - use the VTX identifier from the end of the hit description (MaarjAM),
	1 - use the full hit description,
	2 - use the Genbank accession code (2nd column) (GB/NCBI)
	3 - use the SH identifier is used (UNITE)
	""")
parser.add_argument(
	'-tn', metavar = 'NODE_FILE', type = open, help = """
	specify a taxonomy file, where the full tree node connections of node IDs are provided to 
	build full taxonomy tree (nodes.dmp)
	""")
parser.add_argument(
	'-tt', metavar = 'TAXONOMY_FILE', type = open, help = """
	specify a taxonomy file, where for each node ID the scientific name is provided (names.dmp)
	""")
parser.add_argument(
	'-ti', metavar = 'ID_FILE', type = open, help = """
	specify a taxonomy file, where for each GenBank accession code a node ID is specified (gi_taxid_nucl.dmp)
	""")
parser.add_argument(
	'-zeros', metavar = 'ZEROS', type = int, default = 0, help = """
	specify as 1 to write out zero values, otherwise empty values are used
	""")
parser.add_argument(
	'-commontaxa', metavar = 'COMMON_TAXA', type = int, default = 0, help = """
	specify as 1 to generate a common taxonomy, when using GB/NCBI BLAST+ with multiple best hits. For a given rank at least 51 percent of the hits should have the same name to be accepted as a consensus identity; ranks are checked sequentially from lowest to highest. 
	""")
parser.add_argument(
	'-proportions', metavar = 'PROPORTIONS', type = int, default = 0, help = """
	specify as 1 when using multiple best hits, with proportions calculated based on multiple hits, where all hits for one query constitute 1 total hit
	""")
parser.add_argument(
	'-nohits', metavar = 'NOHITS', type = int, default = 0, help = """
	specify as 1 to write out a nohits file
	""")
parser.add_argument(
	'-lookup', metavar = 'LOOKUP_FILE', type = open, help = """
	specify a lookup file to convert the hit description to lookup values
	""")
	
args = parser.parse_args()

def buildTree(tax_id):
	global taxonomy_tree
	if tax_id in taxonomy_tree:
		return taxonomy_tree[tax_id]
	elif tax_id in ids:
		convert_id = ids[tax_id]
		previous_id = 1
		l = []
		while convert_id != 1 and previous_id != convert_id:
			previous_id = convert_id
			if convert_id in nodes:
				if convert_id in names:
					l.append(names[convert_id])
				else:
					l.append("root")
				convert_id = nodes[convert_id]
				if convert_id == 1:
					l.append("root")
			else:
				l.append(tax_id)
				convert_id = 1
		tmp = "; ".join(l[::-1])
		if len(tmp) == 0:
			tmp = tax_id
		taxonomy_tree[tax_id] = tmp
		return tmp
	else:
		taxonomy_tree[tax_id] = tax_id
		return taxonomy_tree[tax_id]

# TREE
ids = {}
converted_ids = {}
nodes = {}
names = {}
taxonomy_tree = {}
# BLAST
nohits = []
rows = {}
cols = {}
hits = {}
full_description = {}
total = 0
# FASTA
sequences = {}
samples = {}
sequences_total = 0
# CLUSTER
ucs = {}
centroids = {}

if args.tn and args.t == 2 and not args.lookup:
	console.log("Fetching list of GenBank IDs from BLAST\n")
	for f in args.b:
		f = f.strip()
		col = f.split("\t")
		try:
			col[1] = int(col[1])
		except ValueError:
			col[1] = int(col[1].split("|")[1])
		ids[col[1]] = True
	sorted_ids = sorted(ids)
	sorted_ids_i = 0
	console.log("Found %d GenBank IDs\n" % len(ids))
	console.log("Fetching conversion of GenBank IDs to nodes\n")
	i = 0
	for f in args.ti:
		col = f.split("\t")
		_id = int(col[0])
		if _id in ids:
			i += 1
			converted_id = int(col[1].strip())
			ids[_id] = converted_id
			converted_ids[converted_id] = True
			if i % 10 == 0:
				console.log("Found %d nodes for %d GenBank IDs\n" % (i, len(ids)))
			if i >= len(ids) or _id == sorted_ids[-1]:
				break
	console.log("Found %d nodes for %d GenBank IDs\n" % (i, len(ids)))
	console.log("Building taxonomy tree lookup\n")
	for f in args.tn:
		f = f.strip()
		col = f.split("\t|\t")
		nodes[int(col[0])] = int(col[1])
	console.log("Taxonomy tree size %d\n" % len(nodes))
	full_tree = {}
	# build lookup of full tree
	for convert_id in converted_ids:
		full_tree[convert_id] = True
		while convert_id != 1:
			if convert_id in nodes:
				convert_id = nodes[convert_id]
				full_tree[convert_id] = True
			else:
				convert_id = 1
				full_tree[convert_id] = True
	console.log("Selected %d nodes for full tree\n" % len(full_tree))
	console.log("Building taxonomy name lookup\n")
	for f in args.tt:
		f = f.replace("\t|\n", "\n").strip()
		col = f.split("\t|\t")
		if int(col[0]) in full_tree:
			if col[3] == "scientific name":
				names[int(col[0])] = col[1]
	console.log("Taxonomy name lookup size %d\n" % len(names))
	console.log("Writing out taxonomy tree\n")
	fh = open("%s.taxa" % args.b.name.replace(".blast", ""), "w+")
	for _id in sorted(ids):
		taxa = buildTree(_id)
		fh.write("%s\t%s\n" % (_id, taxa))
	fh.close()

if args.uc is not None:
	console.log("Parsing cluster file\n")
	for f in args.uc:
		f = f.strip()
		col = f.split("\t")
		# ignore short rows
		if len(col) < 9:
			continue
		if col[0] == "S":
			_cluster = col[8]
			sequences_total += 1
			index = _cluster.split("|")[0]
			if index not in samples:
				samples[index] = 1
			else:
				samples[index] += 1
			if _cluster not in ucs:
				ucs[_cluster] = [_cluster]
				centroids[_cluster] = True
		elif col[0] == "H":
			_member = col[8]
			_cluster = col[9]
			sequences_total += 1
			index = _member.split("|")[0]
			if index not in samples:
				samples[index] = 1
			else:
				samples[index] += 1
			# if for some reason centroids and members are outputed in random order, add centroid
			if _cluster not in ucs:
				ucs[_cluster] = [_cluster]
			ucs[_cluster].append(_member)

lookup = {}
if args.lookup is not None:
	for f in args.lookup:
		f = f.strip()
		col = f.split("\t")
		if len(col) > 1:
			lookup[col[0]] = col[1]

i = 0
multiple_queries = {}
multiple_hits = {}
hits_commontaxa = {}
hits_commontaxa_counted = {}
if args.proportions or args.commontaxa:
	if args.proportions:
		console.log("Parsing BLAST+ to calculate proportions\n")
	if args.commontaxa:
		console.log("Parsing BLAST+ to calculate common taxa\n")
	fh = open(args.b.name, "r")
	for f in fh:
		f = f.strip()
		col = f.split("\t")
		if int(col[13]) - int(col[11]) < 10 and int(col[12]) - int(col[9]) >= 10:
			# if reference sequence ends, substract the missing part (query sequence runs over refernce sequence)
			mlen = min(int(col[12]), int(col[13])) - (int(col[12]) - int(col[9]))
		else:
			mlen = min(int(col[12]), int(col[13]))
		alen = int(col[6])
		# remove short alignments
		if args.lbp and alen < args.lbp:
			continue
		if float(col[4]) >= args.i and alen >= mlen * (args.l / 100.0):
			if args.vs is not None and args.ve is not None and args.vs < args.ve:
				if col[10] > col[11]:
					col[10], col[11] = col[11], col[10]
				if args.vs < int(col[10]) or args.ve > int(col[11]):
					continue
			if args.vqs is not None and args.vqe is not None and args.vqs < args.vqe:
				if args.vqs < int(col[8]) or args.vqe > int(col[9]):
					continue
			if args.proportions and not args.commontaxa:
				if col[0] in multiple_queries:
					multiple_queries[col[0]] += 1
				else:
					multiple_queries[col[0]] = 1
			if args.commontaxa and args.lookup:
				try:
					col[1] = int(col[1])
				except ValueError:
					col[1] = int(col[1].split("|")[1])
				hit = str(col[1])
				if hit in lookup:
					hit = lookup[hit]
					if col[0] not in multiple_hits:
						multiple_hits[col[0]] = [hit]
					else:
						multiple_hits[col[0]].append(hit)

if args.commontaxa and args.lookup:
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

console.log("Parsing BLAST+\n")
for f in args.b:
	f = f.strip()
	col = f.split("\t")
	if int(col[13]) - int(col[11]) < 10 and int(col[12]) - int(col[9]) >= 10:
		# if reference sequence ends, substract the missing part (query sequence runs over refernce sequence)
		mlen = min(int(col[12]), int(col[13])) - (int(col[12]) - int(col[9]))
	else:
		mlen = min(int(col[12]), int(col[13]))
	alen = int(col[6])
	i += 1
	if i % 100000 == 0:
		console.log("%d/%d hits found/parsed\n" % (total, i))
	# remove short alignments
	if args.lbp and alen < args.lbp:
		continue
	if float(col[4]) >= args.i and alen >= mlen * (args.l / 100.0):
		if args.vs is not None and args.ve is not None and args.vs < args.ve:
			if col[10] > col[10]:
				col[10], col[11] = col[11], col[10]
			if args.vs < int(col[10]) or args.ve > int(col[11]):
				continue
		if args.vqs is not None and args.vqe is not None and args.vqs < args.vqe:
			if args.vqs < int(col[8]) or args.vqe > int(col[9]):
				continue
		# when using common taxa, get hit based on query sequence
		if col[0] in hits_commontaxa:
			hit = hits_commontaxa[col[0]]
			# if we have counted the hit, no need to recount the same hits
			if col[0] in hits_commontaxa_counted:
				continue
			hits_commontaxa_counted[col[0]] = True
			lookup_found = True
		else:
			if args.t == 0:
				hit = col[2].split(" ")[-1]
			elif args.t == 1:
				hit = col[2].split("| ")[-1]
			elif args.t == 2:
				hit = col[1]
				try:
					col[1] = int(col[1])
				except ValueError:
					try:
						col[1] = int(col[1].split("|")[1])
					except:
						col[1] = col[1]
				if args.tn:
					hit = buildTree(col[1])
				if args.lookup:
					hit = str(col[1])
			elif args.t == 3:
				tmp = col[1].split("|")
				if len(tmp) > 2:
					hit = tmp[2]
				else:
					# fallback, not SH entry
					hit = col[1]
			else:
				# fallback to full description
				hit = col[2]
			lookup_found = False
			if hit in lookup:
				hit = lookup[hit]
				lookup_found = True
			full_description[hit] = col[2]
			# each hit as weight of 1, if using proportions argument, divide by number of hits per sequence
		if args.proportions and col[0] in multiple_queries:
			hit_size = 1.0 / multiple_queries[col[0]]
		else:
			hit_size = 1
			# when using clustered based, add also number of sequences per cluster for a hit
		if args.uc is not None:
			if col[0] in ucs:
				_samples = ucs[col[0]]
				if not lookup_found:
					if args.t == 0:
						hit = col[2].split(" ")[-1]
					elif args.t == 1:
						hit = col[2].split("| ")[-1]
					else:
						hit = col[1]
				
				for tmp in _samples:
					sample = tmp.split("|")[0]
					if hit in rows:
						rows[hit] += hit_size
					else:
						rows[hit] = hit_size
					if sample in cols:
						cols[sample] += hit_size
					else:
						cols[sample] = hit_size
					index = "_".join([str(hit), sample])
					if index in hits:
						hits[index] += hit_size
					else:
						hits[index] = hit_size
					total += hit_size
					sequences[tmp] = True
		else:
			sample = col[0].split("|")[0]
			if hit in rows:
				rows[hit] += hit_size
			else:
				rows[hit] = hit_size
			if sample in cols:
				cols[sample] += hit_size
			else:
				cols[sample] = hit_size
			index = "_".join([hit, sample])
			if index in hits:
				hits[index] += hit_size
			else:
				hits[index] = hit_size
			total += hit_size
			sequences[col[0]] = True

console.log("%d/%d hits found/parsed\n" % (total, i))
header = None
fh_clustered = None
i = 0

if args.f is not None:
	samples = {}
	sequences_total = 0
	console.log("Counting sequences per sample\n")
	fh = None
	if args.nohits:
		console.log("Writing out nohits\n")
		fh = open("%s.i%d.a%d.nohits.fasta" % (args.b.name, args.i, args.l), "w+")
		if args.uc:
			fh_clustered = open("%s.i%d.a%d.clustered.nohits.fasta" % (args.b.name, args.i, args.l), "w+")
	write = False
	for f in args.f:
		if len(f) > 0 and f[0] == ">":
			sequences_total += 1
			if sequences_total % 100000 == 0:
				console.log("%d/%d nohits sequences found/parsed\n" % (i, sequences_total))
			index = f[1:].split("|")[0]
			header = f[1:].strip().split("\t")[0].split(" ")[0]
			if index in samples:
				samples[index] += 1
			else:
				# count also samples, which do not have any hits
				if index not in cols:
					cols[index] = 0
				samples[index] = 1
			if header not in sequences:
				i += 1				
				write = True
			else:
				write = False
		if fh is not None and write:
			fh.write(f)
			# check if it is centroid or not
			#print(header)
			if header in centroids:
				fh_clustered.write(f)
	if fh is not None:
		fh.close()
	if fh_clustered is not None:
		fh_clustered.close()
	console.log("%d/%d nohits sequences found/parsed\n" % (i, sequences_total))	

sorted_rows = sorted(rows.items(), key=operator.itemgetter(1), reverse=True)
sorted_cols = sorted(cols.items(), key=operator.itemgetter(1), reverse=True)

console.log("Writing out pivot table with results\n")
fh = open("%s.i%d.a%d.tsv" % (args.b.name, args.i, args.l), "w+")
### print header of the tsv file
fh.write("Hit/Sample\tSamples count\tTotal\t")
### go through samples as columns
for k in sorted_cols:
	fh.write("%s\t" % (k[0]))

fh.write("taxonomy\nTaxa count\t\t%d\t" % (len(rows)))
for k in sorted_cols:
	cnt = 0
	for _k in sorted_rows:
		index = "_".join([str(_k[0]), k[0]])
		if index in hits:
			cnt += 1
	fh.write("%d\t" % (cnt))

# only if fasta file is defined
if len(samples) > 0:
	fh.write("\nCleaned reads\t\t%d\t" % (sequences_total))
	for k in sorted_cols:
		if k[0] in samples:
			fh.write("%d\t" % (samples[k[0]]))
		else:
			fh.write("%s\t" % ("0" if args.zeros == 1 else ""))

fh.write("\nTotal\t%d\t%d\t" % (len(cols), total))
for k in sorted_cols:
	fh.write("%d\t" % (k[1]))

for k in sorted_rows:
	cnt = 0
	for _k in sorted_cols:
		index = "_".join([str(k[0]), _k[0]])
		if index in hits:
			cnt += 1
	fh.write("\n%s\t%d\t%d\t" % (k[0], cnt, k[1]))
	for _k in sorted_cols:
		index = "_".join([str(k[0]), _k[0]])
		if index in hits:
			if args.proportions:
				fh.write("%f\t" % (hits[index]))
			else:
				fh.write("%d\t" % (hits[index]))
		else:
			fh.write("%s\t" % ("0" if args.zeros == 1 else ""))
	if k[0] in full_description:
		if "k__" in full_description[k[0]]:
			full_description[k[0]] = "100%%|%s" % full_description[k[0]]
		fh.write("%s" % full_description[k[0]])

fh.write("\n")

if fh is not None:
	fh.close()


console.log("Writing out transposed pivot table with results\n")
fh = open("%s.i%d.a%d.tp.tsv" % (args.b.name, args.i, args.l), "w+")
### print header of the tsv file
if sequences_total > 0:
	fh.write("Sample/Hit\tTaxa count\tCleaned reads\tTotal\t")
else:
	fh.write("Sample/Hit\tTaxa count\tTotal\t")
### go through taxa as columns
for k in sorted_rows:
	fh.write("%s\t" % (k[0]))

fh.write("\nSamples count\t\t\t%d\t" % (len(cols)))

for k in sorted_rows:
	cnt = 0
	for _k in sorted_cols:
		index = "_".join([str(k[0]), _k[0]])
		if index in hits:
			cnt += 1
	fh.write("%d\t" % (cnt))

if sequences_total > 0:
	fh.write("\nTotal\t%d\t%d\t%d\t" % (len(rows), sequences_total, total))
else:
	fh.write("\nTotal\t%d\t%d\t" % (len(rows), total))

for k in sorted_rows:
	fh.write("%d\t" % (k[1]))

for k in sorted_cols:
	cnt = 0
	for _k in sorted_rows:
		index = "_".join([str(_k[0]), k[0]])
		if index in hits:
			cnt += 1
	cnt_seqs = "0" if args.zeros == 1 else ""
	if sequences_total > 0:
		if k[0] in samples:
			cnt_seqs = samples[k[0]]
		fh.write("\n%s\t%d\t%s\t%d\t" % (k[0], cnt, cnt_seqs, k[1]))
	else:
		fh.write("\n%s\t%d\t%s\t" % (k[0], cnt, k[1]))
	for _k in sorted_rows:
		index = "_".join([str(_k[0]), k[0]])
		if index in hits:
			if args.proportions:
				fh.write("%f\t" % (hits[index]))
			else:
				fh.write("%d\t" % (hits[index]))
		else:
			fh.write("%s\t" % ("0" if args.zeros == 1 else ""))

fh.write("\n")

if fh is not None:
	fh.close()
