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
	Builds TAXONOMY tree based on the BLAST+ results using GB/NCBI taxonomy and node files.
"""
parser = argparse.ArgumentParser(description = """
	Builds TAXONOMY tree based on the BLAST+ results using GB/NCBI taxonomy and node files.
	""")
parser.add_argument(
	'-b', metavar = 'BLAST FILE', required = True, type = open, help = """
	BLAST+ tabulated output that was generated by BLAST+ using custom output
	""")
parser.add_argument(
	'-ti', metavar = 'ID FILE', required = True, type = open, help = """
	Taxonomy file where for each GenBank accession ID lookup node ID is specified (gi_taxid_nucl.dmp or gi_taxid_nucl.dmp.gz)
	""")
parser.add_argument(
	'-tt', metavar = 'TAXONOMY FILE', required = True, type = open, help = """
	Taxonomy file where for each node ID scientific name is provided (names.dmp or names.dmp.gz)
	""")
parser.add_argument(
	'-tn', metavar = 'NODE FILE', required = True, type = open, help = """
	Taxonomy file where full tree node connections are provided to 
	build full taxonomy tree (nodes.dmp or nodes.dmp.gz)
	""")
parser.add_argument(
	'-simplify', metavar = '[0|1]', required = False, type = int, default = 0, help = """
	Use a simplified format of the taxonomy tree by only retrieving the following ranks: 
	kingdom, phylum, class, order, family, genus, species
	""")
parser.add_argument(
	'-ranks', metavar = '[0|1]', required = False, type = int, default = 0, help = """
	For each category shows rank in parenthesis
	""")
	
args = parser.parse_args()

ids = {}
converted_ids = {}
nodes = {}
ranks = {}
names = {}
taxonomy_tree = {}

def packedFile(filename):
	ext = filename.split(".")[-1]
	if ext.lower() == "gz":
		return True
	else:
		return False	

def buildTree(tax_id):
	global taxonomy_tree
	filtered_categories = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
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
					if convert_id in ranks:
						rank = ranks[convert_id]
					else:
						rank = "no rank"
					l.append(names[convert_id] + " (" + rank + ")")
				else:
					l.append("root")
				convert_id = nodes[convert_id]
				if convert_id == 1:
					l.append("root")
			else:
				l.append(tax_id)
				convert_id = 1
		modified = []
		if args.simplify:
			for taxa in l:
				if not isinstance(taxa, int) and (taxa.split(" (")[-1].split(")")[0] in filtered_categories or "Archaea" in taxa):
					modified.append(taxa)
			l = modified
		
		tmp = "; ".join(map(str, l[::-1]))
		if args.ranks == 0:
			tmp = re.sub(r'\([^)]*\)', '', tmp).replace(" ;", ";")
		if len(tmp) == 0:
			tmp = tax_id
		taxonomy_tree[tax_id] = tmp
		return tmp
	else:
		taxonomy_tree[tax_id] = tax_id
		return taxonomy_tree[tax_id]

console.log("Fetching list of GenBank IDs from BLAST+\n")
for f in args.b:
	f = f.strip()
	col = f.split("\t")
	try:
		col[1] = int(col[1])
	except ValueError:
		col[1] = int(col[1].split("|")[1])
	except IndexError:
		# probably empty line, ignore
		continue
	ids[col[1]] = True

sorted_ids = sorted(ids)
sorted_ids_i = 0

console.log("Found %d GenBank IDs\n" % len(ids))

if len(ids) == 0:
	console.log("BLAST+ file is empty or wrong format used")
	sys.exit(1)

console.log("Fetching conversion of GenBank IDs to nodes\n")
i = 0
step = int(len(ids) / 100)
if step < 100:
	step = 100
if packedFile(args.ti.name):
	fti = gzip.open(args.ti.name, "r")
else:
	fti = open(args.ti.name, "r")
for f in fti:
	if packedFile(args.ti.name):
		col = f.decode().split("\t")
	else:
		col = f.split("\t")
	try:
		_id = int(col[0])
	except:
		# skip, probably wrong format?
		continue
	if _id in ids:
		i += 1
		try:
			converted_id = int(col[1].strip())
		except:
			# skip, probably wrong format?
			continue
		ids[_id] = converted_id
		converted_ids[converted_id] = True
		if i % step == 0:
			console.log("Found %d nodes for %d GenBank IDs\n" % (i, len(ids)))
		if i >= len(ids) or _id == sorted_ids[-1]:
			break
fti.close()
console.log("Found %d nodes for %d GenBank IDs\n" % (i, len(ids)))

console.log("Building taxonomy tree lookup\n")
if packedFile(args.tn.name):
	ftn = gzip.open(args.tn.name, "r")
else:
	ftn = open(args.tn.name, "r")
for f in ftn:
	f = f.strip()
	col = f.split("\t|\t")
	try:
		nodes[int(col[0])] = int(col[1])
		ranks[int(col[0])] = col[2]
	except:
		# skip, probably wrong format?
		continue
ftn.close()
console.log("Taxonomy tree size %d\n" % len(nodes))

full_tree = {}
# build lookup of full tree
for convert_id in converted_ids:
	try:
		full_tree[convert_id] = True
		while convert_id != 1:
			if convert_id in nodes:
				convert_id = nodes[convert_id]
				full_tree[convert_id] = True
			else:
				convert_id = 1
				full_tree[convert_id] = True
	except:
		continue

console.log("Selected %d nodes for full tree\n" % len(full_tree))

console.log("Building taxonomy name lookup\n")
if packedFile(args.tt.name):
	ftt = gzip.open(args.tt.name, "r")
else:
	ftt = open(args.tt.name, "r")
for f in ftt:
	f = f.replace("\t|\n", "\n").strip()
	col = f.split("\t|\t")
	try:
		if int(col[0]) in full_tree:
			if col[3] == "scientific name":
				names[int(col[0])] = col[1]
	except:
		continue
		
ftt.close()
console.log("Taxonomy name lookup size %d\n" % len(names))

for _id in sorted(ids):
	taxa = buildTree(_id)
	sys.stdout.write("%s\t%s\n" % (_id, taxa))
