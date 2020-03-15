#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import fileinput
import sys

"""
	Converts a UC cluster file into legacy BLASTCLUST format, where for each cluster all 
	the entries are shown in single line separated by space.
"""

parser = argparse.ArgumentParser(description = """
	Converts a UC cluster file into legacy BLASTCLUST format, where for each cluster all 
	the entries are shown in single line separated by space.
	""")
parser.add_argument(
	'-uc', metavar = 'CLUSTER_FILE', required = True, type = open, help = """
	specify a clustered UC file to be converted
	""")
parser.add_argument(
	'-sort', metavar = 'SORT', required = False, type = bool, help = """
	sort clusters by number of elements (sequences) from largest to smallest
	""")
	
args = parser.parse_args()

clusters = {}
centroids = {}
centroids_i = 0

for f in args.uc:
	col = f.strip().split("\t")
	# ignore short rows
	if len(col) < 9:
		continue
	if col[0] == "S":
		centroids_i += 1
		centroids[col[8]] = centroids_i 
		if col[8] not in clusters:
			clusters[col[8]] = [col[8]]
	if col[0] == "H":
		if col[9] not in clusters:
			centroids_i += 1
			centroids[col[9]] = centroids_i
			clusters[col[9]] = [col[9]]
		clusters[col[9]].append(col[8])

keys = list(clusters.keys())

if args.sort:
	keys = sorted(clusters, key=lambda k: len(clusters[k]), reverse=True)

i = 0
for key in keys:
	if len(clusters[key]) > 1:
		if args.sort:
			i += 1
			sys.stdout.write("OTU%s\t%s\n" % (i, " ".join(clusters[key])))
		else:
			sys.stdout.write("OTU%s\t%s\n" % (centroids[key], " ".join(clusters[key])))

