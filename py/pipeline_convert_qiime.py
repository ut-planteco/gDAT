#!/usr/bin/env python
from __future__ import division

import os
import sys
import argparse

"""
	Converts FASTA file header formats. If it is in QIIME v1 format, (_) underscores will be converted to (|) pipe to work with the gDAT pipeline. If the script finds | it will convert the FASTA to be compatible with QIIME v1.
"""
parser = argparse.ArgumentParser(description = """
	Converts FASTA file header formats. If it is in QIIME v1 format, (_) underscores will be converted to (|) pipe to work with the gDAT pipeline. If the script finds | it will convert the FASTA to be compatible with QIIME v1.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA', required = True, type = open, help = """
	specify a FASTA file to be converted between different pipeline formats
	""")
	
args = parser.parse_args()

qiime = False
first = True
for f in args.f:
	f = f.strip()
	if first:
		if f.find("_") < f.find("|"):
			qiime = True
		first = False
	if len(f) > 0 and f[0] == ">":
		if qiime:
			sys.stdout.write(f.replace("_", "|"))
		else:
			sys.stdout.write(f.replace("|", "_"))
	else:
		sys.stdout.write(f)