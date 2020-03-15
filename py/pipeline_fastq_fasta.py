#!/usr/bin/env python
from __future__ import division

import os
import sys
import argparse
import fileinput

"""
	Converts FASTQ file to FASTA or FASTA+QUAL files.
"""

parser = argparse.ArgumentParser(description = """
	Converts FASTQ file to FASTA or FASTA+QUAL files.
	""")
parser.add_argument(
	'-fq', metavar = 'FASTQ_FILE', type = open, required = True, help = """
	specify a FASTQ file
	""")
parser.add_argument(
	'-output_qual', metavar = 'QUAL_FILE', type = bool, default = False, help = """
	specify whether to generate a QUALITY file with the FASTA file
	""")
parser.add_argument(
	'-phred', metavar = 'PHRED', type = int, default = 33, help = """
	specify the PHRED value for the FASTQ file, 33 or 64
	""")

args = parser.parse_args()

# define default variables
header = ""
sequence = ""
quality = ""
qf = None

name = args.fq.name.replace(".fastq", "").replace(".fq", "")
out = open("%s.fasta" % name, 'w+')
qual_out = None
if args.output_qual:
	qual_out = open("%s.qual" % name, 'w+')

i = 0
firstline = True
for f in args.fq:
	if firstline:
		firstline = False
		if len(f) == 0 or f[0] != "@":
			sys.stderr.write("FASTQ file is corrupt or wrong file is used. First line has to start with '@'\n")
			sys.exit(1)
	# if we have non-standard FASTQ, reset counter if the line is small (not sequence nor quality) and contains @
	if i % 4 != 0 and len(f) > 0 and len(f) < 100 and f[0] == "@":
		i = 0
	f = f.strip()
	if i % 4 == 0:
		header = f.replace("@", ">")
	elif i % 4 == 1:
		# some of the tools do not allow reading dashes for the sequence, replace it with N
		out.write("%s\n%s\n" % (header, f.replace("-", "N")))
	elif i % 4 == 3 and qual_out is not None:
		tmp = []
		for c in f:
			tmp.append("%s" % (ord(c) - args.phred))
		qual_out.write("%s\n%s\n" % (header, " ".join(tmp)))
	i += 1
    
if out is not None:
	out.close()

if qual_out is not None:
	qual_out.close()
