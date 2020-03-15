#!/usr/bin/env python
from __future__ import division

import os
import sys
import argparse
import fileinput

"""
	Converts FASTA/FASTA+QUAL files to FASTQ. Input is taken from a FASTA file, and the output filename is generated from the input filename, with the FASTA extension swapped to FASTQ.
"""

parser = argparse.ArgumentParser(description = """
	Converts FASTA/FASTA+QUAL files to FASTQ. Input is taken from a FASTA file, and the output filename is generated from the input filename, with the FASTA extension swapped to FASTQ.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', type = open, required = True, help = """
	specify a FASTA file
	""")
parser.add_argument(
	'-qf', metavar = 'QUAL_FILE', type = open, help = """
	specify a QUALITY file (optional)
	""")
parser.add_argument(
	'-dq', metavar = 'QUAL_FILE', type = int, default = 40, help = """
	specify the default quality value if a QUALITY file is not provided
	""")
parser.add_argument(
	'-phred', metavar = 'PHRED', type = int, default = 33, help = """
	specify the PHRED value for the FASTQ file, 33 or 64
	""")
	
args = parser.parse_args()

header = ""
sequence = ""
quality = ""
qf = None

out = open(args.f.name.replace(".fna", "").replace(".fasta", "") + '.fastq', 'w+')

def parseSequence(header, sequence, quality, args, out):
	if len(sequence) == 0:
		return False
	output_quality = ""
	if len(quality) == 0:
		for i in range(len(sequence)):
			output_quality += chr(args.dq + args.phred)
	else:
		col = quality.split(" ")
		for i in col:
			if i is not "":
				output_quality += chr(int(i) + args.phred)
	if out is not None:
		out.write("@%s\n%s\n+\n%s\n" % (header, sequence, output_quality))
	return True

firstline = True
for f in args.f:
	f = f.strip()
	if firstline:
		firstline = False
		if len(f) == 0 or f[0] != ">":
			sys.stderr.write("FASTA file is corrupt or wrong file used. First line has to start with '>'\n")
			sys.exit(1)
	if args.qf is not None:
		qf = args.qf.readline().strip()
	if len(f) > 0 and f[0] == ">":
		if qf is not None and f.split(" ")[0] != qf.split(" ")[0]:
			raise ValueError('Fasta and quality file headers do not match')
		parseSequence(header, sequence, quality, args, out)
		header = f[1:]
		sequence = ""
		quality = ""
	else:
		sequence += f
		if len(sequence) >= 1000000:
			sys.stderr.write("Found sequence over 1 million bases, script terminated to avoid running out of memory\n")
			sys.exit(1)
		if qf is not None:
			quality += "%s " % qf

parseSequence(header, sequence, quality, args, out)

if out is not None:
	out.close()
