#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import gzip
import array
import sys

"""
	Demultiplex 454 sequences by barcode and primer into samples, filtering out low quality 
	sequences using average quality threshold.
"""
parser = argparse.ArgumentParser(description = """ 
	Demultiplex 454 sequences by barcode and primer into samples, filtering out low quality 
	sequences using average quality threshold.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', type = open, help = """
	FASTA file of raw sequences, supports gz packed files
	""")
parser.add_argument(
	'-fq', metavar = 'FASTQ_FILE', type = open, help = """
	FASTQ file of raw sequences. This file contains sequences and qualities, supports gz packed files
	""")
parser.add_argument(
	'-qf', metavar = 'QUALITY_FILE', type = open, help = """
	QUALITY file where sequence qualities are stored for each nucleotide. QUALITY
	file headers should match FASTA file headers, supports gz packed files
	""")
parser.add_argument(
	'-b', metavar = 'BARCODE_FILE', required = True, type = open, help = """
	BARCODE file where sample ID, barcodes and primers are stored in tabular 
	file format
	""")
parser.add_argument(
	'-bs', metavar = 'SAMPLE_COLUMN', required = True, type = int, help = """
	sample column in BARCODE file (left-most column starts from 1)
	""")
parser.add_argument(
	'-bb', metavar = 'BARCODE_COLUMN', required = True, type = int, help = """
	barcode column in BARCODE file
	""")
parser.add_argument(
	'-bp', metavar = 'PRIMER_COLUMN', type = int, help = """
	primer column in BARCODE file
	""")
parser.add_argument(
	'-q', metavar = 'AVERAGE_QUALITY', type = float, default = 25, help = """
	specify the lower limit of average quality to exclude the sequence
	(%(default)s)
	""")
parser.add_argument(
	'-min_allowed_base', metavar = 'MIN_BASE_QUALITY', default = 0, type = int, help = """
	specify the minimum allowed base quality of the base to exclude the sequence  
	""")
parser.add_argument(
	'-min_base_trimmed', metavar = 'TRIM_BASE_QUALITY', default = 0, type = int, help = """
	specify the minimum allowed base quality of the base to retain the 3' strand (recommended = 10)
	""")	
parser.add_argument(
	'-trimq', metavar = 'TRIM_QUALITY', type = float, help = """
	specify the lower limit of the average quality of within the trimming window to retain the sequence end (recommended = 20)
	""")
parser.add_argument(
	'-trimw', metavar = 'TRIM_WINDOW', type = int, help = """
	specify the window size to calculate average quality to for trimming the sequence end (recommended = 50)
	""")
parser.add_argument(
	'-ml', metavar = 'MINIMUM_LENGTH', type = int, default = 50, help = """
	specify the minimum allowed sequence length of the sequence to include after trimming  (%(default)s)
	""")
parser.add_argument(
	'-tl', metavar = 'TRUNCATE_LENGTH', type = int, help = """
	truncate sequences longer than the provided length to remove low quality regions, reverse primers or adapter sequences
	""")
parser.add_argument(
	'-phred', metavar = 'PHRED', type = int, default = 33, help = """
	specify PHRED score for the FASTQ file: 33 or 64
	""")
parser.add_argument(
	'-fasta', metavar = 'FASTA_OUTPUT', type = bool, default = False, help = """
	write output in FASTA format
	""")
parser.add_argument(
	'-mismatch', metavar = 'ALLOW_MISMATCH', type = int, default = 0, help = """
	specify mismatches in barcode and/or primer:
		0 - no mismatch
		1 - mismatch in barcode
		2 - mismatch in primer
		3 - mismatch in both
	""")
parser.add_argument(
	'-allow_indel', metavar = 'ALLOW_INDEL', type = int, default = 0, help = """
	specify insertion or deletion errors in primer:
	0 - no indels
	1 - allow indels in primer
	""")
parser.add_argument(
	'-remove_primer', metavar = 'REMOVE_PRIMER', default = True, type = bool, help = """
	specify whether primer is removed from the sequence
	""")
#parser.add_argument(
#	'-allow_primer_position', metavar = "ALLOW_PRIMER_BASES", default = 1, type = int, help = """
#	specify primer starting position within the sequence
#	""")
parser.add_argument(
	'-primer', metavar = 'FORWARD_PRIMER_SEQUENCE', default = "", help = """
	specify forward primer; supports degenerated primers. This will override primer information in the sample sheet file
	""")
parser.add_argument(
	'-homopolymer', metavar = 'HOMOPOLYMER_BASES', default = 0, type = int, help = """
	truncate sequences with repeating bases by reducing homopolymer length to user specific value
	""")
parser.add_argument(
	'-ignore_header', metavar = 'IGNORE_HEADER', default = False, type = bool, help = """
	do not check that FASTA and QUALITY file headers match
	""")
parser.add_argument(
	'-adapter', metavar = 'ADAPTER_SEQUENCE', default = "", help = """
	specify overhang adapter sequences, which can occur when read length is longer than DNA insert size, use commas to separate multiple sequences
	""")
	
args = parser.parse_args()

def calculateAverageQuality(quality, fq, phred):
	if len(quality) > 0:
		if fq is not None:
			if sys.version_info >= (3, 0):
				arr = array.array('b')
				arr.frombytes(quality.encode())
			else:
				arr = array.array('b', quality.encode())
			return (sum(arr) / float(len(quality))) - phred
		else:
			return sum(quality) / len(quality)
	else:
		return 0

def calculateSlidingWindow(quality, threshold, window, fq, phred):
	for i in range(len(quality) - window):
		if fq is not None:
			if sys.version_info >= (3, 0):
				arr = array.array('b')
				arr.frombytes(quality[i:i + window].encode())
			else:
				arr = array.array('b', quality[i:i + window].encode())
			avg = sum(arr) / float(window) - phred
			if avg < threshold:
				return i
		else:
			avg = sum(quality[i:(i + window)]) / float(window)
			if avg < threshold:
				return i
	return len(quality)

def parseSequence(header, sequence, quality):
	global args, lookup, output, homopolymers
	for k in lookup_lens:
		_length = len(sequence)
		key = sequence[:k]
		if key in lookup:
			_tmp = lookup[key]
			# trim based on minimum base criteria
			if args.min_base_trimmed > 0:
				for i in range(len(quality)):
					if quality[i] < args.min_base_trimmed:
						sequence = sequence[:i - 1]
						quality = quality[:i - 1]
						break
			# include based on minimum base criteria
			if args.min_allowed_base > 0:
				for i in quality:
					if i < args.min_allowed_base:
						return False
			if len(homopolymers) > 0:
				for homopolymer in homopolymers:
					if homopolymer in sequence:
						return False
			if args.remove_primer is False:
				cutoff = len(_tmp[1])
			else:
				cutoff = len(key)
			sequence = sequence[cutoff:]

			if adapters:
				for adapter in adapters:
					pos = sequence.find(adapter)
					if pos > 0:
						sequence = sequence[:pos]
						break
			if args.fq is not None:
				quality = quality[cutoff:]
			if args.tl is not None:
				sequence = sequence[0:args.tl]
				if args.fq is not None:
					quality = quality[0:args.tl]
			avg_quality = 40
			if len(quality) > 0:
				if args.trimq is not None and args.trimw is not None:
					trim_len = calculateSlidingWindow(quality, args.trimq, args.trimw, args.fq, args.phred)
					if trim_len < len(quality):
						quality = quality[:trim_len]
						sequence = sequence[:trim_len]
				avg_quality = calculateAverageQuality(quality, args.fq, args.phred)
			if avg_quality >= args.q and len(sequence) >= args.ml:
				if args.fq is not None and args.fasta is False:
					output.write("@%s|%s\t%s\t%s\t%s\t%.2f\t%d\t%d\n%s\n+\n%s\n" % (_tmp[0], header[1:], _tmp[0], _tmp[1], _tmp[2], avg_quality, _length, len(sequence), sequence, quality))
				else:
					output.write(">%s|%s\t%s\t%s\t%s\t%.2f\t%d\t%d\n%s\n" % (_tmp[0], header[1:], _tmp[0], _tmp[1], _tmp[2], avg_quality, _length, len(sequence), sequence))
				return True
	return False

def permutations(primer, conversion):
	permutations = {}
	permutations_total = 1
	arr = [primer]

	for i in range(len(primer)):
		if primer[i] in conversion:
			permutations[i] = conversion[primer[i]]
			permutations_total *= len(conversion[primer[i]])

	if(permutations_total > 10000):
		sys.stderr.write("Barcode/primer '%s' combination contains too many IUPAC code permutations and will be ignored\n" % primer)
		return [primer]

	for i in range(permutations_total):
		base = 0
		newseq = primer
		for pos in permutations:
			if base == 0:
				base_i = i % len(permutations[pos])
			else:
				base_i = i // base % len(permutations[pos])
			newseq = newseq[:pos] + permutations[pos][base_i] + newseq[pos + 1:]
			arr.append(newseq)
			base += len(permutations[pos])
	return arr

def parseBarcodes(line):
	global lookup, lookup_lens, args, maxcols
	col = line.split("\t")
	_sample = col[args.bs - 1].replace(" ", "-")
	_barcode = col[args.bb - 1]
	_primer = ""
	if args.bp is not None:
		_primer = col[args.bp - 1]
	if args.primer is not None:
		_primer = args.primer
	_primers = _primer.upper().split(",")

	_tmp = "%s%s" % (_barcode, _primer)

	# check nucleotide content of the sequence, if wrong chars, ignore
	_miss = 0
	for i in _tmp:
		if i not in nucleotides and i not in conversion:
			_miss += 1
	if _miss > 0:
		sys.stderr.write("Unknown nucleotide(s) used for sample '%s', sample entry ignored\n" % _sample)
		return

	# allow insertion and deletions errors
	if args.allow_indel:
		tmp = {}
		for primer in _primers:
			for i in range(len(primer)):
				# insertions
				tmp["%sN%s" % (primer[0:i], primer[i + 1:])] = True
				# deletions
				if i < len(primer) - 1:
					tmp["%s%sN" % (primer[0:i], primer[i + 2:])] = True
		for primer in tmp:
			_primers.append(primer)

	for _p in _primers:
		index = _barcode + _p.strip()

		if args.mismatch == 0:
			lookup[index] = [_sample, _barcode, _primer]
			lookup_lens[len(index)] = True
		else:
			tmp_list = {}
			for i in range(len(index)):
				if (i < len(_barcode) and (args.mismatch == 1 or args.mismatch == 3)) or (i >= len(_barcode) and (args.mismatch == 2 or args.mismatch == 3)):
					for nucleotide in nucleotides:
						tmp_list["%s%s%s" % (index[0:i], nucleotide, index[i + 1:])] = True
			for k in tmp_list:
				lookup[k] = [_sample, _barcode, _primer]
				lookup_lens[len(k)] = True

def packedFile(filename):
	ext = filename.split(".")[-1]
	if ext.lower() == "gz":
		return True
	else:
		return False
"""
	Read BARCODES into memory
"""
lookup = {}
lookup_lens = {}
adapters = []
nucleotides = ['A', 'C', 'G', 'T']
#IUPAC table
conversion = {"R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"], "K": ["G", "T"],
	"M": ["A", "C"], "B": ["C", "G", "T"], "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
	"N": ["A", "C", "G", "T"]}
homopolymers = []
if args.homopolymer > 0:
	for nucleotide in nucleotides:
		homopolymers.append(nucleotide * args.homopolymer)

maxcols = max([args.bs, args.bb, args.bp])

if packedFile(args.b.name):
	with gzip.open(args.b.name, 'r') as f:
		for line in f:
			parseBarcodes(line.decode())
else:
	with open(args.b.name, 'r') as f:
		for line in f:
			parseBarcodes(line)

# check all the primers and use IUPAC table to convert new set of primers
tmp_dict = {}
for primer in lookup:
	tmp = permutations(primer, conversion)
	if len(tmp) > 0:
		for t in tmp:
			tmp_dict[t] = lookup[primer]
for primer in tmp_dict:
	lookup[primer] = tmp_dict[primer]
tmp_dict = {}

if args.adapter:
	col = args.adapter.split(",")
	for c in col:
		adapters.append(c.strip())

"""
	Read FASTA and QUALITY file simultaneously and write output
"""
if args.f is not None:
	output = open(args.f.name.replace(".fna", "").replace(".fasta", "").replace(".gz", "").replace(".tar", "") + '.cleaned.fasta', 'w+')	# write results into file
elif args.fq is not None:
	output_extension = "fasta" if args.fasta is True else "fastq"
	output = open(args.fq.name.replace(".fq", "").replace(".fastq", "").replace(".gz", "").replace(".tar", "") + '.cleaned.%s' % (output_extension), 'w+')	# write results into file
	args.f = args.fq
else:
	raise ValueError("Argument -f or -fq is required")

header = ""
sequence = ""
quality = []
seq_i = 0
seq_c = 0
k = 0
qf = None
packed_files = [False, False]
# open files depending if they are packed or not
if packedFile(args.f.name):
	f = gzip.open(args.f.name, 'r')
	packed_files[0] = True
else:
	f = open(args.f.name, 'r')
if args.qf is not None:
	if packedFile(args.qf.name):
		qf = gzip.open(args.qf.name, 'r')
		packed_files[1] = True
	else:
		qf = open(args.qf.name, 'r')
with f as fh:
	for r1 in fh:
		r1 = r1.strip()
		if packed_files[0] and sys.version_info >= (3, 0):
			r1 = r1.decode()
		r2 = None
		# FASTQ reading
		if args.fq is not None:
			if k % 4 == 0:
				seq_i += 1
				if seq_i % 100000 == 0:
					console.log("%d/%d sequences cleaned/parsed\n" % (seq_c, seq_i))
				header = r1.split(" ")[0].split("\t")[0]
			elif k % 4 == 1:
				sequence = r1
			elif k % 4 == 3:
				quality = r1
				if parseSequence(header, sequence, quality):
					seq_c += 1
			k += 1
		else:
			if qf is not None:
				r2 = qf.readline().strip()
				if packed_files[1] and sys.version_info >= (3, 0):
					r2 = r2.decode()
			if r1.startswith(">"):
				seq_i += 1
				if seq_i % 100000 == 0:
					console.log("%d/%d sequences cleaned/parsed\n" % (seq_c, seq_i))
					
				if not args.ignore_header and r2 is not None and r1.split(" ")[0] != r2.split(" ")[0]:
					raise ValueError('Fasta and quality file headers do not match')
				if parseSequence(header, sequence, quality):
					seq_c += 1
				header = r1.split(" ")[0].split("\t")[0]
				sequence = ""
				quality = []
			else:
				sequence += r1
				if r2 is not None:
					# add together average quality
					col = r2.split(" ")
					for i in col:
						if i is not '':
							quality.append(int(i))

# add final sequence
if args.f is not None and parseSequence(header, sequence, quality):
	seq_c += 1
console.log("%d/%d sequences cleaned/parsed\n" % (seq_c, seq_i))
if f is not None:
	f.close()
