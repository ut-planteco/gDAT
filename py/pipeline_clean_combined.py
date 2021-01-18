from __future__ import division

import os
import argparse
import re
import sys
import itertools
import gzip
import array

"""
	Clean and filter demultiplexed sequences (combined into single file) based on primer and quality thresholds.
"""
parser = argparse.ArgumentParser(description = """
	Clean and filter demultiplexed sequences (combined into single file) based on primer and quality thresholds.
	""")
parser.add_argument(
	'-f', metavar = 'COMBINED_READ', required = True, type = open, help = """
	specify a FASTQ file, supports gz packed files
	""")
parser.add_argument(
	'-q', metavar = 'AVERAGE_QUALITY', type = float, default = 30, help = """
	specify lower limit of average quality to retain the sequence
	(%(default)s)
	""")
parser.add_argument(
	'-min_allowed_base', metavar = 'MIN_BASE_QUALITY', default = 0, type = int, help = """
	specify minimum allowed base quality to include the sequence
	""")
parser.add_argument(
	'-trimq', metavar = 'TRIM_QUALITY', type = float, help = """
	specify minimum average quality in the trimming window to trim the sequence end (recommended = 20)
	""")
parser.add_argument(
	'-trimw', metavar = 'TRIM_WINDOW', type = int, help = """
	specify window size to calculate average quality for sequence trimming (recommended = 50)
	""")
parser.add_argument(
	'-ml', metavar = 'MINIMUM_LENGTH', type = int, default = 50, help = """
	specify minimum allowed sequence length of the sequence to include after trimming (%(default)s)
	""")
parser.add_argument(
	'-phred', metavar = 'PHRED', type = int, default = 33, help = """
	specify the PHRED score for the FASTQ file, 33 or 64
	""")
parser.add_argument(
	'-fasta', metavar = 'OUTPUT_FASTA', type = bool, default = False, help = """
	write output in FASTA format
	""")
parser.add_argument(
	'-mismatch', metavar = 'ALLOW_MISMATCH', type = int, default = 0, help = """
	specify whether mismatches in the primer are allowed:
		0 - no mismatch
		1 - allow mismatch
	""")
parser.add_argument(
	'-allow_indel', metavar = 'ALLOW_INDEL', type = int, default = 0, help = """
	specify whether insertion or deletion errors in primer are allowed:
	0 - no indels
	1 - allow indels in primer
	""")
parser.add_argument(
	'-remove_primer', metavar = 'REMOVE_PRIMER', default = True, type = bool, help = """
	specify whether the primer is removed from the sequence
	""")
parser.add_argument(
	'-forward_primer', metavar = 'FORWARD_PRIMER_SEQUENCE', default = "", help = """
	specify forward primers, use commas to separate multiple primers, supports degenerated primers
	""")
parser.add_argument(
	'-reverse_primer', metavar = 'REVERSE_PRIMER_SEQUENCE', default = "", help = """
	specify reverse primers, use commas to separate multiple primers, supports degenerated primers
	""")
#parser.add_argument(
#	'-allow_primer_position', metavar = "ALLOW_PRIMER_BASES", default = 1, type = int, help = """
#	specify allowed primer starting position within the sequence
#	""")
parser.add_argument(
	'-primers_mixed', metavar = 'PRIMERS_MIXED', default = False, type = bool, help = """
	specify whether primers are mixed between forward and reverse reads, allows to check both primers and outputs sequences in correct order for pairing
	""")
parser.add_argument(
	'-homopolymer', metavar = 'HOMOPOLYMER_BASES', default = 0, type = int, help = """
	truncate sequences with repeating bases by reducing homopolymer length to user specific value
	""")
	
args = parser.parse_args()

nucleotides = ['A', 'C', 'G', 'T']
extras = [' ', ',']
#IUPAC table
conversion = {"R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"], "K": ["G", "T"],
	"M": ["A", "C"], "B": ["C", "G", "T"], "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
	"N": ["A", "C", "G", "T"]}
homopolymers = []
if args.homopolymer > 0:
	for nucleotide in nucleotides:
		homopolymers.append(nucleotide * args.homopolymer)

def avgQuality(qual, phred):
	if not qual:
		return 0
	if sys.version_info >= (3, 0):
		arr = array.array('b')
		try:
			arr.frombytes(qual.encode())
		except:
			arr.frombytes(qual)
	else:
		try:
			arr = array.array('b', qual.encode())
		except:
			arr = array.array('b', qual)
	return (sum(arr) / float(len(qual))) - phred

def calculateSlidingWindow(quality, threshold, window, phred):
	for i in range(len(quality) - window):
		if sys.version_info >= (3, 0):
			arr = array.array('b')
			try:
				arr.frombytes(quality[i:i + window].encode())
			except:
				arr.frombytes(quality[i:i + window])
		else:
			try:
				arr = array.array('b', quality[i:i + window].encode())
			except:
				arr = array.array('b', quality[i:i + window])
		avg = sum(arr) / float(window) - phred
		if avg < threshold:
			return i
	return len(quality)

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

def rev_comp(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(nn[n] for n in reversed(st))

def packedFile(filename):
	ext = filename.split(".")[-1]
	if ext.lower() == "gz":
		return True
	else:
		return False

out = open("%s.cleaned.%s" % (args.f.name.replace(".fastq", "").replace(".fq", ""), "fasta" if args.fasta is True else "fastq"), "w+")
header = ""
sequence = ""
sequence1 = ""
sequence2 = ""
quality = ""
seq_i = 0
seq_sel = 0
files_f = 0
files_r = 0
flookup_lens = {}
rlookup_lens = {}

flookup = {}
rlookup = {}
fprimer = args.forward_primer.upper()
rprimer = args.reverse_primer.upper()

fprimers = fprimer.split(",")
rprimers = rprimer.split(",")
# strip empty spaces
for k, v in enumerate(fprimers):
	fprimers[k] = v.strip()
for k, v in enumerate(rprimers):
	rprimers[k] = v.strip()

_tmp = "%s%s" % (fprimer, rprimer)

# check nucleotide content of the sequence, if wrong chars, ignore
_miss = 0
for i in _tmp:
	if i not in nucleotides and i not in conversion and i not in extras:
		_miss += 1
if _miss > 0:
	sys.stderr.write("Unknown nucleotide(s) used for primer\n" % sample)
	sys.exit(1)		

# allow insertion and deletions errors
if args.allow_indel:
	tmp = {}
	for fprimer in fprimers:
		for i in range(len(fprimer)):
			# insertions
			tmp["%sN%s" % (fprimer[0:i], fprimer[i + 1:])] = True
			# deletions
			if i < len(fprimer) - 1:
				tmp["%s%sN" % (fprimer[0:i], fprimer[i + 2:])] = True
	for primer in tmp:
		fprimers.append(primer)
	tmp = {}
	for rprimer in rprimers:
		for i in range(len(rprimer)):
			# insertions
			tmp["%sN%s" % (rprimer[0:i], rprimer[i + 1:])] = True
			# deletions
			if i < len(rprimer) - 1:
				tmp["%s%sN" % (rprimer[0:i], rprimer[i + 2:])] = True
	for primer in tmp:
		rprimers.append(primer)	

if args.mismatch == 0:
	for primer in fprimers:
		flookup_lens[len(primer)] = True
		flookup[primer] = True
	for primer in rprimers:
		rlookup_lens[len(primer)] = True
		rlookup[primer] = True
else:
	for primer in fprimers:
		flookup_lens[len(primer)] = True
		for i in range(len(primer)):
			for nucleotide in nucleotides:
				flookup["%s%s%s" % (primer[0:i], nucleotide, primer[i + 1:])] = True
	for primer in rprimers:
		rlookup_lens[len(primer)] = True
		for i in range(len(primer)):
			for nucleotide in nucleotides:
				rlookup["%s%s%s" % (primer[0:i], nucleotide, primer[i + 1:])] = True
# check all the primers and use IUPAC table to convert new set of primers
tmp_dict = {}
for primer in flookup:
	tmp = permutations(primer, conversion)
	if len(tmp) > 0:
		for t in tmp:
			tmp_dict[t] = True
for primer in tmp_dict:
	flookup[primer] = True
tmp_dict = {}
for primer in rlookup:
	tmp = permutations(primer, conversion)
	if len(tmp) > 0:
		for t in tmp:
			tmp_dict[t] = True
for primer in tmp_dict:
	rlookup[primer] = True
tmp_dict = {}

i = 0
file1_packed = False
if packedFile(args.f.name):
	fh = gzip.open(args.f.name, "r")
	file1_packed = True
else:
	fh = open(args.f.name, "r")

for line in fh:
	line = line.strip()
	if file1_packed:
		line = line.decode()
	if i % 4 == 0:
		seq_i += 1
		if seq_i % 100000 == 0:
			sys.stderr.write("%d/%d sequences cleaned/parsed\n" % (seq_sel, seq_i, files_f, files_r))
		header = line[1:].split(" ")[0]
	elif i % 4 == 1:
		sequence = line
	elif i % 4 == 3:
		if args.trimq is not None and args.trimw is not None:
			trim_len = calculateSlidingWindow(line, args.trimq, args.trimw, args.phred)
			if trim_len < len(line):
				line = line[:trim_len]
				sequence = sequence[:trim_len]
		fprimer_len = 0
		selected = False
		selection = 0
		selectiont = 0
		correct_strand = False
		if len(fprimer) > 0:
			selectiont += 1
			for l in flookup_lens:
				if sequence[:l] in flookup:
					selected = True
					selection += 1
					fprimer_len = l
					correct_strand = True
					break
			# no primer found, try reverse complement?
			if not correct_strand and args.primers_mixed:
				sequence = rev_comp(sequence)
				line = reversed(line)
				for l in flookup_lens:
					if sequence[:l] in flookup:
						selected = True
						selection += 1
						fprimer_len = l
						correct_strand = False
						break
		fprimer_len = 0
		rprimer_len = 0
		if len(rprimer) > 0:
			selectiont += 1
			tmp = rev_comp(sequence)
			for l in rlookup_lens:
				if tmp[:l] in rlookup:
					selection += 1
					rprimer_len = l
					break
		if selectiont == selection:
			selected = True
		if selected and args.remove_primer:
			if correct_strand:
				if fprimer_len > 0:
					sequence = sequence[fprimer_len:]
					line = line[fprimer_len:]
				if rprimer_len > 0:
					sequence = sequence[:-rprimer_len]
					line = line[:-rprimer_len]
			else:
				if fprimer_len > 0:
					sequence = sequence[rprimer_len:]
					line = line[rprimer_len:]
				if rprimer_len > 0:
					sequence = sequence[:-fprimer_len]
					line = line[:-fprimer_len]
		# include based on minimum base criteria
		if selected and args.min_allowed_base > 0:
			for _i in line:
				if ord(_i) - args.phred < args.min_allowed_base:
					selected = False
					break
		if selected and avgQuality(line, args.phred) < args.q:
			selected = False
		if selected and args.ml is not None and len(sequence) < args.ml:
			selected = False
		if selected and len(homopolymers) > 0:
			for homopolymer in homopolymers:
				if homopolymer in sequence:
					# truncate homopolymer reads and quality
					_s = ""
					_q = ""
					_chr = ""
					_chr_i = 0
					_pos = -1
					for _ in sequence:
						_pos += 1
						if _chr == _:
							_chr_i += 1
						else:
							_chr_i = 0
							_chr = _
						if _chr_i >= args.homopolymer:
							continue
						_s += _
						if not args.fasta:
							_q += line[_pos]
					sequence = _s
					if not args.fasta:
						line = _q
					break
		if selected:
			seq_sel += 1
			if args.fasta:
				out.write(">%s\n%s\n" % (header, sequence))
			else:
				out.write("@%s\n%s\n+\n%s\n" % (header, sequence, line))
	i += 1

sys.stderr.write("%d/%d sequences cleaned/parsed\n" % (seq_sel, seq_i))

if out:
	out.close()
