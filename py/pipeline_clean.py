#!/usr/bin/env python
from __future__ import division

import os
import argparse
import gzip
import console
import sys
import array

"""
	Cleaning of Illumina raw single-end or paired-end reads by barcode and primer checking and 
	quality filtering.
"""
parser = argparse.ArgumentParser(description = """
	Cleaning of Illumina raw single-end or paired-end reads by barcode and primer checking and 
	quality filtering.
	""")
parser.add_argument(
	'-fr', metavar = 'FORWARD_READ', required = True, type = open, help = """
	specify a FASTQ file containing forward raw reads (R1), supports gz packed files
	""")
parser.add_argument(
	'-rr', metavar = 'REVERSE_READ', type = open, help = """
	specify a FASTQ file containing reverse raw reads (R2 or R4), supports gz packed files
	""")
parser.add_argument(
	'-fo', metavar = 'FORWARD_OLIGO', type = open, help = """
	specify a FASTQ file containing forward oligo/index reads (barcodes), (R2), supports gz packed files
	""")
parser.add_argument(
	'-ro', metavar = 'REVERSE_OLIGO', type = open, help = """
	specify a FASTQ file containing reverse oligo/index reads (barcodes), (R3), supports gz packed files
	""")
parser.add_argument(
	'-b', metavar = 'BARCODE_FILE', required = True, type = open, help = """
	specify sample sheet file, where sample ID, barcodes and primers are stored in tabular format
	""")
parser.add_argument(
	'-bs', metavar = 'SAMPLE_COLUMN', required = True, type = int, help = """
	specify sample column in the BARCODE file (the left-most column is numbered 1)
	""")
parser.add_argument(
	'-bfb', metavar = 'FORWARD_BARCODE_COLUMN', required = True, type = int, help = """
	specify forward barcode column number in the BARCODE file
	""")
parser.add_argument(
	'-bfp', metavar = 'FORWARD_PRIMER_COLUMN', type = int, help = """
	specify forward primer column number in the BARCODE file
	""")
parser.add_argument(
	'-brb', metavar = 'REVERSE_BARCODE_COLUMN', type = int, help = """
	specify reverse barcode column number in the BARCODE file (if forward and reverse barcodes are concatenated, use same number as for the forward barcode). When using only forward reads and also marking down reverse barcode, reverse barcode will be generated as reverse complement and will be searched at the end of the sequence
	""")
parser.add_argument(
	'-brp', metavar = 'REVERSE_PRIMER_COLUMN', type = int, help = """
	specify reverse primer column number in BARCODE file
	""")
parser.add_argument(
	'-quality', metavar = 'AVERAGE_QUALITY', default = 30, type = int, help = """
	specify the lower limit of average quality required to retain the sequence
	""")
parser.add_argument(
	'-min_allowed_base', metavar = 'MIN_BASE_QUALITY', default = 0, type = int, help = """
	specify the minimum allowed base quality required to retain the sequence
	""")
parser.add_argument(
	'-min_base_trimmed', metavar = 'TRIM_BASE_QUAILTY', default = 0, type = int, help = """
	specify the minimum allowed base quality below which the to trim 3' strand should be trimmed (recommended = 10)
	""")
parser.add_argument(
	'-trimq', metavar = 'TRIM_QUALITY', type = float, help = """
	specify the minimum average quality of the trimming window below which the sequence end should be trimmed (recommended = 20)
	""")
parser.add_argument(
	'-trimw', metavar = 'TRIM_WINDOW', type = int, help = """
	specify the window size used to calculate average quality for sequence end trimming (recommended = 50)
	""")
parser.add_argument(
	'-ml', metavar = 'MINIMUM_LENGTH', type = int, default = 50, help = """
	specify the minimum sequence length to include after trimming (%(default)s)
	""")
parser.add_argument(
	'-mismatch', metavar = 'ALLOW_MISMATCH', type = int, default = 0, help = """
	specify whether mismates in the barcode and primer are allowed:
	0 - no mismatch
	1 - mismatch in barcode and primer
	""")
parser.add_argument(
	'-allow_indel', metavar = 'ALLOW_INDEL', type = int, default = 0, help = """
	specify whether insertion or deletion errors in primer are allowed:
	0 - no indels
	1 - allow indels in primer
	""")
parser.add_argument(
	'-ignore_bases', metavar = 'IGNORE_BASES', default = 0, type = int, help = """
	ignores the number of first bases for quality and primer checking, can be used when sequencer has calibration issues to allow inclusion of low quality sequences.
	""")
parser.add_argument(
	'-phred', metavar = 'PHRED', type = int, default = 33, help = """
	specify PHRED score for the FASTQ file, 33 or 64
	""")
parser.add_argument(
	'-forward_trim', metavar = 'FORWARD_TRIM_BASES', default = 0, type = int, help = """
	specify the position at which to trim forward reads; use 220 for ITS2
	""")
parser.add_argument(
	'-reverse_trim', metavar = 'REVERSE_TRIM_BASES', default = 0, type = int, help = """
	specify the position at which to trim reverse reads; use 180 for ITS2
	""")
parser.add_argument(
	'-forward_primer', metavar = 'FORWARD_PRIMER_SEQUENCE', default = "", help = """
	specify the forward primers, use commas to separate multiple primers, supports degenerated primers
	""")
parser.add_argument(
	'-reverse_primer', metavar = 'REVERSE_PRIMER_SEQUENCE', default = "", help = """
	specify the reverse primers, use commas to separate multiple primers, supports degenerated primers
	""")
parser.add_argument(
	'-remove_primer', metavar = 'REMOVE_PRIMER', default = True, type = bool, help = """
	specify if primer is removed from the sequence
	""")
parser.add_argument(
	'-primers_mixed', metavar = 'PRIMERS_MIXED', default = False, type = bool, help = """
	specify if primers are mixed between forward and reverse reads, allows to check both primers and outputs sequences in correct order for pairing
	""")
#parser.add_argument(
#	'-allow_primer_position', metavar = "ALLOW_PRIMER_BASES", default = 1, type = int, help = """
#	specify allowed primer starting position within the sequence
#	""")
parser.add_argument(
	'-fasta', metavar = 'FASTA_OUTPUT', type = bool, default = False, help = """
	write output in FASTA format
	""")
parser.add_argument(
	'-homopolymer', metavar = 'HOMOPOLYMER_BASES', default = 0, type = int, help = """
	truncate sequences with repeating bases by reducing homopolymer length to user specific value
	""")
parser.add_argument(
	'-adapter', metavar = 'ADAPTER_SEQUENCE', default = "", help = """
	specify whether overhang adapter sequences may be present, which can occur when read length is longer than DNA insert size, use commas to separate multiple sequences
	""")
	
args = parser.parse_args()

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

def startsWith(full, sub):
	if len(full) >= len(sub) and sub == full[0:len(sub)]:
		return True
	else:
		return False

def permutations(primer, conversion):
	permutations = {}
	permutations_total = 1
	arr = [primer]

	for i in range(len(primer)):
		if primer[i] in conversion:
			permutations[i] = conversion[primer[i]]
			permutations_total *= len(conversion[primer[i]])

	if(permutations_total > 100000):
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

def reverseComplement(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
		'Y': 'R', 'R': 'Y', 'S': 'S', 'W': 'W', 'K': 'M',
		'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
		'N': 'N'}
	return "".join(complement.get(base, base) for base in reversed(seq))

def packedFile(filename):
	ext = filename.split(".")[-1]
	if ext.lower() == "gz":
		return True
	else:
		return False	
# allowed nucleotide swaps
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

# build lookup for barcodes as hashmap to speed up fetching the correct sample
# if mismatch allowed, build additional hashes with 1 nucleotide swap
lookup = {}
adapters = []
flookup_lens = {}
rlookup_lens = {}
barcode_lens = {}
barcode_lookup = {}
fprimer_list = {}
rprimer_list = {}
sys.stderr.write("Reading barcodes from samplesheet\n")
#TODO: if same barcode uses different primers, make a primer list to go through to see the CORRECT one
if packedFile(args.b.name):
	fh = gzip.open(args.b.name, "r")
else:
	fh = open(args.b.name, "r")

def castint(str):
	try:
		return int(str)
	except:
		return 0

if args.adapter:
	col = args.adapter.split(",")
	for c in col:
		adapters.append(c.strip())

#TODO: barcode and primer lookup should be split from each other, so it would be using smaller memory print
# maximum column
maxcols = max([castint(args.bs), castint(args.bfb), castint(args.brb), castint(args.bfp), castint(args.brp)])

for r in fh:
	if len(r) == 0 or (len(r) > 0 and r[0] == "#"):
		continue
	# if CSV file is used
	if "csv" in args.b.name.lower():
		# check if , or ; is used
		col = r.strip().replace("\"", "").split(",")
		if len(col) == 1:
			col = re.strip().replace("\"", "").split(";")
	else:
		# otherwise use tabs
		col = r.strip().split("\t")
	if len(col) < maxcols - 1:
		continue
	# check if correct columns are used
	_sample = col[args.bs - 1].replace(" ", "-")
	_fbarcode = col[args.bfb - 1].upper()
	_rbarcode = ""
	_fprimer = ""
	_rprimer = ""
	if args.brb:
		_rbarcode = col[args.brb - 1].upper()
		if args.bfb == args.brb:
			# split the barcode
			tmp = _fbarcode
			_fbarcode = tmp[:int(len(tmp) / 2)]
			_rbarcode = tmp[int(len(tmp) / 2):]
	if args.bfp:
		_fprimer = col[args.bfp - 1].upper()
	if args.brp:
		_rprimer = col[args.brp - 1].upper()
	# user input will override primer information
	if args.forward_primer:
		_fprimer = args.forward_primer.upper()
	if args.reverse_primer:
		_rprimer = args.reverse_primer.upper()
	
	if not args.rr:
		if _rbarcode:
			_rbarcode = reverseComplement(_rbarcode)
		if _rprimer:
			_rprimer = reverseComplement(_rprimer)

	_tmp = "%s%s%s%s" % (_fbarcode, _rbarcode, _fprimer, _rprimer)

	# check nucleotide content of the sequence, if wrong chars, ignore
	_miss = 0
	for i in _tmp:
		if i not in nucleotides and i not in conversion and i not in extras:
			_miss += 1
	if _miss > 0:
		sys.stderr.write("Unknown nucleotide(s) used for sample '%s', sample entry ignored\n" % _sample)
		#continue
	
	fprimers = _fprimer.split(",")
	rprimers = _rprimer.split(",")
	# strip empty spaces
	for k, v in enumerate(fprimers):
		fprimers[k] = v.strip()
	for k, v in enumerate(rprimers):
		rprimers[k] = v.strip()
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
	"""
		_index = _fbarcode + _rbarcode
		barcode_lens[len(_index)] = True
		barcode_lookup[_index] = [_sample, _fbarcode, _fprimer, _rbarcode, _rprimer]

		for fprimer in fprimers:
			flookup_lens[len(fprimer)] = True
			if args.mismatch == 0:
				lookup[fprimer] = _fprimer
			else:
				for i in range(len(index)):
					for nucleotide in nucleotides:
						lookup["%s%s%s" % (index[0:i], nucleotide, index[i + 1:])] = [_sample, _fbarcode, _fprimer, _rbarcode, _rprimer]
	"""

	if len(fprimers) > 1 or len(rprimers) > 1:
		for fprimer in fprimers:
			for rprimer in rprimers:
				if args.ignore_bases:
					if len(fprimer) > args.ignore_bases:
						fprimer = fprimer[args.ignore_bases:]
					if len(rprimer) > args.ignore_bases:
						rprimer = rprimer[args.ignore_bases:]
				findex = _fbarcode + fprimer
				if not args.rr:
					rindex = rprimer + _rbarcode
				else:
					rindex = _rbarcode + rprimer
				index = findex + rindex

				if len(findex) > 0:
					flookup_lens[len(findex)] = True
				if len(rindex) > 0:
					rlookup_lens[len(rindex)] = True

				if args.mismatch == 0:
					lookup[index] = [_sample, _fbarcode, _fprimer, _rbarcode, _rprimer]
				else:
					for i in range(len(index)):
						for nucleotide in nucleotides:
							lookup["%s%s%s" % (index[0:i], nucleotide, index[i + 1:])] = [_sample, _fbarcode, _fprimer, _rbarcode, _rprimer]
	else:
		if args.ignore_bases:
			if len(_fprimer) > args.ignore_bases:
				_fprimer = _fprimer[args.ignore_bases:]
			if len(_rprimer) > args.ignore_bases:
				_rprimer = _rprimer[args.ignore_bases:]
		findex = _fbarcode + _fprimer
		if not args.rr:
			rindex = _rprimer + _rbarcode
		else:
			rindex = _rbarcode + _rprimer
		index = findex + rindex

		if len(findex) > 0:
			flookup_lens[len(findex)] = True
		if len(rindex) > 0:
			rlookup_lens[len(rindex)] = True

		if args.mismatch == 0:
			lookup[index] = [_sample, _fbarcode, _fprimer, _rbarcode, _rprimer]
		else:
			for i in range(len(index)):
				for nucleotide in nucleotides:
					lookup["%s%s%s" % (index[0:i], nucleotide, index[i + 1:])] = [_sample, _fbarcode, _fprimer, _rbarcode, _rprimer]

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

#for k in flookup:
#	print("%s:%s,%s,%s" % (k, flookup[k][0], flookup[k][1], flookup[k][2]))
sys.stderr.write("Built %d hashes for barcodes\n" % (len(lookup)))

# define some helper variables
seq_i = 0
seq_sel = 0
header = ""
fseq = ""
foli = ""
rseq = ""
roli = ""
fqual = ""
rqual = ""
r1 = None
r2 = None
r3 = None
r4 = None
packed_files = [False, False, False, False]

output_extension = "fasta" if args.fasta is True else "fastq"
out = open(args.fr.name.replace(".fq", "").replace(".fastq", "").replace(".gz", "").replace(".tar", "") + '.cleaned.%s' % (output_extension), 'w+')	# write results into file
l = 0
if packedFile(args.fr.name):
	frh = gzip.open(args.fr.name, "r")
	packed_files[0] = True
else:
	frh = open(args.fr.name, "r")
if args.rr is not None:
	if packedFile(args.rr.name):
		rrh = gzip.open(args.rr.name, "r")
		packed_files[3] = True
	else:
		rrh = open(args.rr.name, "r")
if args.fo is not None:
	if packedFile(args.fo.name):
		foh = gzip.open(args.fo.name, "r")
		packed_files[1] = True
	else:
		foh = open(args.fo.name, "r")
if args.ro is not None:
	if packedFile(args.ro.name):
		roh = gzip.open(args.ro.name, "r")
		packed_files[2] = True
	else:
		roh = open(args.ro.name, "r")


# start reading forward file (R1)
for r1 in frh:
	if packed_files[0] and sys.version_info >= (3, 0):
		r1 = r1.decode()
	r1 = r1.strip()
	# try to read other files if they are presented, R4 is reverse here, R2 and R3 are oligos
	if args.rr is not None:
		r4 = rrh.readline().strip()
		if packed_files[3] and sys.version_info >= (3, 0):
			r4 = r4.decode()
	if args.fo is not None:
		r2 = foh.readline().strip()
		if packed_files[1] and sys.version_info >= (3, 0):
			r2 = r2.decode()
	if args.ro is not None:
		r3 = roh.readline().strip()
		if packed_files[2] and sys.version_info >= (3, 0):
			r3 = r3.decode()
	if l % 4 == 0:
		header = r1[1:].split(" ")[0]
		seq_i += 1
		if seq_i % 100000 == 0:
			sys.stderr.write("%d/%d sequences cleaned/parsed\n" % (seq_sel, seq_i))
	elif l % 4 == 1:
		fseq = r1
		if r2 is not None:
			foli = r2
		if r3 is not None:
			roli = r3
		if r4 is not None:
			rseq = r4
		if args.ignore_bases:
			if len(fseq) > args.ignore_bases:
				fseq = fseq[args.ignore_bases:]
			if len(rseq) > args.ignore_bases:
				rseq = rseq[args.ignore_bases:]
	elif l % 4 == 3:
		fqual = r1
		if r4 is not None:
			rqual = r4
		if args.ignore_bases:
			if len(fqual) > args.ignore_bases:
				fqual = fqual[args.ignore_bases:]
			if len(rqual) > args.ignore_bases:
				rqual = rqual[args.ignore_bases:]
		sample = None
		if r2:
			fseq = foli + fseq
			fqual = r2 + fqual
		if r3 and rseq:
			rseq = roli + rseq
			rqual = r3 + rqual
		for len1 in flookup_lens:
			if not sample:
				if rlookup_lens:
					for len2 in rlookup_lens:
						if not sample:
							if not rseq:
								# only use forward read, use primers at the front and at the end
								index = fseq[:len1] + fseq[-len2:]
							else:
								# use both forward and reverse reads, use primers at the front
								index = fseq[:len1] + rseq[:len2]
							if index in lookup:
								sample = lookup[index]
								break
							if args.primers_mixed:
								index = rseq[:len1] + fseq[:len2]
								if index in lookup:
									sample = lookup[index]
									fseq, rseq = rseq, fseq
									fqual, rqual = rqual, fqual
									break
				else:
					index = fseq[:len1]
					if index in lookup:
						sample = lookup[index]
						break
					if args.primers_mixed and rseq:
						index = rseq[:len1]
						if index in lookup:
							sample = lookup[index]
							break
		if sample:
			avg_quality = 40
			select = True
			if args.trimq is not None and args.trimw is not None:
				if len(fqual) > 0:
					trim_len = calculateSlidingWindow(fqual, args.trimq, args.trimw, args.phred)
					if trim_len < len(fqual):
						fqual = fqual[:trim_len]
						fseq = fseq[:trim_len]
				if len(rqual) > 0:
					trim_len = calculateSlidingWindow(rqual, args.trimq, args.trimw, args.phred)
					if trim_len < len(rqual):
						rqual = rqual[:trim_len]
						rseq = rseq[:trim_len]
			if args.remove_primer:
				fseq = fseq[len(sample[1] + sample[2]):]
				fqual = fqual[len(sample[1] + sample[2]):]
				if rseq:
					rseq = rseq[len(sample[3] + sample[4]):]
					rqual = rqual[len(sample[3] + sample[4]):]
				elif rlookup_lens:
					fseq = fseq[:-len(sample[3] + sample[4])]
					fqual = fqual[:-len(sample[3] + sample[4])]
			if args.forward_trim > 0:
				fseq = fseq[:args.forward_trim]
				fqual = fqual[:args.forward_trim]
			if rseq and args.reverse_trim > 0:
				rseq = rseq[:args.reverse_trim]
				rqual = rqual[:args.reverse_trim]
			if adapters:
				found = False
				for adapter in adapters:
					pos = fseq.find(adapter)
					if pos > 0:
						fseq = fseq[:pos]
						fqual = fqual[:pos]
						found = True
					if rseq:
						pos = rseq.find(adapter)
						if pos > 0:
							rseq = rseq[:pos]
							rqual = rqual[:pos]
							found = True
					if found:
						break
			# trim based on minimum base criteria
			if args.min_base_trimmed > 0:
				for i in range(len(fqual)):
					if ord(fqual[i]) - args.phred < args.min_base_trimmed:
						fseq = fseq[:i - 1]
						fqual = fqual[:i - 1] 
						break
				if rseq:
					for i in range(len(rqual)):
						if ord(rqual[i]) - args.phred < args.min_base_trimmed:
							rseq = rseq[:i - 1]
							rqual = rqual[:i - 1]
							break
			# include based on minimum base criteria
			if args.min_allowed_base > 0:
				for _i in fqual:
					if ord(_i) - args.phred < args.min_allowed_base:
						select = False
						break
				if rseq and select:
					for _i in rqual:
						if ord(_i) - args.phred < args.min_allowed_base:
							select = False
							break
			# check if minimum length is fulfilled for forward read
			if args.ml > 0 and len(fseq) < args.ml:
				select = False
			# check if minimum length is fulfilled for reverse read
			if rseq and args.ml > 0 and len(rseq) < args.ml:
				select = False
			# check average quality of the remaining sequence for forward/reverse read
			if rqual:
				if (avgQuality(fqual, args.phred) + avgQuality(rqual, args.phred)) / 2 < args.quality:
					select = False
			else:
				if avgQuality(fqual, args.phred) < args.quality:
					select = False
			# check homopolymers and remove sequences if they occur in forward or reverse reads
			if len(homopolymers) > 0:
				for homopolymer in homopolymers:
					if homopolymer in fseq:
						# truncate homopolymer with quality
						_s = ""
						_q = ""
						_chr = ""
						_chr_i = 0
						_pos = -1
						for _ in fseq:
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
								_q += fqual[_pos]
						fseq = _s
						if not args.fasta:
							fqual = _q
						break
				for homopolymer in homopolymers:
					if rseq and homopolymer in rseq:
						# truncate homopolymer with quality
						_s = ""
						_q = ""
						_chr = ""
						_chr_i = 0
						_pos = -1
						for _ in rseq:
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
								_q += rqual[_pos]
						rseq = _s
						if not args.fasta:
							rqual = _q
						break
			if select:
				seq_sel += 1
				if r1 and r4:
					if args.fasta:
						out.write(">%s|%d|%s\n%s\n" % (sample[0], seq_sel, header, fseq))
						out.write(">%s|%d|%s\n%s\n" % (sample[0], seq_sel, header, rseq))
					else:
						out.write("@%s|%d|%s\n%s\n+\n%s\n" % (sample[0], seq_sel, header, fseq, fqual))
						out.write("@%s|%d|%s\n%s\n+\n%s\n" % (sample[0], seq_sel, header, rseq, rqual))
				else:
					if args.fasta:
						out.write(">%s|%d|%s\n%s\n" % (sample[0], seq_sel, header, fseq))
					else:
						out.write("@%s|%d|%s\n%s\n+\n%s\n" % (sample[0], seq_sel, header, fseq, fqual))

	l += 1
if out:
	out.close()

sys.stderr.write("%d/%d sequences cleaned/parsed\n" % (seq_sel, seq_i))
