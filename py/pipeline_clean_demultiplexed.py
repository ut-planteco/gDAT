from __future__ import division

import os
import argparse
import re
import sys
import itertools
import gzip
import array

"""
	Clean and filter demultiplexed sequences based on primer and quality thresholds.
"""
parser = argparse.ArgumentParser(description = """
	Clean and filter demultiplexed sequences based on primer and quality thresholds.
	""")
parser.add_argument(
	'-folder', metavar = 'FOLDER', type = str, required = True, help = """
	specify folder containing FASTQ files, supports gz packed files
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
	'-min_base_trimmed', metavar = 'TRIM_BASE_QUAILTY', default = 0, type = int, help = """
	specify minimum allowed base quality to trim the 3' strand (recommended = 10)
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
	'-forward_trim', metavar = 'FORWARD_TRIM_BASES', default = 0, type = int, help = """
	specify the position at which to trim forward reads; use 220 for ITS2
	""")
parser.add_argument(
	'-reverse_trim', metavar = 'REVERSE_TRIM_BASES', default = 0, type = int, help = """
	specify the position at which to trim reverse reads; use 180 for ITS2
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
	'-ignore_bases', metavar = 'IGNORE_BASES', default = 0, type = int, help = """
	ignores the number of first bases for quality and primer checking, can be used when sequencer has calibration issues to allow inclusion of low quality sequences.
	""")
parser.add_argument(
	'-remove_primer', metavar = 'REMOVE_PRIMER', default = True, type = bool, help = """
	specify whether the primer is removed from the sequence
	""")
parser.add_argument(
	'-subfolders', metavar = 'SUBFOLDERS', default = False, type = bool, help = """
	specify whether to search subfolders for input files
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
	'-force_forward', metavar = 'FORCE_FORWARD', default = False, help = """
	specify whether to use and output only forward reads
	""")
parser.add_argument(
	'-force_reverse', metavar = 'FORCE_REVERSE', default = False, help = """
	specify whether to use and output only reverse reads
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

def packedFile(filename):
	ext = filename.split(".")[-1]
	if ext.lower() == "gz":
		return True
	else:
		return False

files = []
adapters = []
tmp = {}
cnt = [0, 0, 0]
folders = []

if args.adapter:
	col = args.adapter.split(",")
	for c in col:
		adapters.append(c.strip())

if args.subfolders:
	folders = [x[0] for x in os.walk(args.folder)]
else:
	folders = [args.folder]
for folder in folders:
	for f in os.listdir(folder):
		ext = f.lower().split(".")
		if ext[-1] == "fq" or ext[-1] == "fastq" or (ext[-1] == "gz" and ("fq" in ext or "fastq" in ext)):
			if "_R1_" in f:
				cnt[0] += 1
				cnt[2] += 1
				sample = f.split("_")[0]
				if f in tmp:
					tmp[f][0] = "%s/%s" % (folder, f)
				else:
					tmp[f] = ["%s/%s" % (folder, f), "", sample]
			if "_R2_" in f:
				cnt[1] += 1
				cnt[2] += 1
				r1 = f.replace("_R2_", "_R1_")
				sample = f.split("_")[0]
				if r1 in tmp:
					tmp[r1][1] = "%s/%s" % (folder, f)
				else:
					tmp[r1] = ["", "%s/%s" % (folder, f), sample]
sys.stderr.write("Found %d files: %d forward reads, %d reverse reads\n" % (cnt[2], cnt[0], cnt[1]))
for k in tmp:
	files.append(tmp[k])
	
out = open("%s/demultiplexed.cleaned.%s" % (args.folder, "fasta" if args.fasta is True else "fastq"), "w+")
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

for col in files:
	sample = col[2]
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
		sys.stderr.write("Unknown nucleotide(s) used for sample '%s', sample entry ignored\n" % sample)
		continue		

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
			if args.ignore_bases:
				if len(primer) > args.ignore_bases:
					primer = primer[args.ignore_bases:]
			flookup_lens[len(primer)] = True
			flookup[primer] = True
		for primer in rprimers:
			if args.ignore_bases:
				if len(primer) > args.ignore_bases:
					primer = primer[args.ignore_bases:]
			rlookup_lens[len(primer)] = True
			rlookup[primer] = True
	else:
		for primer in fprimers:
			if args.ignore_bases:
				if len(primer) > args.ignore_bases:
					primer = primer[args.ignore_bases:]
			flookup_lens[len(primer)] = True
			for i in range(len(primer)):
				for nucleotide in nucleotides:
					flookup["%s%s%s" % (primer[0:i], nucleotide, primer[i + 1:])] = True
		for primer in rprimers:
			if args.ignore_bases:
				if len(primer) > args.ignore_bases:
					primer = primer[args.ignore_bases:]
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
	if len(col[0]) > 0 and len(col[1]) > 0 and args.force_forward is False and args.force_reverse is False:
		if args.folder in col[0]:
			name1 = col[0]
			name2 = col[1]
		else:
			name1 = args.folder + "/" + col[0]
			name2 = args.folder + "/" + col[1]
		i = 0
		files_f += 1
		files_r += 1
		if packedFile(name1):
			fh1 = gzip.open(name1, "r")
			file1_packed = True
		else:
			fh1 = open(name1, "r")
			file1_packed = False
		if packedFile(name2):
			fh2 = gzip.open(name2, "r")
			file2_packed = True
		else:
			fh2 = open(name2, "r")
			file2_packed = False
		with fh1 as file1, fh2 as file2:
			for f1, f2 in zip(file1, file2):
				# if it is byte type object in Python 3+ for packed files, decode to string
				if file1_packed:
					f1 = f1.decode()
				if file2_packed:
					f2 = f2.decode()
				f1 = f1.strip()
				f2 = f2.strip()
				if i % 4 == 0:
					seq_i += 1
					if seq_i % 100000 == 0:
						sys.stderr.write("%d/%d sequences cleaned/parsed from %d/%d forward/reverse files\n" % (seq_sel, seq_i, files_f, files_r))
					header = f1[1:].split(" ")[0]
				elif i % 4 == 1:
					if args.ignore_bases:
						f1 = f1[args.ignore_bases:]
						f2 = f2[args.ignore_bases:]
					sequence1 = f1
					sequence2 = f2
				elif i % 4 == 3:
					if args.ignore_bases:
						f1 = f1[args.ignore_bases:]
						f2 = f2[args.ignore_bases:]
					selected = True
					if args.trimq is not None and args.trimw is not None:
						trim_len = calculateSlidingWindow(f1, args.trimq, args.trimw, args.phred)
						if trim_len < len(f1):
							f1 = f1[:trim_len]
							sequence1 = sequence1[:trim_len]
						trim_len = calculateSlidingWindow(f2, args.trimq, args.trimw, args.phred)
						if trim_len < len(f2):
							f1 = f1[:trim_len]
							sequence2 = sequence2[:trim_len]
					fprimer_len = 0
					rprimer_len = 0
					selected = False
					if len(fprimer) == 0 and len(rprimer) == 0:
						selected = True
					if len(fprimer) > 0:
						for l in flookup_lens:
							if sequence1[:l] in flookup:
								selected = True
								fprimer_len = l
								break
							if args.primers_mixed:
								if sequence2[:l] in flookup:
									selected = True
									fprimer_len = l
									# swap sequences around to print them out correctly
									sequence1, sequence2 = sequence2, sequence1
									f1, f2 = f2, f1
									break
					if len(rprimer) > 0:
						for l in rlookup_lens:
							if sequence2[:l] in rlookup:
								selected = True
								rprimer_len = l
								break
					if args.remove_primer:
						sequence1 = sequence1[fprimer_len:]
						f1 = f1[fprimer_len:]
						sequence2 = sequence2[rprimer_len:]
						f2 = f2[rprimer_len:]
					if args.forward_trim > 0:
						sequence1 = sequence1[0:args.forward_trim]
						f1 = f1[0:args.forward_trim]
					if args.reverse_trim > 0:
						sequence2 = sequence2[0:args.reverse_trim]
						f2 = f2[0:args.reverse_trim]
					if adapters:
						found = False
						for adapter in adapters:
							pos = sequence1.find(args.adapter)
							if pos > 0:
								sequence1 = sequence1[:pos]
								f1 = f1[:pos]
								found = True
							pos = sequence2.find(args.adapter)
							if pos > 0:
								sequence2 = sequence2[:pos]
								f2 = f2[:pos]
								found = True
							if found:
								break
					# trim based on minimum base criteria
					if args.min_base_trimmed > 0:
						for _i in range(len(f1)):
							if ord(f1[_i]) - args.phred < args.min_base_trimmed:
								sequence1 = sequence1[:_i - 1]
								f1 = f1[:_i - 1] 
								break
						if f2:
							for _i in range(len(f2)):
								if ord(f2[_i]) - args.phred < args.min_base_trimmed:
									sequence2 = sequence2[:_i - 1]
									f2 = f2[:_i - 1]
									break
					# include based on minimum base criteria
					if args.min_allowed_base > 0:
						for _i in f1:
							if ord(_i) - args.phred < args.min_allowed_base:
								selected = False
								break
						if f2 and selected:
							for _i in f2:
								if ord(_i) - args.phred < args.min_allowed_base:
									selected = False
									break						
					if selected and (avgQuality(f1, args.phred) + avgQuality(f2, args.phred)) / 2 < args.q:
						selected = False
					if args.ml is not None and (len(sequence1) < args.ml or len(sequence2) < args.ml):
						selected = False
					if len(homopolymers) > 0:
						for homopolymer in homopolymers:
							if homopolymer in sequence1:
								# truncate homopolymer with quality
								_s = ""
								_q = ""
								_chr = ""
								_chr_i = 0
								_pos = -1
								for _ in sequence1:
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
										_q += f1[_pos]
								sequence1 = _s
								if not args.fasta:
									f1 = _q
								break
						for homopolymer in homopolymers:
							if homopolymer in sequence2:
								# truncate homopolymer with quality
								_s = ""
								_q = ""
								_chr = ""
								_chr_i = 0
								_pos = -1
								for _ in sequence2:
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
										_q += f2[_pos]
								sequence2 = _s
								if not args.fasta:
									f2 = _q
								break
					if selected:
						seq_sel += 1
						if args.fasta:
							out.write(">%s|%s\n%s\n" % (col[2], header, sequence1))
							out.write(">%s|%s\n%s\n" % (col[2], header, sequence2))
						else:
							out.write("@%s|%s\n%s\n+\n%s\n" % (col[2], header, sequence1, f1))
							out.write("@%s|%s\n%s\n+\n%s\n" % (col[2], header, sequence2, f2))
				i += 1
	else:
		if col[1] and args.force_reverse:
			name = col[1]
		else:
			name = col[0]
		i = 0
		files_f += 1
		if name is '':
			continue
		file1_packed = False
		if packedFile(name):
			fh = gzip.open(name, "r")
			file1_packed = True
		else:
			fh = open(name, "r")
		
		for line in fh:
			line = line.strip()
			if file1_packed:
				line = line.decode()
			if i % 4 == 0:
				seq_i += 1
				if seq_i % 100000 == 0:
					sys.stderr.write("%d/%d sequences cleaned/parsed from %d/%d forward/reverse files\n" % (seq_sel, seq_i, files_f, files_r))
				header = line[1:].split(" ")[0]
			elif i % 4 == 1:
				if args.ignore_bases:
					line = line[args.ignore_bases:]
				sequence = line
			elif i % 4 == 3:
				if args.ignore_bases:
					line = line[args.ignore_bases:]
				if args.trimq is not None and args.trimw is not None:
					trim_len = calculateSlidingWindow(line, args.trimq, args.trimw, args.phred)
					if trim_len < len(line):
						line = line[:trim_len]
						sequence = sequence[:trim_len]
				fprimer_len = 0
				selected = False
				if len(fprimer) > 0:
					for l in flookup_lens:
						if sequence[:l] in flookup:
							selected = True
							fprimer_len = l
							break
				else:
					selected = True
				if args.remove_primer:
					sequence = sequence[len(args.forward_primer):]
					line = line[len(args.forward_primer):]
				if args.forward_trim > 0:
					sequence = sequence[:args.forward_trim]
					line = line[:args.forward_trim]
				if adapters:
					for adapter in adapters:
						pos = sequence.find(adapter)
						if pos > 0:
							sequence = sequence[:pos]
							line = line[:pos]
							break
				if args.min_base_trimmed > 0:
					for _i in range(len(line)):
						if ord(line[_i]) - args.phred < args.min_base_trimmed:
							sequence = sequence[:_i - 1]
							line = line[:_i - 1] 
							break
				# include based on minimum base criteria
				if args.min_allowed_base > 0:
					for _i in line:
						if ord(_i) - args.phred < args.min_allowed_base:
							selected = False
							break
				if selected and avgQuality(line, args.phred) < args.q:
					selected = False
				if args.ml is not None and len(sequence) < args.ml:
					selected = False
				if len(homopolymers) > 0:
					for homopolymer in homopolymers:
						if homopolymer in sequence:
							# truncate homopolymer with quality
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
						out.write(">%s|%s\n%s\n" % (col[2], header, sequence))
					else:
						out.write("@%s|%s\n%s\n+\n%s\n" % (col[2], header, sequence, line))
			i += 1

sys.stderr.write("%d/%d sequences cleaned/parsed from %d/%d forward/reverse files\n" % (seq_sel, seq_i, files_f, files_r))

if out:
	out.close()
