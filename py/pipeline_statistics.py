#!/usr/bin/env python
from __future__ import division

import os
import argparse
import gzip
import console
import fileinput
import sys
import operator
import array

"""
	Summarize FASTQ files that are located in a folder.
"""
parser = argparse.ArgumentParser(description = """
	Summarize FASTQ files that are located in a folder.
	""")
parser.add_argument(
	'-f', metavar = 'BASE_FOLDER', required = False, type = str, help = """
	specify a folder location for FASTA/FASTQ files, selects all files with fastq, fq, fasta or their packed gz version
	""")
parser.add_argument(
	'-i', metavar = 'BASE_FILES', required = False, type = str, help = """
	specify multiple FASTA/FASTQ files separating by commas
	""")
parser.add_argument(
	'-kmerlen', metavar = 'KMER_LENGTH', required = False, type = int, default = 0, help = """
	specify the kmer length to calculate barcode and/or primer distribution
	""")
parser.add_argument(
	'-kmershow', metavar = 'KMER_COUNT', required = False, type = int, default = 50, help = """
	specify a number of kmers to show in the output, use 0 to show all of them
	""")
parser.add_argument(
	'-phred', metavar = 'PHRED_SCORE', required = False, type = int, default = 33, help = """
	specify the PHRED score used for FASTQ files
	""")
parser.add_argument(
	'-separatestats', metavar = '[0|1]', required = False, type = int, default = 0, help = """
	specify as 1 to write separate statistics for each file
	""")
parser.add_argument(
	'-skip', metavar = '[0|1]', required = False, type = int, default = 0, help = """
	specify as 1 to count only every 10th sequence to speed up the process
	""")
parser.add_argument(
	'-quality', metavar = '[0|1]', required = False, type = int, default = 0, help = """
	specify as 1 to calculate the quality distribution
	""")
parser.add_argument(
	'-kmerfront', metavar = '[0|1]', required = False, type = int, default = 0, help = """
	specify as 1 to calculate kmers only at the front of the sequence
	""")
parser.add_argument(
	'-kmerpos', metavar = 'KMER_POS', required = False, type = int, default = 1, help = """
	specify at which position kmers should be counted. This is useful when having barcode at the front of the sequence and want to ignore barcode information and check the primer information
""")

args = parser.parse_args()


def packedFile(filename):
	ext = filename.split(".")[-1]
	if ext.lower() == "gz":
		return True
	else:
		return False

def sumQual(qual):
	if sys.version_info >= (3, 0):
		arr = array.array('b')
		try:
			arr.frombytes(qual.encode())
		except:
			arr.frombytes(qual)
	else:
		arr = array.array('b', qual)
	return sum(arr)

def calcQual(qual):
	if sys.version_info >= (3, 0):
		arr = array.array('b')
		try:
			arr.frombytes(qual.encode())
		except:
			arr.frombytes(qual)
	else:
		arr = array.array('b', qual)
	return (sum(arr) / float(len(qual))) - args.phred

KMER = args.kmerlen
QSTEP = 10
kmers = {}
skip = False
if args.skip:
	skip = True
skipsize = 40
avgqsum = 0
avgqi = 0
qptotal = []
qpcount = []
maxlen = 0
wformat = False

if args.f:
	files = []
	for filename in os.listdir(args.f):
		if ".txt" in filename or ".zip" in filename:
			continue
		if ".fq" in filename or ".fastq" in filename:
			files.append("%s/%s" % (args.f, filename))
	args.i = ",".join(files)
if args.i:
	files = args.i.split(",")
	for f in files:
		if packedFile(f):
			fh = gzip.open(f, "r")
		else:
			fh = open(f, "r")
		l = 0
		avgqsum_single = 0
		avgqi_single = 0
		qptotal_single = []
		kmers_perfile = {}
		qptotal_perfile = []
		qpcount_perfile = []
		maxlen_perfile = 0
		for line in fh:
			if wformat or l % 4 == 0:
				# check if proper FASTQ file is used
				if len(line) > 0 and line[0] == "@":
					wformat = False
				else:
					wformat = True
			if l % 4 == 1:
				# build kmer information
				if skip and l % skipsize <= 4 or not skip:
					if args.kmerfront:
						index = line[args.kmerpos - 1:KMER + args.kmerpos - 1]
						if index not in kmers_perfile:
							kmers_perfile[index] = 0
						if index not in kmers:
							kmers[index] = 0
						kmers[index] += 1
						kmers_perfile[index] += 1
					else:
						for i in range(args.kmerpos - 1, len(line) - KMER, KMER):
							index = line[i:i + KMER]
							if index not in kmers_perfile:
								kmers_perfile[index] = 0
							if index not in kmers:
								kmers[index] = 0
							kmers[index] += 1
							kmers_perfile[index] += 1
			elif l % 4 == 3:
				if skip and l % skipsize <= 4 or not skip:
					# build quality information
					if args.quality:
						t = sumQual(line)
						avgqsum += t
						avgqi += len(line)
						avgqsum_single += t
						avgqi_single += len(line)
						if len(line) > maxlen:
							maxlen = len(line)
							for i in range(0, maxlen, QSTEP):
								if i // QSTEP not in qptotal:
									qptotal.append(0)
									qpcount.append(0)
						if len(line) > maxlen_perfile:
							maxlen_perfile = len(line)
							for i in range(0, maxlen_perfile, QSTEP):
								if i // QSTEP not in qptotal_perfile:
									qptotal_perfile.append(0)
									qpcount_perfile.append(0)
						for i in range(0, len(line), QSTEP):
								q = calcQual(line[i:i + QSTEP])
								qptotal[i // QSTEP] += q
								qpcount[i // QSTEP] += 1
								qptotal_perfile[i // QSTEP] += q
								qpcount_perfile[i // QSTEP] += 1
			#if l % 10000 == 0:
			#	sys.stderr.write("%s: %s\r" % (f, l))
			l += 1
		q = 0
		tmp_kmer = ""
		tmp_kmer_count = 0
		if fh:
			fh.close()
			if args.separatestats:
				fh = open("%s.kmer%d.stats.txt" % (f, KMER), "wb+")
				fh.write("Kmers distribution using kmer = %d\n" % KMER)
				keys = sorted(kmers_perfile.items(), key = lambda x : x[1], reverse = True)
				i = 0
				for k in keys:
					if i == 0:
						tmp_kmer = k[0]
						tmp_kmer_count = k[1] * (10 if args.skip else 1)
					i += 1
					if args.kmershow > 0:
						if i > args.kmershow:
							break
					fh.write("%s\t%s\n" % (k[0], k[1] * (10 if args.skip else 1)))
				if args.quality:
					q = (avgqsum_single / avgqi_single - args.phred)
					fh.write("Average quality: %.3f\n" % q)
					fh.write("Pos\tAvg. q\tTotal\tCount\n")
					for i in range(0, maxlen_perfile, QSTEP):
						qa = qptotal_perfile[i // QSTEP] / qpcount_perfile[i // QSTEP]
						if qptotal_perfile[i // QSTEP] > 0:
							fh.write(("%d-%d\t%.3f\t%d\t%d\n" % (i + 1, i + QSTEP, qa, qptotal_perfile[i // QSTEP], qpcount_perfile[i // QSTEP] * (10 if args.skip else 1))))
				fh.close()
		if q > 0:
			sys.stderr.write("%s: %s sequences, %.3f avg, %s %d \n" % (f, l, q, tmp_kmer, tmp_kmer_count))
		else:
			sys.stderr.write("%s: %s sequences\n" % (f, l))

if len(kmers) > 0:
	print("Kmers distribution using kmer = %d" % args.kmerlen)
	keys = sorted(kmers.items(), key = lambda x : x[1], reverse = True)
	i = 0
	for k in keys:
		i += 1
		if args.kmershow > 0:
			if i > args.kmershow:
				break
		print("%s\t%s" % (k[0], k[1] * (10 if args.skip else 1)))
	if args.quality:
		print("Average quality: %.3f" % (avgqsum / avgqi - args.phred))
		print("Pos\tAvg. q\tTotal\tCount")
		for i in range(0, maxlen, QSTEP):
			qa = qptotal[i // QSTEP] / qpcount[i // QSTEP]
			if qptotal[i // QSTEP] > 0:
				print("%d-%d\t%.3f\t%d\t%d" % (i + 1, i + QSTEP, qa, qptotal[i // QSTEP], qpcount[i // QSTEP] * (10 if args.skip else 1)))
else:
	sys.stderr.write("Could not find FASTQ files or files are empty\n")

