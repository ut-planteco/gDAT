from __future__ import division

import os
import argparse
import re
import sys
import itertools
import gzip
import array
import console
import subprocess

"""
	Combine reads using FLASh or vsearch inside a folder containing Illumina R1 and R2 reads.
"""
parser = argparse.ArgumentParser(description = """
	Combine reads using FLASh or vsearch inside a folder containing Illumina R1 and R2 reads.
	""")
parser.add_argument(
	'-o', metavar = 'ALLOW_OUTIES', default = 0, type = int, help = """
	allow overlap section to be on another end
	""")
parser.add_argument(
	'-m', metavar = 'MINIMUM_LENGTH', default = 10, type = int, help = """
	specify minimum overlap length
	""")
parser.add_argument(
	'-M', metavar = 'MAXIMUM_LENGTH', default = 300, type = int, help = """
	specify maximum overlap length
	""")
parser.add_argument(
	'-x', metavar = 'MISMATCHES', default = 0.25, type = float, help = """
	specify maximum allowed ratio between the number of 
	mismatched base pairs and the overlap length
	""")
parser.add_argument(
	'-t', metavar = 'THREADS', default = 1, type = int, help = """
	specify number of threads
	""")
parser.add_argument(
	'-fasta', metavar = 'FASTA', default = 0, type = int, help = """
	specify if generate output is FASTA or FASTQ
	""")
parser.add_argument(
	'-p', metavar = "PROGRAM", default = 1, type = int, help = """
	specify program to combine reads together:
		1 - FLASh
		2 - vsearch
	""")
parser.add_argument(
	'-folder', metavar = 'FOLDER', type = str, required = True, help = """
	specify folder containing FASTQ files, supports gz packed files
	""")
	
args = parser.parse_args()


folders = []

folders = [x[0] for x in os.walk(args.folder)]
cnt = [0, 0, 0]
tmp = {}
files = []

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
	
out = open("%s/combined.%s" % (args.folder, "fasta" if args.fasta is 1 else "fastq"), "w+")

for col in files:
	if len(col[0]) == 0:
		console.log("Found %s sample but forward read is missing\n" % col[2])
	elif len(col[1]) == 0:
		console.log("Found %s sample but reverse read is missing\n" % col[2])
	else:
		console.log("Found %s sample, combining %s and %s reads\n" % (col[2], col[0], col[1]))
		# check if packed files, unpack them
		if ".gz" in col[0]:
			fh = open("%s/tmp1" % args.folder, "w")
			first = True
			with gzip.open(col[0], "r") as f:
				for r in f:
					if sys.version_info[0] >= 3:
						r = r.decode("utf-8")
					if first:
						first = False
						# check if file name is first row
						if ".fastq" in r:
							r = "@%s" % r.split("@")[-1]
					if ord(r[0]) > 0:
						fh.write("%s" % r)
			fh.close()
			col[0] = "%s/tmp1" % args.folder

		if ".gz" in col[1]:
			fh = open("%s/tmp2" % args.folder, "w")
			first = True
			with gzip.open(col[1], "r") as f:
				for r in f:
					if sys.version_info[0] >= 3:
						r = r.decode("utf-8")
					if first:
						first = False
						# check if file name is first row
						if ".fastq" in r:
							r = "@%s" % r.split("@")[-1]
					if ord(r[0]) > 0:
						fh.write("%s" % r)
			fh.close()
			col[1] = "%s/tmp2" % args.folder

		# construct command
		if args.p == 1:
			command = 'flash %s -c -m %d -M %d -x %f -t %s "%s" "%s" > "%s/tmp"' % (
						"-O" if args.o else "",
						args.m,
						args.M,
						args.x,
						args.t,
						col[0],
						col[1],
						args.folder
						)
		else:
			command = 'vsearch --fastq_mergepairs "%s" --reverse "%s" --threads %s --fastq_minovlen %s %s --fastqout "%s/tmp"' % (
						col[0],
						col[1],
						args.t,
						args.m,
						"--fastq_allowmergestagger" if args.o else "",
						args.folder
						)
		if os.name != "nt":
			command = "./%s" % command
		console.log("%s\n" % command)
		if os.name == "nt" or sys.version_info[0] >= 3:
			process = subprocess.Popen(command, 
				stdout=subprocess.PIPE, 
				stderr=subprocess.PIPE,
				#encoding="utf8",
				shell=True
				#start_new_session=True
			)
		else:
			process = subprocess.Popen(command, 
				stdout=subprocess.PIPE, 
				stderr=subprocess.PIPE,
				shell=True,
				preexec_fn=os.setsid
			)
		process.wait()
		if process.returncode == 0:
			# combine results into file
			l = 0
			for line in open("%s/tmp" % args.folder):
				if args.fasta:
					if l % 4 == 0:
						out.write(">%s|%s\n" % (col[2], line.strip()[1:]))
					elif l % 4 == 1:
						out.write("%s\n" % line.strip())
				else:
					if l % 4 == 0:
						out.write("@%s|%s\n" % (col[2], line.strip()[1:]))
					else:
						out.write("%s\n" % line.strip())
				l += 1
			console.log("Process finished successfully\n")
		else:
			console.log("Process finished with error code %s\n" % process.returncode)
out.close()
