#!/usr/bin/env python
from __future__ import division

import os
import argparse
import console
import fileinput
import sys
import glob
import subprocess

"""
	Runs multiple BLAST+ programs sequentially to overcome problem with long command parameter error under Windows and track the progress when partitioned databases are used.
"""

parser = argparse.ArgumentParser(description = """
	Runs multiple BLAST+ programs sequentially to overcome problem with long command parameter error under Windows and track the progress when partitioned databases are used.
	""")
parser.add_argument(
	'-f', metavar = 'FASTA_FILE', required = True, type = str, help = """
	specify an input FASTA for the BLAST+ program
	""")
parser.add_argument(
	'-db', metavar = 'DATABASE_FILE', required = True, type = str, help = """
	specifya  BLAST+ database with partition numbers, partitions are run sequentially. Use a higher partition than 00 to continue a terminated BLAST+ 
	""")
parser.add_argument(
	'-hits', metavar = 'HITS', required = False, type = int, default = 10, help = """
	specify a number of best hits for each query for each database partition using BLAST+ 
	""")
parser.add_argument(
	'-evalue', metavar = 'EVALUE', required = False, type = str, default = "1e-50", help = """
	specify a BLAST+ e-value threshold, smaller values allow longer alignments
	""")
parser.add_argument(
	'-t', metavar = 'THREADS', required = False, type = int, help = """
	specify the number of threads to be used by the BLAST+ program
	""")
parser.add_argument(
	'-program', metavar = 'PROGRAM', required = False, type = str, default = "blastn", help = """
	specify BLAST+ program binary file location. If installed as system wide variable or in same folder, use "blastn"
	""")
	
args = parser.parse_args()

# run BLASTs in partitions
file_name = os.path.basename(args.db)
directory_name = os.path.dirname(args.db)
tmp = file_name.replace(".nhr", "").split(".")
tmp.pop()   # remove index
partitions = glob.glob("%s/%s.*.nhr" % (directory_name, ".".join(tmp)))
ids = {}
ids_sorted = []
for partition in partitions:
	partition_i = int(partition.split(".")[-2])
	ids[partition_i] = partition.replace(".nhr", "")
	ids_sorted.append(partition_i)
ids_sorted.sort()
# check at what partition user wants to continue BLAST+
predefined_id = int(args.db.replace(".nhr", "").split(".")[-1])
console.log("Found %s partitions\n" % len(ids))
commands = []
for i in ids:
	if predefined_id <= i:
		console.log("Executing BLAST+ against partition %s\n" % i)
		# construct BLAST command
		command = '%s -query "%s" -dust no -evalue %s -max_target_seqs %s -num_threads %s -db "%s" -outfmt "6 qseqid sseqid stitle evalue pident nident length frames qstart qend sstart send qlen slen score" > "%s.%d.blast"' % (args.program,
			args.f,
			args.evalue,
			args.hits,
			args.t,
			ids[i],
			args.f,
			i
		)
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
			console.log("Process finished successfully\n")
		else:
			console.log("Process finished with error code %s\n" % process.returncode)