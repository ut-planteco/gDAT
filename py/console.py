#!/usr/bin/env python
from __future__ import division

import sys
import datetime
import time
import tarfile

start = time.time()

def log(output):
	global start
	now = datetime.datetime.utcfromtimestamp(time.time() - start)
	sys.stderr.write("[%s] %s" % (now.strftime("%H:%M:%S"), output))
	sys.stderr.flush()

def unpack(name, folder):
	tar = tarfile.open(name)
	tar.extractall(folder)
	tar.close()

