#!/usr/bin/env python

import os, sys
from optparse import OptionParser, OptionGroup
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_SeqIO import *
import subprocess
import shlex
import gzip
import shutil
from random import *
import math

##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


dirname=sys.argv[1]
tmpname=sys.argv[2]
output_file_name=sys.argv[3]
number_of_blast_files=int(sys.argv[4])
remove_self_matches=sys.argv[5]

output=gzip.open(output_file_name+".crunch.gz", "w")

subject_db={}
for line in open(dirname+"/"+tmpname+".subjectinfo", "rU"):
	words=line.strip().split("\t")
	subject_db[words[0]]=[words[1], int(words[2])]

query_db={}
for line in open(dirname+"/"+tmpname+".queryinfo", "rU"):
	words=line.strip().split("\t")
	query_db[words[0]]=[words[1], int(words[2])]


for blast_number in xrange(1,number_of_blast_files+1):
	blast_filename=dirname+"/"+tmpname+".blast"+str(blast_number)
	for line in open(blast_filename, "rU"):
		words=line.strip().split()
		
		if remove_self_matches!="True" or query_db[words[0]][0].split()[0]!=subject_db[words[1]][0].split()[0]:
			print >> output, words[11], words[2], query_db[words[0]][1]+int(words[6]), query_db[words[0]][1]+int(words[7]), query_db[words[0]][0].split()[0], subject_db[words[1]][1]+int(words[8]), subject_db[words[1]][1]+int(words[9]), subject_db[words[1]][0].split()[0]


#if os.path.isdir(options.tmpdir):
#	shutil.rmtree(options.tmpdir)
#elif os.path.exists(options.tmpdir):
#	os.remove(options.tmpdir)

output.close()
