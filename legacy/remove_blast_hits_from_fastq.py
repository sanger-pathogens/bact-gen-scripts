#!/usr/bin/env python

import os, sys
from optparse import OptionParser



fastq=open(sys.argv[1])

if len(sys.argv)>2:
	evalue=sys.argv[2]

os.system("blastall -p blastn -v 1 -b 1 -m 8 -o "+fastq+".blast -e "+evalue


blast=open(sys.argv[2])

blastlist=set()
for line in blast:
	if len(blastlist)==0 or line.split()[0][:-2] not in blastlist:
		blastlist.add(line.split()[0][:-2])


output=open(sys.argv[1].replace(".fastq", "_filtered.fastq"),"w")

count=0
for line in fastq:
	if line.strip().replace("@","")[:-2] not in blastlist:
		print >> output, line.strip()
		for x in range(0,3):
			print >> output, fastq.next().strip()
		count+=1
	else:
		for x in range(0,3):
			output, fastq.next().strip()
			
print count, "reads remaining"