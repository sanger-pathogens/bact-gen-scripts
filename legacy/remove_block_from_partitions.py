#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=5) or '-h' in sys.argv[1:]:
	print "remove_block_from_partitions.py <partitions file> <first base in region to remove> <last base in region to remove> <output file name>"
	sys.exit()


lines=open(sys.argv[1], "rU").readlines()

start=int(sys.argv[2])

end=int(sys.argv[3])

output=open(sys.argv[4], 'w')

for line in lines:
	newline= ' '.join(line.strip().split()[:3])
	words=line.strip().split()[3:]
	toadd=[]
	for word in words:
		if int(word.replace(',',''))<start:
			toadd.append(word.replace(',',''))
		elif int(word.replace(',',''))>end:
			toadd.append(str(int(word.replace(',',''))-((end+1)-start)))
			
	newline=newline+', '.join(toadd)
	
	print >> output, newline