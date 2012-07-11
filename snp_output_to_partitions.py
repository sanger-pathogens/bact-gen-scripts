#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=3) or '-h' in sys.argv[1:]:
	print "snp_output_to_partitions.py <snp output file> <output file name>"
	sys.exit()

lines=open(sys.argv[1], "rU").readlines()[1:]



syn=[]
non=[]
inter=[]

for x, line in enumerate(lines):
	y=x+1
	if line.split()[3].replace('/','').replace('1','') in ['N','2']:
		non.append(str(y))
	elif line.split()[3].replace('/','').replace('1','')=='S':
		syn.append(str(y))
	elif line.split()[3].replace('/','').replace('1','')=='-':
		inter.append(str(y))
	
output=open(sys.argv[2], "w")

print >> output, "DNA, synonymous = "+', '.join(syn)
print >> output, "DNA, nonsynonymous = "+', '.join(non)
print >> output, "DNA, intergenic = "+', '.join(inter)

output.close()