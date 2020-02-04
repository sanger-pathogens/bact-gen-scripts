#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=4) or '-h' in sys.argv[1:]:
	print "location_to_mapping.py <text files iwth lines to extract> <alignment file> <output file name>"
	sys.exit()

locations=open(sys.argv[1], "rU").readlines()

lines=open(sys.argv[2], "rU").read().split('>')[1:]

seqs={}

for line in lines:
	seqs[line.split('\n')[0].split()[0]]=''.join(line.split('\n')[1:])


output=open(sys.argv[3], "w")

names=seqs.keys()

names.sort()

print >> output, "position\t"+'\t'.join(names)

for location in locations:
	outstring=location.strip()
	for name in names:
		outstring=outstring+'\t'+seqs[name][int(location)-1]
	print >> output, outstring
	

output.close()