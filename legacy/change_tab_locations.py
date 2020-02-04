#!/usr/bin/env python

import os, sys

locs=sys.argv[1]
tabfile=sys.argv[2]

loc_dict={}
for line in open(locs, "rU"):
	words=line.strip().split()
	if len(words)!=2:
		continue
	loc_dict[words[0]]=words[1]
#
#print loc_dict
#sys.exit()

for line in open(tabfile, "rU"):
	line=line.strip()
	words=line.split()
	if len(words)==3 and words[0]=="FT" and words[1]=="SNP":
		if not words[2] in loc_dict:
			print "Could not find", words[2]
			sys.exit()
		print "FT   SNP             "+loc_dict[words[2]]
	else:
		print line
