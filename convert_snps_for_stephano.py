#!/usr/bin/env python


##################
# Import modules #
##################

import string, re, os, sys

lines=open(sys.argv[1],"rU").readlines()

headings=lines[0].strip().split("\t")
taxanames=headings[10:]
print '\t'.join(map(str,["position", "isolate", "ref", "SNP"]))
for line in lines[1:]:
	words=line.strip().split("\t")
	ref=words[7]
	pos=words[1]
	taxa=words[10:]
	for x, taxon in enumerate(taxa):
		if taxon!=ref and taxon not in ["?", "-", "N", "."]:
			print '\t'.join(map(str,[pos, taxanames[x], ref, taxon]))
	