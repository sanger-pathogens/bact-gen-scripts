#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math

if len (sys.argv)!=2 or '-h' in sys.argv[1:]:
	print "MUMmer_snps_to_tab.py <tiling file from Abacus> > <outfile name>"

snps=sys.argv[1]	

lines=open(snps, "rU").readlines()

print "ID   SNPs"

for line in lines:
	words=line.strip().split()
	
	if words[1]=='.':
		print 'FT   insertion       '+words[0]+'..'+words[0]
		print 'FT                   /contig="'+words[14]+'"'
		print 'FT                   /colour="4"'
	elif words[2]=='.':
		print 'FT   deletion        '+words[0]+'..'+words[0]
		print 'FT                   /contig="'+words[14]+'"'
		print 'FT                   /colour="3"'
	else:
		print 'FT   SNP             '+words[0]+'..'+words[0]
		print 'FT                   /contig="'+words[14]+'"'
		print 'FT                   /colour="2"'
