#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math

if len (sys.argv)!=2 or '-h' in sys.argv[1:]:
	print "MUMmer_tiling_to_tab.py <tiling file from Abacus> > <outfile name>"
	sys.exit()

crunch=sys.argv[1]	

lines=open(crunch, "rU").readlines()

print "ID   Contigs"

reflen=int(lines[0].strip().split()[1])

for line in lines[1:]:
	words=line.strip().split()
	
	if int(words[0])<0:
		start=str(reflen-int(words[0][1:]))
		words[0]='1'
		if words[6]=='-':
			print "FT   contig          complement("+start+'..'+str(reflen)+')'
		else:
			print "FT   contig          "+start+'..'+end
		print 'FT                   /systematic_id="'+words[7]+'b"'
		print 'FT                   /colour="2"'
		print 'FT                   /method="mummer"'
	elif int(words[1])>reflen:
		end=str(int(words[1])-reflen)
		words[1]=str(reflen)
		if words[6]=='-':
			print 'FT   contig          complement(1..'+end+')'
		else:
			print 'FT   contig          1..'+end
		print 'FT                   /systematic_id="'+words[7]+'a"'
		print 'FT                   /colour="2"'
		print 'FT                   /method="mummer"'

	if words[6]=='-':
		print 'FT   contig          complement('+words[0]+'..'+words[1]+')'
	else:
		print 'FT   contig          '+words[0]+'..'+words[1]
	print 'FT                   /systematic_id="'+words[7]+'"'
	print 'FT                   /colour="2"'
	print 'FT                   /method="mummer"'
