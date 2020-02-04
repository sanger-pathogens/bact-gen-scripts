#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time

if len(sys.argv)!=3:
	print 'dindel_2_tab.py <dindel vcf file> <output tab file name>'

try:
	infile=open(sys.argv[1], 'rU')
except StandardError:
	print 'Cannot open file', sys.argv[1]
	sys.exit()

output=open(sys.argv[2], 'w')

print >> output, 'ID   SV'

for line in infile:
	if len(line)==0 or line[0]=='#':
		continue
	
	
	words=line.strip().split()
	
	if len(words[3])==1:
		indeltype='I'
	elif len(words[4])==1:
		indeltype='D'
	else:
		indeltype="?"
	
	if indeltype=='I':
		location=str(int(words[1]))+'..'+str(int(words[1])+1)
	else:
		location=str(int(words[1])+1)+'..'+str((int(words[1]))+(len(words[3])-len(words[4])))
	
	
	print >> output, 'FT   '+indeltype+' '*(16-len(indeltype))+location
	
	print >> output, 'FT                   /ref="'+words[3]+'"'
	print >> output, 'FT                   /alt="'+words[4]+'"'
	print >> output, 'FT                   /quality="'+words[5]+'"'
	print >> output, 'FT                   /filter="'+words[6]+'"'
	print >> output, 'FT                   /info="'+words[7]+'"'
	print >> output, 'FT                   /format="'+words[8]+'"'
	print >> output, 'FT                   /sample="'+words[9]+'"'
	
	if indeltype=='D':
		print >> output, 'FT                   /colour=6'
	elif indeltype=='I':
		print >> output, 'FT                   /colour=4'
	
	


output.close()
	

	