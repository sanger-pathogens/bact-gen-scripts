#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


if (len (sys.argv)!=4) or '-h' in sys.argv[1:]:
	print "parse_paup_site_reconstructions.py <snp file> <paup log file> <outfile prefix>"
	sys.exit()


snplines=open(sys.argv[1],'rU').readlines()

posnlist=[]

for x, line in enumerate(snplines):
	posnlist.append(line.split()[0])



pauplines=open(sys.argv[2],'rU').readlines()

ambiguous={}
unambiguous={}
prevCI=''

for x, line in enumerate(pauplines):
	words=line.split()
	if (words[5]=='==>' and words[1]=='1.000'):
		if not unambiguous.has_key(words[3]+words[5]+words[7]):
			unambiguous[words[3]+words[5]+words[7]]=[]
		unambiguous[words[3]+words[5]+words[7]].append(posnlist[int(words[0])])
		prevCI=words[1]
		prevchar=words[0]
	elif (words[3]=='==>' and prevCI=='1.000'):
		if not unambiguous.has_key(words[1]+words[3]+words[5]):
			unambiguous[words[1]+words[3]+words[5]]=[]
		unambiguous[words[1]+words[3]+words[5]].append(posnlist[int(prevchar)])


keys=unambiguous.keys()
keys.sort()

output =open(sys.argv[3],'w')

for key in keys:
	print >> output, key+'\t'+','.join(unambiguous[key])
	
output.close()