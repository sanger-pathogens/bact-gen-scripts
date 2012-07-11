#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math


if (len (sys.argv)!=3 and len (sys.argv)!=5) or '-h' in sys.argv[1:]:
	print "Find_conserved_block_from_covcount.py <covcount plot file> <number of taxa> > <output file>"
	sys.exit()


covcountfile=sys.argv[1]
taxacutoff=int(sys.argv[2])

if covcountfile.split('.')[-1]=='gz':
	lines=gzip.open(covcountfile, "r").readlines()
else:
	lines=open(covcountfile, "rU").readlines()

incore='a'
mintochangetocore=100
mintochangetoexcluded=1000

excludecount=0
corecount=0
curcorestart=0
curcoreend=0

output=open('core_to_whole.csv','w')

print "ID   CORE"

coreposn=0
for x, line in enumerate(lines):
	line=int(line.strip())
	#print line, incore
	if line>=taxacutoff:
		corecount=corecount+1
		if corecount==1 and incore!='y':
			curcorestart=x+1
		excludecount=0
		if incore!='y' and corecount>=mintochangetocore:

			incore='y'
	else:
		excludecount=excludecount+1
		if excludecount==1 and incore!='n':
			curcoreend=x+1
		corecount=0
		if incore!='n' and excludecount>=mintochangetoexcluded:
			incore='n'
			print "FT   CORE            "+str(curcorestart)+'..'+str(curcoreend)
			
			for y in range(curcorestart, curcoreend+1):
				coreposn=coreposn+1
				print >> output, str(coreposn)+","+str(y)
			
			#sys.exit()

if incore=='y':
	print "FT   CORE            "+str(curcorestart)+'..'+str(x+1)
	for y in range(curcorestart, curcoreend+1):
		coreposn=coreposn+1
		print >> output, str(coreposn)+","+str(y)


output.close()






