#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math

if len(sys.argv)!=5 or '-h' in sys.argv[1:]:
	print "plot2tab.py <plot file> <cutoff value> <make tab from values above (a) or below (b) the cutoff> <output file>"
	sys.exit()

try: cutoff=float(sys.argv[2])
except ValueError:
	print 'Non-numeric value %s found in argument 2' % sys.argv[2]

abovebelow=sys.argv[3]
if abovebelow not in ['a','b']:
	print 'argument 3 must be a or b'
	sys.exit()

plotfile=sys.argv[1]
if not os.path.isfile(plotfile):
	print "Cannot find file:", plotfile

if plotfile.split('.')[-1]=='gz':
	gfeil=gzip.open(plotfile, 'r')
	lines=gfeil.readlines()
else:
	lines=open(plotfile, 'rU').readlines()

outfile=sys.argv[4]

overwrite='y'

while overwrite!='y' and os.path.isfile(outfile):
	overwrite=raw_input(outfile+' exists! Overwrite? (y/n): ')
	overwrite.lower()
	if overwrite=='n':
		outfile=raw_input('Enter a new output file name: ')


output=open(outfile,"w")

print >> output, "ID   misc"
count=1
inblock=['n', 'n', 'n', 'n']
featureline=['', '', '', '']

for line in lines:
	for x, lane in enumerate(line.strip().split()):
		if x>3:
			continue
		value=float(lane)
		
		if abovebelow=='a' and value>cutoff and inblock[x]=='n':
			featureline[x]="FT   misc_feature    "+str(count)+".."
			inblock[x]='y'
		if abovebelow=='b' and value<cutoff and inblock[x]=='n':
			featureline[x]="FT   misc_feature    "+str(count)+".."
			inblock[x]='y'
		elif abovebelow=='a' and value<=cutoff and inblock[x]=='y':
			featureline[x]=featureline[x]+str(count-1)
			print >> output, featureline[x]
			print >> output, 'FT                   /colour=2'
			#print >> output, 'FT                   /colour='+str(4-x)
			#print >> output, 'FT                   /note="above '+str(cutoff)+'"'
			inblock[x]='n'
		elif abovebelow=='b' and value>=cutoff and inblock[x]=='y':
			featureline[x]=featureline[x]+str(count-1)
			print >> output, featureline[x]
			print >> output, 'FT                   /colour=2'
			#print >> output, 'FT                   /colour='+str(4-x)
			#print >> output, 'FT                   /note="below '+str(cutoff)+'"'
			inblock[x]='n'
	
	count=count+1
for x, block in enumerate(inblock):
	if block=='y':
		featureline[x]=featureline[x]+str(count-1)
		print >> output, 'FT                   /colour='+str(4-x)
		print >> output, featureline[x]

output.close()
	
