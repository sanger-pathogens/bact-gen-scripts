#!/usr/bin/env python
import string, re, gzip
import os, sys, math


if len(sys.argv)!=3 or '-h' in sys.argv[1:]:
	print "cigar2_odd_plot.py <cigar file (can be zipped)> <output file name> <max allowed insert size>"
	sys.exit()


print "Reading cigar file...",
sys.stdout.flush()
if sys.argv[1].split('.')[-1]=='gz':
	lines=gzip.open(sys.argv[1], 'r').readlines()
else:
	lines=open(sys.argv[1], 'r').readlines()
print "Done."
print "Calculating insert sizes..."
sys.stdout.flush()

numlines=len(lines)
count=0
total=0

inserts=[]

for loc, line in enumerate(lines):
	count=count+1
	if count>1000:
		total=total+count
		count=0
		print str(int((float(total)/numlines)*100))+"% complete\r",
	
	words=line.split()
	if len(words)<5:
		continue
	#Identify pairs that map
	if loc>0 and words[1].replace(".R","").replace("/2","")==lines[loc-1].split()[1].replace(".F","").replace("/1",""):
		#if insert is too big or too small add positions to the appropriate list
		if words[4]=="-" and lines[loc-1].split()[4]=="+":
#			if (int(words[6])-int(lines[loc-1].split()[7]))<100000 and int(words[6])-int(lines[loc-1].split()[7])>0:
#				inserts.append(int(words[6])-int(lines[loc-1].split()[7]))
			inserts.append(int(words[6])-int(lines[loc-1].split()[7]))

		elif words[4]=="+" and lines[loc-1].split()[4]=="-":
#			if (int(lines[loc-1].split()[6])-int(words[7]))<100000 and int(lines[loc-1].split()[6])-int(words[7])>0:
#				inserts.append(int(lines[loc-1].split()[6])-int(words[7]))
			inserts.append(int(lines[loc-1].split()[6])-int(words[7]))
				
print "100% complete"	
print "Printing output file...",
sys.stdout.flush()	


output=open(sys.argv[2],'w')

for x in inserts:
	print >> output, x

output.close()
print "Done."


