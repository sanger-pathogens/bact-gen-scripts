#!/usr/bin/env python
import string, re, sets
import os, sys, getopt, random, math


if (len (sys.argv)!=7) or '-h' in sys.argv[1:]:
	print "location_to_mapping.py <text files with lines to extract> <length of reference (TW20)> <core tab file> <ref name(TW20)> <alignment> <output file name>"
	sys.exit()

locations=open(sys.argv[1], "rU").readlines()

lines=open(sys.argv[3], "rU").read().split('\n')[1:]

noncore=sets.Set()

start=1
end=1
for line in lines:
	if len(line.strip().split())<3:
		continue
	start=int(line.strip().split()[2].split("..")[0])
	if start>end:
		for x in range(end,start):
			noncore.add(x)
	print start, end
	end=int(line.strip().split()[2].split("..")[1])

if start>end:
	for x in range(end,int(sys.argv[2])+1):
		noncore.add(x)

print len(noncore)


lines=open(sys.argv[5], "rU").read().split('>')[1:]

seqs={}
nonref=''

for line in lines:
	seqs[line.split('\n')[0].split()[0]]=''.join(line.split('\n')[1:])
	if line.split('\n')[0].split()[0]!=sys.argv[4]:
		nonref=line.split('\n')[0].split()[0]


sitetranslation=[]
posinref=0


for x, site in enumerate(seqs[nonref]):
	if seqs[sys.argv[4]][x]!='-':
		posinref=posinref+1
		if posinref-1 in noncore:
			if site!='-':
				sitetranslation.append("n")
		else:
			if site!='-':
				sitetranslation.append("c")
	elif site!='-':
		sitetranslation.append('-')



output=open(sys.argv[6], "w")

print len(sitetranslation), len(seqs[nonref].replace('-',''))


for location in locations:
	location=location.strip()
	if sitetranslation[int(location)-1]=='n':
		print >>output, location, "noncore"
	elif sitetranslation[int(location)-1]=='c':
		print >>output, location, "core"
	elif sitetranslation[int(location)-1]=='-':
		print >>output, location, "unmapped"
		
	

output.close()