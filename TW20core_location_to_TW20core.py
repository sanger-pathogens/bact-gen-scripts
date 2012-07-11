#!/usr/bin/env python
import string, re, sets
import os, sys, getopt, random, math


if (len (sys.argv)!=5) or '-h' in sys.argv[1:]:
	print "location_to_mapping.py <text files with lines to extract> <length of reference (TW20)> <core tab file> <output file name>"
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


output=open(sys.argv[4], "w")



for location in locations:
	location=location.strip()
	if int(location) in noncore:
		print >>output, location, "noncore"
	else:
		print >>output, location, "core"
		
	

output.close()