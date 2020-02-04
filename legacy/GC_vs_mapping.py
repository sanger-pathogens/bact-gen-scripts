#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math


if (len (sys.argv)!=3) or '-h' in sys.argv[1:]:
	print "GC_vs_mapping.py <mapping file> <GC plot file>"
	sys.exit()

print sys.argv[1].split("/")[-1].split(".")[0]



if sys.argv[1].split('.')[-1]=='gz':
	mapping=gzip.open(sys.argv[1],'r').read().split('\n')[:-1]
else:
	mapping=open(sys.argv[1],'r').read().split('\n')[:-1]
	
GC=open(sys.argv[2],'rU').read().split('\n')[:-1]


GCmap={}

for x in range(0,len(mapping)):
	
	if not GCmap.has_key(int(GC[x])):
		GCmap[int(GC[x])]=[int(mapping[x]),1]
	else:
		GCmap[int(GC[x])][0]=GCmap[int(GC[x])][0]+int(mapping[x])
		GCmap[int(GC[x])][1]=GCmap[int(GC[x])][1]+1
		

keys=GCmap.keys()
keys.sort()

for key in keys:
	print key,  float(GCmap[key][0])/GCmap[key][1]
	
