#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math

for x in sys.argv[1:]:
	print x
	lines=open(x, 'rU').readlines()
	output=open(x+'.new', 'w')
	count=1
	for line in lines:
		if count==2:
			print >> output, line.strip()[:36]
		elif count==4:
			print >> output, line.strip()[:36]
			count=0
		else:
			print >> output, line.strip()
		
			
		count=count+1
	output.close()
	