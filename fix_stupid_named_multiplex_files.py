#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math


for f in sys.argv[1:]:
	
	os.system("mv "+f+" "+f+".old")
	

for f in sys.argv[1:]:
	print f,
	#get suffix (should be fastq or fastq.gz)
	
	suffix = ".".join(f.split(".")[1:])
	
	#print "suffix =", suffix
	
	#get prefix
	
	prefix = f.split(".")[0]
	
	#print "prefix =", prefix
	
	#split names by underscore
	
	name=prefix.split("_")
	
	#check file is multiplexed
	
	if len(name)!=4:
		print "File doesn't seem to be multiplexed!"
		os.system("mv "+f+".old "+f)
		continue
	
	try:
		tag=int(name[3])
	except StandardError:
		print "Tag doesn't seem to be an integer!"
		os.system("mv "+f+".old "+f)
		continue
		
	
	try:
		fr=int(name[2])
	except StandardError:
		print "Forward/reverse key doesn't seem to be an integer!"
		os.system("mv "+f+".old "+f)
		continue
	
	if fr>2:
		print "Forward/reverse must be 1 or 2!"
		os.system("mv "+f+".old "+f)
		continue
	
	if fr==tag:
		print "->", f
		os.system("mv "+f+".old "+f)
		continue
	
	#move file	
	print "->", ".".join(["_".join([name[0], name[1], name[3], name[2]]), suffix])
	os.system("mv "+f+".old "+".".join(["_".join([name[0], name[1], name[3], name[2]]), suffix]))
	#print ".".join(["_".join([name[0], name[1], name[3], name[2]]), suffix])
	


	
		
	






