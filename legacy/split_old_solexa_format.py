#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time

def splitfile(readfile, run):

	strains={}

	total=0
	for x in readfile:
		lines=[x.strip()]
		
		for y in range(0,3):
			lines.append(readfile.next().strip())
#		if len(lines[1].split("#"))<2:
#			continue
		
		tagseq=lines[1].strip()[-6:]
		lines[1]=lines[1].strip()[:-6]
		lines[3]=lines[3].strip()[:-6]

		
		if tags.has_key(tagseq):
			tag=tags[tagseq]
		else:
			tag="0"
		
		#tag=lines[0].split("#")[1].split('/')[0]
		
		if not strains.has_key(tag):
			if sys.argv[3].lower()=='p':
				if stupid=='y':
					strains[tag]=open(run+"_"+fr+"_"+tag+".fastq","w")
				else:
					strains[tag]=open(run+"_"+tag+"_"+fr+".fastq","w")
					
			else:
				strains[tag]=open(run+"_"+tag+".fastq","w")
		
		for line in lines:
			print >> strains[tag], line
	
	for strain in strains.keys():		
		strains[strain].close()
	print "Done."#, total+count,'sequences extracted'



if ((len (sys.argv)!=4) and (len(sys.argv)!=5)) or '-h' in sys.argv[1:]:
	print "split_solexa.py <file containing tags> <fastq file to split> <single/paired (s/p)>"
	print "Tag file should be comma separated with the tag sequence followed by the tag name or number"
	print "E.g.:"
	print "ATCACG,HGSA142"
	print "CGATGT,3HK"
	print "TTAGGC,CHI59"
	print "or"
	print "ATCACG,1"
	print "CGATGT,2"
	print "TTAGGC,3"
	sys.exit()


lines=open(sys.argv[1], "rU").readlines()

tags={}
for line in lines:
	words=line.strip().split(",")
	if len(words)<2:
		continue
	tags[words[0]]=str(words[1])

if len(tags)<1:
	print "Haven't found any valid tags. Exiting"
	sys.exit()


start=time.time()
fr=sys.argv[2].split('_')[-1].split('.')[0]

if sys.argv[3].lower() not in ['s','p']:
	print "second argument must be s or p"
	sys.exit()
if len(sys.argv)==5 and sys.argv[4].lower()=='s':
	stupid='y'
else:
	stupid='n'
	

readfile=open(sys.argv[2], "rU")

#nlines=0
#print "Counting lines..."
##for line in readfile:
##	nlines += 1
#
#nlines=int(os.popen("wc -l "+sys.argv[1]).read().strip().split()[0])
#
#
#print nlines/4, "sequences to extract..."

print sys.argv[2]+":"
run="_".join(sys.argv[2].split('/')[-1].split('_')[:2])

splitfile(readfile,run)

if sys.argv[3].lower()=='p':
	splitname=sys.argv[2].split("_")
	if fr=='1':
		fr='2'
		splitname[-1]="2"+splitname[-1][1:]
		print "_".join(splitname)+":"
		readfile=open("_".join(splitname), "rU")
	else:
		fr='1'
		splitname[-1]="1"+splitname[-1][1:]
		print "_".join(splitname)+":"
		readfile=open("_".join(splitname), "rU")
	splitfile(readfile,run)
		
print time.time()-start

