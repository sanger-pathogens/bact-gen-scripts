#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time

def splitfile(readfile, run):

	strains={}
#
	
	count=0
	total=0
	for x in readfile:
		#lines=[x.strip()]
		
		count=count+1
		if count>=10000:
			total=total+count
			count=0
			#print str(int((float(total)/nlines)*100)),'% extracted\r',
			print total,'sequences extracted\r',
			sys.stdout.flush()
		
		tag=x.split("#")[1].split('/')[0]
		
		if not strains.has_key(tag):
			strains[tag]=[]
		
		strains[tag].append(x.strip())
		
		for y in range(0,3):
			strains[tag].append(readfile.next().strip())
	
	for strain in strains.keys():		
		strains[strain].close()
	print "Done.", total+count,'sequences extracted'
	
	for strain in strains.keys():
		if badtags=="n" and strain=="0":
			continue
		if sys.argv[2].lower()=='p':
			output=open(run+"_"+str(strain)+"_"+fr+".fastq","w")	
		else:
			output=open(run+"_"+str(strain)+".fastq","w")
		
		for line in x:
			print >> output, line
		output.close()



if ((len (sys.argv)!=3) and (len(sys.argv)!=4)) or '-h' in sys.argv[1:]:
	print "split_solexa.py <fastq file to split> <single/paired (s/p)> <print badtag file (b)>"
	sys.exit()

start=time.time()
fr=sys.argv[1].split('_')[-1].split('.')[0]

if sys.argv[2].lower() not in ['s','p']:
	print "second argument must be s or p"
	sys.exit()
if len(sys.argv)==4 and sys.argv[3].lower()=='b':
	badtags='y'
else:
	badtags='n'
	

readfile=open(sys.argv[1], "rU")

#nlines=0
#print "Counting lines..."
##for line in readfile:
##	nlines += 1
#
#nlines=int(os.popen("wc -l "+sys.argv[1]).read().strip().split()[0])
#
#
#print nlines/4, "sequences to extract..."

print sys.argv[1]+":"
run="_".join(sys.argv[1].split('/')[-1].split('_')[:2])

splitfile(readfile,run)

if sys.argv[2].lower()=='p':
	splitname=sys.argv[1].split("_")
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

