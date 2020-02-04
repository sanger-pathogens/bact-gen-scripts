#!/usr/bin/env python


import os, sys
from optparse import OptionParser

##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options]"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	

	parser.add_option("-q", "--query", action="store", dest="queries", help="multifasta containing query sequences used in blast", default="", metavar="FILE")
	parser.add_option("-b", "--blastout", action="store", dest="blastout", help="blast output file in m8 format", default="", metavar="FILE")
	

	return parser.parse_args()

(options, args)=get_user_options()

queryfile=options.queries

querydata=open(queryfile,"rU").read()
querylengths={}
for line in querydata.split(">")[1:]:
	bits=line.split("\n")
	querylengths[bits[0].split()[0]]=len(''.join(bits[1:]))


filename=options.blastout

lines=open(filename,"rU").readlines()

contig=""
taxon=""
data=[]

print "\t".join(["Contig", "Top hit", "gi", "accession", "%id of hits to top hit", "total length of hits to top hit"])

for line in lines:
	words=line.strip().split()
	if words[0]=="":
		continue
	if words[0]!=contig:
		if not words[0] in querylengths:	
			print "blast file and query file don't match"
			sys.exit()
		
		if contig!="":
			totallength=0
			totalid=0.0
			for datum in data:
				totallength=totallength+datum[1]
				totalid=totalid+(datum[0]*datum[1])
			totalid=totalid/totallength
			print "\t".join([contig, '"'+taxon+'"', gi, accession, str(totalid), str(totallength)])
		covered=[0]*querylengths[words[0]]
		percentid=[0.0]*querylengths[words[0]]
		accession=words[1].split("|")[3]
		gi=words[1].split("|")[1]
		taxon=os.popen("pfetch -D "+accession).read().strip()
		contig=words[0]
		for x in xrange(int(words[7])-1,int(words[8])):
			covered[x]=1
		data=[(float(words[2]),int(words[3]))]
	else:
		data.append((float(words[2]),int(words[3])))


if contig!="":
	totallength=0
	totalid=0.0 
	for datum in data:
		totallength=totallength+datum[1]
		totalid=totalid+(datum[0]*datum[1])
	totalid=totalid/totallength
	print contig, taxon, gi, accession, totalid, totallength

