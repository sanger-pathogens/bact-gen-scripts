#!/usr/bin/env python


import os, sys

filename=sys.argv[1]

lines=open(filename,"rU").readlines()

contig=""
taxon=""
data=[]

print "\t".join(["Contig", "Top hit", "%id of hits to top hit", "total length of hits to top hit"])

for line in lines:
	words=line.strip().split()
	if words[0]!=contig:
		if contig!="":
			totallength=0
			totalid=0.0
			for datum in data:
				totallength=totallength+datum[1]
				totalid=totalid+(datum[0]*datum[1])
			totalid=totalid/totallength
			print "\t".join([contig, '"'+taxon+'"', str(totalid), str(totallength)])
		accession=words[1].split("|")[3]
		taxon=os.popen("pfetch -D "+accession).read().strip()
		contig=words[0]
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
	print contig, taxon, totalid, totallength

