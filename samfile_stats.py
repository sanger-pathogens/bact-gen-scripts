#!/usr/bin/env python
#/usr/bin/env python
import string
import os, sys
import pysam
from random import *

print "File", "Total", "Mapped", "%_Mapped", "In_proper_pair", "%_In_proper_pair", "Read_and_Mate_Mapped", "Mate_on_different_contig", "Singletons", "%_Singletons", "Unmapped", "%_Unmapped"
for samfilename in sys.argv[1:]:
	
	if samfilename.split(".")[-1]=="bam":
		try:
			samfile = pysam.Samfile( samfilename, "rb" )
		except StandardError:
			print "Cannot open",samfilename, "skipping"
			continue
	elif samfilename.split(".")[-1]=="sam":
		try:
			samfile = pysam.Samfile( samfilename, "r" )
		except StandardError:
			print "Cannot open",samfilename, "skipping"
			continue
	else:
		print samfilename, "Not a bam or sam file. Must end in .sam or .bam"
		continue
		
	refs=samfile.references
	lengths=samfile.lengths
	
	Total=0
	Properpair=0
	Mapped=0
	Mappedwithmate=0
	Unmapped=0
	Singletons=0
	Mateondiffcontig=0
	
	if len(refs)>1:
		Totalchr={}
		Properpairchr={}
		Mappedwithmatechr={}
		Mappedchr={}
		Mateondiffcontigchr={}
		
		Singletonschr={}
		
		for ref in refs:
			Totalchr[ref]=0
			Properpairchr[ref]=0
			Mappedwithmatechr[ref]=0
			Mappedchr[ref]=0
			Mateondiffcontigchr[ref]=0
			Singletonschr[ref]=0
		
	
	for read in samfile:
		Total+=1
		if not read.is_unmapped:
			Mapped+=1
			if len(refs)>1:
				Mappedchr[samfile.getrname(read.rname)]+=1
		
			if not read.mate_is_unmapped:
				Mappedwithmate+=1
				if len(refs)>1:
					Mappedwithmatechr[samfile.getrname(read.rname)]+=1
			else:
				Singletons+=1
				if len(refs)>1:
					Singletonschr[samfile.getrname(read.rname)]+=1
			
		else:
			Unmapped+=1
		
		if read.is_proper_pair:
			Properpair+=1
			if len(refs)>1:
				Properpairchr[samfile.getrname(read.rname)]+=1
		
		if len(refs)>1 and read.rname!=read.mrnm and not read.mate_is_unmapped and not read.is_unmapped:
			Mateondiffcontig+=1
			Mateondiffcontigchr[samfile.getrname(read.rname)]+=1
		
		
	print samfilename, Total, Mapped, 100*(float(Mapped)/Total), Properpair, 100*(float(Properpair)/Total), Mappedwithmate, Mateondiffcontig, Singletons, 100*(float(Singletons)/Total), Unmapped, 100*(float(Unmapped)/Total)
		
		
		
		