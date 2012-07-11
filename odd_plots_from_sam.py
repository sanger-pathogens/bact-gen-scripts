#!/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time


def setBit( x, bitNum ):
  return x | 1 << bitNum
  
def clearBit( x, bitNum ):
  return x & ~( 1 << bitNum )
  
def checkBit( x, bitNum ):
  return (x & (1<<bitNum)) != 0



#function to convert extended cigar string to a list
def readExtendedCIGAR(cigarstring):

	cigarlist=[]
	prevfind=0
	
	for x in range(0, len(cigarstring)):
		if cigarstring[x] in ["M", "I", "D", "N", "S", "H", "P"]:
			cigarlist.append([cigarstring[x], int(cigarstring[prevfind:x])])
			prevfind=x+1
	

	return cigarlist


#function to add extended cigar indels to a plot list

def addinsertionstoplotlist(cigarlist, contig, startingbase, plotlist):

	position=startingbase-1

	for x in cigarlist:
		if x[0]=="M":
			for y in range(0, x[1]):
				position=position+1
		elif x[0]=="D":
			for y in range(0, x[1]):
				position=position+1
		
		elif x[0]=="I":
			if position<len(plotlist[contig]):
				plotlist[contig][position]=plotlist[contig][position]+1

	
	return plotlist
	

#function to add extended cigar indels to a plot list

def adddeletionstoplotlist(cigarlist, contig, startingbase, plotlist):

	position=startingbase-1

	for x in cigarlist:
		if x[0]=="M":
			for y in range(0, x[1]):
				position=position+1
		elif x[0]=="D":
			for y in range(0, x[1]):
				if position<len(plotlist[contig]):
					plotlist[contig][position]=plotlist[contig][position]+1
				position=position+1

	
	return plotlist



#function to add extended cigar list to a plot list

def addcigarlisttoplotlist(cigarlist, contig, startingbase, plotlist):

	position=startingbase-1

	for x in cigarlist:
		if x[0]=="M":
			for y in range(0, x[1]):
				if position<len(plotlist[contig]):
					plotlist[contig][position]=plotlist[contig][position]+1
				position=position+1
		elif x[0]=="D":
			for y in range(0, x[1]):
				position=position+1
	
	return plotlist
	
	
#function to add inserts between reads that are too far apart to a plot list

def addlargeinsertstoplotlist(cigarlist, contig, startingbase, endbase, plotlist):

	position=startingbase-1

	for x in cigarlist:
		if x[0]=="M":
			for y in range(0, x[1]):
				position=position+1
		elif x[0]=="D":
			for y in range(0, x[1]):
				position=position+1
	
	for x in range(position,endbase):
		if x<len(plotlist[contig]):
			plotlist[contig][x]=plotlist[contig][x]+1
	
	return plotlist

#function to print a plot list to file			

def print_plot(contiglist, plotlist, outputfilename):
	output=open(outputfilename, "w")
	for contig in contiglist:
		for x in plotlist[contig]:
			print >> output, x
	output.close()
	
	
#function to print a list of plot list to file			

def print_multiple_plots(contiglist, plotlistlist, outputfilename):
	output=open(outputfilename, "w")
	for contig in contiglist:
		for x in range(0,len(plotlistlist[0][contig])):
			for plotlist in plotlistlist:
				print >> output, plotlist[contig][x],
			print >> output, "\n",
	output.close()


#bothneg=[0]*reflen
#bothpos=[0]*reflen
onlyforward={}
onlyreverse={}
complete_circle={}
Mapping={}
Smallinsertions={}
Smalldeletions={}
Reads_too_far_apart={}
inversions={}
#reverse=[0]*reflen


filex=open(sys.argv[1], "rU")
#outputx=open(sys.argv[1].replace(".sam","_straddle.plot"),"w")

outfile=sys.argv[2]

paired=True
if len(sys.argv)==4 and sys.argv[3]=="s":
	paired=False

contigsx={}
contigorder=[]

for x in filex:
	x=x.strip()
	while x[0]=='@':
		if x.split()[0]=="@SQ":
			contigsx[x.split()[1].split(":")[1]]=int(x.split()[2].split(":")[1])
			contigorder.append(x.split()[1].split(":")[1])
			
			onlyreverse[x.split()[1].split(":")[1]]=[0]*int(x.split()[2].split(":")[1])
			onlyforward[x.split()[1].split(":")[1]]=[0]*int(x.split()[2].split(":")[1])
			complete_circle[x.split()[1].split(":")[1]]=[0]*int(x.split()[2].split(":")[1])
			Mapping[x.split()[1].split(":")[1]]=[0]*int(x.split()[2].split(":")[1])
			Smallinsertions[x.split()[1].split(":")[1]]=[0]*int(x.split()[2].split(":")[1])
			Smalldeletions[x.split()[1].split(":")[1]]=[0]*int(x.split()[2].split(":")[1])
			Reads_too_far_apart[x.split()[1].split(":")[1]]=[0]*int(x.split()[2].split(":")[1])
			inversions[x.split()[1].split(":")[1]]=[0]*int(x.split()[2].split(":")[1])
			
		#read the next line
		x=filex.next().strip()
	
	#add forward and reverse pairs to a list
	linesx=[x]
	if paired:
		linesx.append(filex.next().strip())
		extendedcigars=[readExtendedCIGAR(linesx[0].split()[5]), readExtendedCIGAR(linesx[1].split()[5])]
	else:
		extendedcigars=[readExtendedCIGAR(linesx[0].split()[5])]
	
	#print linesx[0].split()

	#add both read and mate to mapping plot if they both map and if they include an insert/deletion, add to small insert/deletions file
	if not checkBit(int(linesx[0].split()[1]), 2):
		Mapping=addcigarlisttoplotlist(extendedcigars[0], linesx[0].split()[2], int(linesx[0].split()[3]), Mapping)
		
		typelist=[]
		for x in extendedcigars[0]:
			typelist.append(x[0])
			
		if "I" in typelist:
			Smallinsertions=addinsertionstoplotlist(extendedcigars[0], linesx[0].split()[2], int(linesx[0].split()[3]), Smallinsertions)
		if "D" in typelist:
			Smalldeletions=adddeletionstoplotlist(extendedcigars[0], linesx[0].split()[2], int(linesx[0].split()[3]), Smalldeletions)
	
	if paired:	
		if not checkBit(int(linesx[1].split()[1]), 2):
			Mapping=addcigarlisttoplotlist(extendedcigars[1], linesx[1].split()[2], int(linesx[1].split()[3]), Mapping)
			typelist=[]
			for x in extendedcigars[1]:
				typelist.append(x[0])
				
			if "I" in typelist:
				Smallinsertions=addinsertionstoplotlist(extendedcigars[1], linesx[1].split()[2], int(linesx[1].split()[3]), Smallinsertions)
			if "D" in typelist:
				Smalldeletions=adddeletionstoplotlist(extendedcigars[1], linesx[1].split()[2], int(linesx[1].split()[3]), Smalldeletions)
				
			
	
	
	#if reads do not map in a proper pair
		if not checkBit(int(linesx[0].split()[1]), 1):
			if checkBit(int(linesx[0].split()[1]), 2):
				#print "read is not mapped"
				continue
			#else:
				#print linesx[0].split()[0], linesx[0].split()[2], linesx[0].split()[3], linesx[0].split()[8], linesx[0].split()[5],
				
			
			#if read maps, but mate does not
			if (checkBit(int(linesx[0].split()[1]), 3) and not checkBit(int(linesx[0].split()[1]), 2)):
				#print "mate is unmapped... potential insertion in query",
				#if read is on reverse strand
				if checkBit(int(linesx[0].split()[1]), 4):
					#print "read is on the reverse strand"
					onlyreverse=addcigarlisttoplotlist(extendedcigars[0], linesx[0].split()[2], int(linesx[0].split()[3]), onlyreverse)
				#else if read is on forward strand
				else:
					#print "read is on the forward strand"
					onlyforward=addcigarlisttoplotlist(extendedcigars[0], linesx[0].split()[2], int(linesx[0].split()[3]), onlyforward)
					
			elif (checkBit(int(linesx[0].split()[1]), 2) and not checkBit(int(linesx[0].split()[1]), 3)):
				#print "read is unmapped... potential insertion in query",
				#if read is on reverse strand
				if checkBit(int(linesx[1].split()[1]), 4):
					#print "mate is on the reverse strand"
					onlyreverse=addcigarlisttoplotlist(extendedcigars[1], linesx[1].split()[2], int(linesx[1].split()[3]), onlyreverse)
				#else if read is on forward strand
				else:
					#print "mate is on the forward strand"
					onlyforward=addcigarlisttoplotlist(extendedcigars[1], linesx[1].split()[2], int(linesx[1].split()[3]), onlyforward)
					
			#else if both read and mate map, but not properly
			else:
				#if both reads are on the same strand
				if (checkBit(int(linesx[0].split()[1]), 4) and checkBit(int(linesx[0].split()[1]), 5)) or (not checkBit(int(linesx[0].split()[1]), 4) and not checkBit(int(linesx[0].split()[1]), 5)):
					print "Reads are on same strand... potential inversion in query"
					inversions=addcigarlisttoplotlist(extendedcigars[1], linesx[1].split()[2], int(linesx[1].split()[3]), inversions)
				else:
					if int(linesx[0].split()[8])>(contigsx[linesx[0].split()[2]]-1000):
						#print "Reads probably straddle contig end"
						complete_circle=addcigarlisttoplotlist(extendedcigars[0], linesx[0].split()[2], int(linesx[0].split()[3]), complete_circle)
						#print extendedcigars[1], linesx[1].split()[2], int(linesx[1].split()[3])
						complete_circle=addcigarlisttoplotlist(extendedcigars[1], linesx[1].split()[2], int(linesx[1].split()[3]), complete_circle)
					elif int(linesx[0].split()[8])<100000 and int(linesx[0].split()[8])>100:
						print "Insert too big... Possible deletion in query"
						Reads_too_far_apart=addlargeinsertstoplotlist(extendedcigars[0], linesx[0].split()[2], int(linesx[0].split()[3]), int(linesx[0].split()[7]), Reads_too_far_apart)
						
	#				elif int(linesx[0].split()[8])<0:
	#					print "Reads overlap"
	#				elif int(linesx[0].split()[8])<100:
	#					print "Insert too small"
	#				else:
	#					print "Possible relocation of DNA??"
			
#print plots
print_plot(contigorder, Mapping, outfile+"_Mapped.plot")
if paired:
	print_plot(contigorder, inversions, outfile+"_Inversions.plot")

	print_multiple_plots(contigorder, [onlyforward, onlyreverse, complete_circle], outfile+"_No_mates.plot")
	print_multiple_plots(contigorder, [Smallinsertions, Smalldeletions, Reads_too_far_apart], outfile+"_Indels.plot")


