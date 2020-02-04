#!/usr/bin/env python
#/usr/bin/env python
import string, re
import os, sys, getopt, random, math
import time
import pysam

MAX_INSERT=700

#CIGAR operators: 0=match, 1=insertion, 2=deletion, 3=skipped region from reference, 4=soft clip on read (i.e. sequence still in SEQ, 5=hard clip on read (sequence not in SEQ), 6=Padding (silent deletion from padded reference - for multiple alignment)
def calculate_cigar_length(cigar_sequence):
	length=0
	
	for f in cigar_sequence:
		if f[0] in [0,1,3,4,5]:
			length+=f[1]
	
	return length


#open the two SAM/BAM format input files (must be sorted by name)

if sys.argv[1].split(".")[-1]=="bam":
	samfile1 = pysam.Samfile( sys.argv[1], "rb" )
	#pysam.sort( "-n", sys.argv[1], "sort1" )
elif sys.argv[1].split(".")[-1]=="sam":
	samfile1 = pysam.Samfile( sys.argv[1], "r" )

if sys.argv[2].split(".")[-1]=="bam":
	samfile2 = pysam.Samfile( sys.argv[2], "rb" )
	#pysam.sort( "-n", sys.argv[2], "sort2" )
elif sys.argv[2].split(".")[-1]=="sam":
	samfile2 = pysam.Samfile( sys.argv[2], "r" )

samout = pysam.Samfile("test.bam", "wb", template=samfile1)


#print samfile1.header
#print samfile2.header


perfect_match_limits={}
smallindels={}

for x in range(len(samfile1.references)):
	perfect_match_limits[x]={}
	smallindels[x]={'insertions':{}, 'deletions': {}}

#reads_iter1=samfile1
reads_iter1=samfile1

reads_iter2=samfile2


potential_joins={}

for count, read1 in enumerate(reads_iter1):
	
	read2 = reads_iter2.next()
	read1mate = reads_iter1.next()
	read2mate = reads_iter2.next()
	
	
	if read1.is_reverse:
		readtmp=read1mate
		read1mate=read1
		read1=readtmp
		readtmp=read2mate
		read2mate=read2
		read2=readtmp
	
	
#
#	
	###########################################################
	# Deal with reads that map perfectly and the same in both #
	###########################################################

	if  read1.is_proper_pair and read2.is_proper_pair and read1.isize==read2.isize:
	
#		for n in range (read1.pos, read1.mpos):#+read1.rlen):
#			
#			if not refseqs[samfile1.references[read1.rname]][n].has_key(read2.rname):
#				refseqs[samfile1.references[read1.rname]][n][read2.rname]=0
#			refseqs[samfile1.references[read1.rname]][n][read2.rname]+=1
		
		added=False
		
		if not perfect_match_limits[read1.rname].has_key(read2.rname):
			perfect_match_limits[read1.rname][read2.rname]=[[-1,-1]]

		
		blocks=perfect_match_limits[read1.rname][read2.rname]
		#print blocks
		newblocks=[]
		newblockstart=-1
		newblockend=-1
		
		#Find the end of a mapped pair by adding the length of the mapped read to its position
		
		end = read1.mpos
		for f in read2mate.cigar:
			if f[0]!=1:
				end += f[1]
				
		#This bit fixes the bwa bug that adds mapping past the end of a contig! Can be removed with fixed bwa		
		
		if end>=samfile1.lengths[read1.rname]:
			end=samfile1.lengths[read1.rname]-1
		
		for block in blocks:
		
			#if it's the first block, add it
		
			if block[0]==-1:
				newblocks.append([read1.pos,end])
				added=True
		
			#if we've already added our read
			
			elif added==True:
				newblocks.append(block)
				
		
			#if we've already found the start of the new block and the mate position is less than the next previously built block, print both blocks to the new list
		
			elif newblockend!=-1 and block[0]>newblockend:
				newblocks.append([newblockstart,newblockend])
				added=True
				newblocks.append(block)
				newblockstart=-1
				newblockend=-1
			
			#if we've haven't found the start of the new block and the mate position is less than the next previously built block, print both blocks to the new list
			
			elif block[0]>end:
				newblocks.append([read1.pos,end])
				added=True
				newblocks.append(block)
			
			#if we've already found the start of the new block and the end is greater than the next previously built block, continue
			
			elif newblockend>block[1]:
				continue
		
			#if the start of the read is within a previous block and the mate is outside
		
			elif read1.pos>=block[0] and read1.pos<=block[1] and end>block[1]:
				#set the start of the new block to the start of the previous one, and the end to the mate position
				newblockstart=block[0]
				newblockend=end
			
			#if the mate is within a previous block and the read outside
			
			elif read1.pos<block[0] and end<=block[1] and end>=block[0]:
				#if we've not already found the start of the new block
				if newblockstart==-1:
					newblocks.append([read1.pos,block[1]])
					added=True
				#if we have
				else:
					newblocks.append([newblockstart,block[1]])
					newblockstart=-1
					added=True
			
			#if the reads span an entire block, set the start and end of the new block to the read and mate positions respectively
			
			elif read1.pos<block[0] and end>block[1]:
				newblockstart=read1.pos
				newblockend=end
			
			
			#otherwise, add the old block to the list of new blocks
			elif read1.pos>=block[0] and end<=block[1]:
				newblocks.append(block)
				added=True
			
			#otherwise, add the old block to the list of new blocks
			else:
				newblocks.append(block)

		#if the reads haven't been added to the new list of blocks, add them now
		if added==False:
			newblocks.append([read1.pos,end])
		#print newblocks
		perfect_match_limits[read1.rname][read2.rname]=newblocks
	
	
	
	
	#############################################################
	# Identify reads at eands of contigs that are circularising #
	#############################################################
	
	
	
	elif not read1.is_proper_pair and not read1.is_unmapped and not read1.mate_is_unmapped and read1.rname==read1.mrnm and samfile1.lengths[read1.rname]>2*MAX_INSERT:
		if not read1.is_reverse and read1.pos+MAX_INSERT>samfile1.lengths[read1.rname] and read1.mate_is_reverse and read1.mpos-MAX_INSERT<0:
			print samfile1.references[read1.rname], read1.pos, read1.mpos
			
		if not read1.mate_is_reverse and read1.mpos+MAX_INSERT>samfile1.lengths[read1.rname] and read1.is_reverse and read1.pos-MAX_INSERT<0:
			print samfile1.references[read1.rname], read1.pos, read1.mpos
	
	
	elif not read2.is_proper_pair and not read2.is_unmapped and not read2.mate_is_unmapped and read2.rname==read2.mrnm and samfile2.lengths[read2.rname]>2*MAX_INSERT:
		if not read2.is_reverse and read2.pos+MAX_INSERT>samfile2.lengths[read2.rname] and read2.mate_is_reverse and read2.mpos-MAX_INSERT<0:
			print samfile2.references[read2.rname], read2.pos, read2.mpos
			
		if not read2.mate_is_reverse and read2.mpos+MAX_INSERT>samfile2.lengths[read2.rname] and read2.is_reverse and read2.pos-MAX_INSERT<0:
			print samfile2.references[read2.rname], read2.pos, read2.mpos
	
	
	
	############################################################################################
	# Identify where reads mapping in both are not the same distance apart (i.e. small indels) #
	############################################################################################
	
	
	elif read1.is_proper_pair and read2.is_proper_pair:
		if read1.rname==0:
#			print read1.pos, read2.pos, read1.isize, read2.isize
#			if read1.cigar!=read2.cigar:
#				print read1.cigar, read2.cigar
#			if read1mate.cigar!=read2mate.cigar:
#				print read1mate.cigar, read2mate.cigar
	
			insertions_in_read=0
			deletions_in_read=0
			baseposn=read1.pos
			for f in read1.cigar:
				#print f
				if f[0]==1:
					#print "here"
					insertions_in_read += f[1]
					for x in range(baseposn,baseposn+f[1]):
						if not smallindels[read1.rname]['insertions'].has_key(x):
							smallindels[read1.rname]['insertions'][x]=0
						smallindels[read1.rname]['insertions'][x]+=1
				elif f[0]==2:
					deletions_in_read += f[1]
					for x in range(baseposn,baseposn+f[1]):
						if not smallindels[read1.rname]['deletions'].has_key(x):
							smallindels[read1.rname]['deletions'][x]=0
						smallindels[read1.rname]['deletions'][x]+=1
				baseposn+=f[1]
		
			#print insertions_in_read, deletions_in_read
	
	
	##################################
	# Identify reads joining contigs #
	##################################
				
	elif read1.is_proper_pair and read2.is_read1 and not read2.is_unmapped and not read2.mate_is_unmapped:# and read2.rname!=read2.mrnm:
		if ((not read2.is_reverse and (samfile2.lengths[read2.rname]-(read2.pos+1))<1000) or (read2.is_reverse and (read2.pos+1)<1000)) and ((not read2.mate_is_reverse and (samfile2.lengths[read2.mrnm]-(read2.mpos+1))<1000) or (read2.mate_is_reverse and (read2.mpos+1)<1000)):
			
#			if not read2.is_reverse:
#				print read1.qname, "may join contigs", read2.rname, "and", read2.mrnm
#			else:
#				print read1.qname, "may join contigs", read2.mrnm, "and", read2.rname
			
			if read2.rname<read2.mrnm:
				if not potential_joins.has_key((read2.rname,read2.mrnm)):
					potential_joins[(read2.rname,read2.mrnm)]=0
				potential_joins[(read2.rname,read2.mrnm)]+=1
			else:
				if not potential_joins.has_key((read2.mrnm,read2.rname)):
					potential_joins[(read2.mrnm,read2.rname)]=0
				potential_joins[(read2.mrnm,read2.rname)]+=1
			
	else:
		samout.write(read1)
		samout.write(read1mate)
		
	
#	if count > 1000000:
#		break	

print smallindels[0]
#sys.exit()
	
outfile=open("test.tab","w")	
donelength=1
for refnum, ref in enumerate(samfile1.references):
	
	for contig in perfect_match_limits[refnum].keys():
	
		for block in perfect_match_limits[refnum][contig]:
		
			print >> outfile, "FT   contig          "+str(block[0]+donelength)+".."+str(block[1]+donelength)
			print >> outfile, "FT                   /primary_name="+samfile2.references[contig]
	
	donelength += samfile1.lengths[refnum]

for join in potential_joins.keys():
	print samfile2.references[join[0]], samfile2.references[join[1]], potential_joins[join]

