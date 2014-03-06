#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################
import string, re
import os, sys
from optparse import OptionParser, OptionGroup
import pysam



################################
# Get the command line options #
################################


def main():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-b", "--bam", action="store", dest="bamfile", help="Input bam file", type="string", metavar="FILE", default="")
	
	parser.add_option("-q", "--base_qual_filter", action="store", dest="base_qual_filter", help="Base quality filter for bam file mapping plots [default= %default]", default=0, type="float")
	parser.add_option("-Q", "--mapping_qual_filter", action="store", dest="mapping_qual_filter", help="Mapping quality filter for bam file plots [default= %default]", default=0, type="float")
	
	return parser.parse_args()


def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp

EITHER_PERCENT_CUTOFF=80
HIGH_PERCENT_CUTOFF=10
MIN_READS=4
MAX_MISMATCHES_ALLOWED=2

################
# Main program #
################		

if __name__ == "__main__":

	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	
	filename=options.bamfile
	
	if not os.path.isfile(filename):
		print 'Cannot find file', filename
		sys.exit()
	
	if options.base_qual_filter<0 or options.base_qual_filter>100:
		print 'Base quality filter must be between 0 and 100'
		sys.exit()
	if options.mapping_qual_filter<0 or options.mapping_qual_filter>100:
		print 'Mapping quality filter must be between 0 and 100'
		sys.exit()
	

	print "Calculating coverage for", filename
	sys.stdout.flush()
	
	try:
		samfile = pysam.Samfile( filename, "rb" )
	except StandardError:
		print 'Failed to open '+filename+'. Is it in bam format?'
		sys.exit()
	
	refs=samfile.references
	lengths=samfile.lengths
	
	contiglocs={}
	totallength=0
	indelreads={}
	startcov={}
	endcov={}	
	starts={}
	ends={}
	
# This code identifies clipped reads and the position they are clipped at	
	for x, ref in enumerate(refs):
		contiglocs[ref]=totallength
		startcov[ref]=[]
		endcov[ref]=[]
		for y in xrange(lengths[x]):
			startcov[ref].append(0.0)
			endcov[ref].append(0.0)
		totallength+=lengths[x]
		indelreads[x]={"insertion":0, "deletion":0, "clipped":{"start":{}, "end":{}}}
	
	
	for read in samfile:
		
		if not read.is_unmapped:
			refpos=read.pos
			if read.cigar[0][0]==4 and read.cigar[-1][0]==4:
				#print "Both ends clipped", read.tid
				continue
			for x, cig in enumerate(read.cigar):
				if cig[0]==0:
					refpos+=cig[1]
				elif cig[0]==1:
					indelreads[read.tid]["insertion"]+=1
				elif cig[0]==2:
					indelreads[read.tid]["deletion"]+=1
					refpos+=cig[1]
				elif cig[0]==4:
					#indelreads[read.tid]["clipped"]+=1
					#if refpos<lengths[read.tid] and refpos!=0:
					if x==0:
						if cig[1]>10:
							print "start", read.tid, refpos-cig[1], cig[1], refpos, lengths[read.tid]#, read
							if refpos==0:
								rpos=1
							else:
								rpos=refpos
							startcov[refs[read.tid]][rpos-1]+=1
							starts[read.qname]=revcomp(read.seq[:cig[1]])
						if not refpos in indelreads[read.tid]["clipped"]["start"]:
							indelreads[read.tid]["clipped"]["start"][refpos]=0
							
#							if not read.is_proper_pair:
#								print  "start", refs[read.tid], refpos, lengths[read.tid], "unpaired"
						indelreads[read.tid]["clipped"]["start"][refpos]+=1
					elif x==len(read.cigar)-1:
						if cig[1]>10:
							print "end", refs[read.tid], refpos, lengths[read.tid]#, read
							endcov[refs[read.tid]][refpos-1]+=1
							ends[read.qname]=read.seq[-1*cig[1]:]
						if not refpos-1 in indelreads[read.tid]["clipped"]["end"]:
							indelreads[read.tid]["clipped"]["end"][refpos-1]=0
							
#							if not read.is_proper_pair:
#								print "end", refs[read.tid], refpos, lengths[read.tid], "unpaired"
						indelreads[read.tid]["clipped"]["end"][refpos-1]+=1
					else:
						print "middle", read.tid, refpos-cig[1], cig[1], refpos, lengths[read.tid]
						
					#refpos+=cig[1]
		
#	startfile=open("starts.fa", "w")	
#	for start in starts:
#		print >> startfile, ">"+start
#		print >> startfile, starts[start]
#	startfile.close()
#	
#	endfile=open("ends.fa", "w")	
#	for end in ends:
#		print >> endfile, ">"+end
#		print >> endfile, ends[end]
#	endfile.close()
	print startcov, endcov
	startplot=open("starts.plot", "w")
	print >> startplot, "#Base left_clipped right_clipped"
	for ref in refs:
		for x in xrange(len(startcov[ref])):
			print >> startplot, x+1, startcov[ref][x], endcov[ref][x]
	startplot.close()
	
	pairwise_start={}
	sequence_mismatches=[]
	
	for s1 in starts:
		pairwise_start[s1]={}
		seq_mismatches=0
		for s2 in starts:
			if s1==s2:
				pairwise_start[s1][s2]=0
				continue
			seq1=starts[s1]
			seq2=starts[s2]
			matches=0
			mismatches=0
			for y, base1 in enumerate(seq1):
				if y+1>len(seq2):
					continue
				base2=seq2[y]
				if base1==base2:
					matches+=1
				else:
					mismatches+=1
			pairwise_start[s1][s2]=mismatches
			if mismatches>MAX_MISMATCHES_ALLOWED:
				seq_mismatches+=1
				#print "start", s1, s2, matches, mismatches
		
		sequence_mismatches.append([seq_mismatches, s1])
	
	sequence_mismatches.sort()
	
	sclusters=[]
	
	for i, x in enumerate(sequence_mismatches):
		inc=False
		for c in sclusters:
			for y in c:
				inc=True
				if pairwise_start[x[1]][y]>MAX_MISMATCHES_ALLOWED:
					inc=False
			if inc:
				c.append(x[1])
				continue
		if not inc:
			sclusters.append([x[1]])
	
	for clust in sclusters:
		print "s", clust
	
	
	pairwise_end={}
	sequence_mismatches=[]
	for e1 in ends:
		pairwise_end[e1]={}
		seq_mismatches=0
		for e2 in ends:
			if e1==e2:
				pairwise_end[e1][e2]=0
				continue
			seq1=ends[e1]
			seq2=ends[e2]
			matches=0
			mismatches=0
			for y, base1 in enumerate(seq1):
				if y+1>len(seq2):
					continue
				base2=seq2[y]
				if base1==base2:
					matches+=1
				else:
					mismatches+=1
			pairwise_end[e1][e2]=mismatches
			if mismatches>MAX_MISMATCHES_ALLOWED:
				seq_mismatches+=1
				#print "end", e1, e2, matches, mismatches
		sequence_mismatches.append([seq_mismatches, e1])
	
	
	sequence_mismatches.sort()
#	print sequence_mismatches, ends
	eclusters=[]
	
	for i, x in enumerate(sequence_mismatches):
		inc=False
		for c in eclusters:
			for y in c:
				inc=True
				if pairwise_end[x[1]][y]>MAX_MISMATCHES_ALLOWED:
					inc=False
			if inc:
				c.append(x[1])
				continue
		if not inc:
			eclusters.append([x[1]])
	
	for clust in eclusters:
		print "e", clust
#	
#	for cluster in clusters:
#		for x in cluster:
#			for y in cluster:
#				print pairwise_start[x][y],
#			print
	
	
#	for x in ends:
#		for y in ends:
#			print pairwise_end[x][y],
#		print
	
	
#	for refnum in indelreads:
#		ref=refs[refnum]
#		reflen=lengths[refnum]
#		#print ref,  indelreads[refnum]
#		if len(indelreads[refnum]["clipped"]["start"])>0:
#			for pos in indelreads[refnum]["clipped"]["start"]:
#				if pos!=0 and indelreads[refnum]["clipped"]["start"][pos]>4:
#					#print pos, indelreads[refnum]["clipped"]["start"][pos], reflen
#					for pileupcolumn in samfile.pileup(ref, pos, pos+1, truncate=True):
#						if float(indelreads[refnum]["clipped"]["start"][pos])/pileupcolumn.n>0.1:
#							print ref, "start", reflen, pileupcolumn.pos, pileupcolumn.n, indelreads[refnum]["clipped"]["start"][pos], float(indelreads[refnum]["clipped"]["start"][pos])/pileupcolumn.n
#		
#		if len(indelreads[refnum]["clipped"]["end"])>0:
#			for pos in indelreads[refnum]["clipped"]["end"]:
#				#print pos, indelreads[refnum]["clipped"]["end"][pos], reflen
#				if pos+1!=reflen and indelreads[refnum]["clipped"]["end"][pos]>4:
#					for pileupcolumn in samfile.pileup(ref, pos, pos+1, truncate=True):
#						if float(indelreads[refnum]["clipped"]["end"][pos])/pileupcolumn.n>0.1:
#							print ref, "end", reflen, pileupcolumn.pos+1, pileupcolumn.n, indelreads[refnum]["clipped"]["end"][pos], float(indelreads[refnum]["clipped"]["end"][pos])/pileupcolumn.n
