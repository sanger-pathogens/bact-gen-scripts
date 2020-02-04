#!/usr/bin/env python


##################
# Import modules #
##################

import string, re, os, sys

from optparse import OptionParser, OptionGroup
from random import *
from numpy import mean, max, min, median, std, sum
import pysam


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


try: samfile = pysam.Samfile( sys.argv[1], "rb" )
except StandardError:
	print tmpname+".bam not a bam file"
	sys.exit() 

refs=samfile.references
lengths=samfile.lengths

joins={}

to_align={}

for read in samfile:
	if not read.is_unmapped and not read.mate_is_unmapped and not read.is_proper_pair:
		if read.tid != read.rnext:
			
				
			insertsize=int(read.qname.split(":")[-1])
			
			if not read.is_reverse:
				read_to_end=lengths[read.tid]-read.pos
			else:
				read_to_end=read.pos+read.qend
			
			if not read.is_reverse:
				mate_to_end=lengths[read.rnext]-read.pnext
			else:
				mate_to_end=read.pnext+read.rlen
				
			#print insertsize, read_to_end, mate_to_end, read_to_end+mate_to_end
			
			if read_to_end+mate_to_end>insertsize:
				continue
			
			if refs[read.tid] in ["7648_1#3_shuffled_94_cov_37", "7648_1#3_shuffled_116_cov_42"] and refs[read.rnext] in ["7648_1#3_shuffled_94_cov_37", "7648_1#3_shuffled_116_cov_42"]:
				
				if refs[read.tid]=="7648_1#3_shuffled_94_cov_37":
					if not read.is_reverse:
						to_align[read.qname.split(":")[0]]="f"
					else:
						to_align[read.qname.split(":")[0]]="r"
				else:
					if read.mate_is_reverse:
						to_align[read.qname.split(":")[0]]="f"
					else:
						to_align[read.qname.split(":")[0]]="r"
			
			if refs[read.tid] in joins:
				if refs[read.rnext] in joins[refs[read.tid]]:
					joins[refs[read.tid]][refs[read.rnext]][0]+=1
					joins[refs[read.tid]][refs[read.rnext]][1].append(insertsize-(read_to_end+mate_to_end))
				elif refs[read.rnext] in joins and refs[read.tid] in joins[refs[read.rnext]]:
					joins[refs[read.rnext]][refs[read.tid]][0]+=1
					joins[refs[read.rnext]][refs[read.tid]][1].append(insertsize-(read_to_end+mate_to_end))
				else:
					joins[refs[read.tid]][refs[read.rnext]]=[1,[insertsize-(read_to_end+mate_to_end)]]
			elif refs[read.rnext] in joins:
				if refs[read.tid] in joins[refs[read.rnext]]:
					joins[refs[read.rnext]][refs[read.tid]][0]+=1
					joins[refs[read.rnext]][refs[read.tid]][1].append(insertsize-(read_to_end+mate_to_end))
				else:
					joins[refs[read.rnext]][refs[read.tid]]=[1,[insertsize-(read_to_end+mate_to_end)]]
			else:
				joins[refs[read.tid]]={}
				joins[refs[read.tid]][refs[read.rnext]]=[1,[insertsize-(read_to_end+mate_to_end)]]

rank=[]
for join in joins:
	for joinb in joins[join]:
		rank.append([joins[join][joinb][0], join, joinb, mean(joins[join][joinb][1]), std(joins[join][joinb][1]), median(joins[join][joinb][1]), max(joins[join][joinb][1]), min(joins[join][joinb][1])])
		#print join, joinb, joins[join][joinb]
rank.sort()
#rank.reverse()

for r in rank:
	print "\t".join(map(str,r))
	
outfile=open("tmp.fasta", "w")

readfile=open(sys.argv[2],"rU")

for line in readfile:
	newlines=[line.strip()]
	for y in range(0,3):
		newlines.append(readfile.next().strip())
	
	
	if newlines[0][1:] in to_align:
		print >> outfile, ">"+newlines[0][1:]
		if to_align[newlines[0][1:]]=="f":
			print >> outfile, newlines[1]
		else:
			print >> outfile, revcomp(newlines[1])
	


outfile.close()