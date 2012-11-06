#!/usr/bin/env python
import os, sys, string
from numpy import mean
from random import sample
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from optparse import OptionParser
#from Bio.SeqIO import PairedFastaQualIterator

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-g", "--genome_size", action="store", dest="genome_size", help="Expected size of genome", default=3000000, type="int")
	parser.add_option("-c", "--coverage", action="store", dest="coverage", help="Target coverage", default=30, type="int")
	parser.add_option("-f", "--forward_fastq", action="store", dest="ffastq", help="Input forward fastq file name", default="")
	parser.add_option("-r", "--reverse_fastq", action="store", dest="rfastq", help="Input reverse fastq file name", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file name prefix", default="Subsample")

	
	
	return parser.parse_args()



def check_input_options(options, args):
	
	if not os.path.isfile(options.ffastq):
		print "Cannot find input file:", options.ffastq
	if options.rfastq!="" and not os.path.isfile(options.rfastq):
		print "Cannot find input file:", options.rfastq

(options, args) = main()
check_input_options(options, args)


try:
	reads=SeqIO.parse(open(options.ffastq), "fastq")
except StandardError:
	print "Could not open file", options.ffastq
	sys.exit()


quals=[]
if options.rfastq!="":
	target_length=(options.genome_size*(float(options.coverage)/2))
	print "Aiming for a total length of", target_length, "bases per fastq"
else:
	target_length=(options.genome_size*options.coverage)
	print "Aiming for a total length of", target_length, "bases"
	

for read in reads:
	readlen=len(read.seq)
	break


try:
	reads=SeqIO.parse(open(options.ffastq), "fastq")
except StandardError:
	print "Could not open file", options.ffastq
	sys.exit()

for x, read in enumerate(reads):
	quals.append([mean(read.letter_annotations["phred_quality"]), x])
#	if x>100000:
#		break
#	for nuc, qual in zip(read,read.letter_annotations["phred_quality"]):
#		print nuc, qual

if options.rfastq!="":
	
	try:
		reads=SeqIO.parse(open(options.rfastq), "fastq")
	except StandardError:
		print "Could not open file", options.rfastq
		
	for x, read in enumerate(reads):
		quals[x][0]=(quals[x][0]+mean(read.letter_annotations["phred_quality"]))/2
#		if x>100000:
#			break


quals.sort(reverse=True)
#quals.reverse()

#print quals[:10], quals[-10:]

#print readlen
total=0
minqual=quals[0]
got_cov=False
readnum=0
num_reads_wanted=len(quals)
for x, qual in enumerate(quals):
	if not got_cov and total>=target_length:
		got_cov=True
		num_reads_wanted=x
	if got_cov and qual[0]<minqual:
		break
	total+=readlen
	minqual=qual[0]
num_reads_found=x+1

tokeep=[]
ties=[]
for y in xrange(0,num_reads_found):
	if quals[y][0]>minqual:
		tokeep.append(quals[y][1])
	elif quals[y][0]==minqual:
		ties.append(quals[y][1])
quals=[]
#print len(tokeep), num_reads_found, num_reads_wanted
if (num_reads_wanted-len(tokeep))>=ties:
	tokeep=tokeep+ties
else:
	tie_sample=sample(ties,num_reads_wanted-len(tokeep))
	tokeep=tokeep+tie_sample

#print len(tie_sample), len(tokeep), len(tokeep)*readlen
	
print "Minimum quality found is", minqual

tokeep.sort()


try:
	reads=SeqIO.parse(open(options.ffastq), "fastq")
except StandardError:
	print "Could not open file", options.ffastq
	sys.exit()

curread=0
foutput=open(options.output+"_1.fastq", "w")
for x, read in enumerate(reads):
	if curread>=len(tokeep):
		break
	if x==tokeep[curread]:
		curread+=1
		print >> foutput, read.format("fastq")

foutput.close()

if options.rfastq!="":
	try:
		reads=SeqIO.parse(open(options.rfastq), "fastq")
	except StandardError:
		print "Could not open file", options.rfastq
		sys.exit()
	
	curread=0
	seqs=[]
	routput=open(options.output+"_2.fastq", "w")
	for x, read in enumerate(reads):
		if curread>=len(tokeep):
			break
		if x==tokeep[curread]:
			curread+=1
			print >> routput, read.format("fastq")
	
	routput.close()
