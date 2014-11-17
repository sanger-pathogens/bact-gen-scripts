#!/usr/bin/env python
import os, sys, string, numpy
from Bio.Seq import Seq
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from optparse import OptionParser

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-m", "--min_contig_len", action="store", dest="minlength", help="Minimum contig length to include in stats", default=0, type="int")
	
	
	return parser.parse_args()



def check_input_options(options, args):
	
	if options.minlength<0:
		print "Minimum length of contigs to include must be >0"



(options, args) = main()
check_input_options(options, args)

print '\t'.join(["File", "Total length", "No. contigs", "Mean length", "Stdev lengths", "Median length", "Max length", "Min length", "N50", "N50n", "MeanGC", "GC% stdev", "Median GC%", "Max GC%", "MinGC%"])

for filename in args:
	try:
		contigs=SeqIO.parse(open(filename), "fasta")
	except StandardError:
		print "Could not open file", filename
		continue
	
	
	
#	print filename
	
	lengths=[]
	GCs=[]
	test=[]
	
	for contig in contigs:
	
		
		seq=str(contig.seq)
		length=len(seq)
		
		if length<options.minlength:
			continue
		
		lengths.append(length)
		test.append([length, contig.name])
		
		lengthnons=len(seq.upper().replace("N",""))
		if lengthnons==0:
			GC=0
		else:
			GC=(float(len(seq.upper().replace("N","").replace("A","").replace("T","")))/lengthnons)*100
		
		GCs.append(GC)
	
	
	test.sort()
	test.reverse()
	print test[0]
		
	if len(lengths)==0:
		print '\t'.join([filename, '0', '0', '0', '0', '0', '0', '0', '0', '0', "-", "-", "-", "-", "-"])
		continue
		
	
#	print "Total length =", numpy.sum(lengths)
#	print "Number of contigs =", len(lengths)
#	print "Mean length =", numpy.mean(lengths)
#	print "Standard deviation of lengths =", numpy.std(lengths)
#	print "Maximum length =", numpy.max(lengths)
#	print "Minimum length =", numpy.min(lengths)
	lengths.sort()
	lengths.reverse()
	fifty=float(numpy.sum(lengths))/2
	count=0
	sum=0
	print lengths[0]
	
	while sum<fifty:
		N50=lengths[count]
		sum+=lengths[count]
		count+=1
	try:
		print '\t'.join(map(str,[filename, numpy.sum(lengths), len(lengths), numpy.mean(lengths), numpy.std(lengths), numpy.median(lengths), numpy.max(lengths), numpy.min(lengths), N50, count, numpy.mean(GCs), numpy.std(GCs), numpy.median(GCs), numpy.max(GCs), numpy.min(GCs)]))
	except StandardError:
		print filename+"\tfailed"
#	print "N50 =", N50
#	print "N50n =", count
#	print "Mean GC% =", numpy.mean(GCs)
#	print "GC% standard deviation =", numpy.mean(GCs)
#	print "Median GC% =", numpy.mean(GCs)
#	print "Maximum GC% =", numpy.max(GCs)
#	print "Minimum GC% =", numpy.min(GCs)
	
	
	
