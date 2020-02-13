#!/usr/bin/env python
import os, sys, string, numpy, scipy.stats
from Bio.Seq import Seq
from optparse import OptionParser
from modules.Si_SeqIO import *

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-m", "--min_contig_len", action="store", dest="minlength", help="Minimum contig length to include in stats", default=0, type="int")
	parser.add_option("-s", "--split_contigs", action="store_true", dest="split", help="Print stats for each contig as well as each file", default=False)
	
	
	return parser.parse_args()



def check_input_options(options, args):
	
	if options.minlength<0:
		print "Minimum length of contigs to include must be >0"



(options, args) = main()
check_input_options(options, args)

print '\t'.join(["File", "Total length", "No. contigs", "Mean length", "Stdev lengths", "Median length", "Max length", "Min length", "Skewness", "Kurtosis", "N50", "N50n", "Mean GC%", "Stdev GC%", "Median GC%", "Max GC%", "MinGC%", "N count", "gap count", "Length without Ns and gaps"])

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
	gaps=[]
	Ns=[]
	lenacgt=[]
	
	for contig in contigs:
		seq=str(contig.seq)
		length=len(seq)
		
		if length<options.minlength:
			continue
		
		lengths.append(length)
		
		lengthnons=len(seq.upper().replace("N",""))
		lengthnogaps=len(seq.upper().replace("-",""))
		lengthnogapnons=len(seq.upper().replace("-","").replace("N",""))
		if lengthnons==0:
			GC=0
		else:
			GC=(float(len(seq.upper().replace("N","").replace("A","").replace("T","")))/lengthnogapnons)*100
		
		Ns.append(length-lengthnons)
		gaps.append(length-lengthnogaps)
		lenacgt.append(lengthnogapnons)
		GCs.append(GC)
		
		if options.split:
			try:
				print '\t'.join(map(str,[contig.name, length, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", GC, "-", "-", "-", "-", Ns[-1], gaps[-1], lenacgt[-1]]))
			except StandardError:
				print filename+"\tfailed"
	
	
	#test.sort()
	#test.reverse()
	#print test[0]
		
	if len(lengths)==0:
		print '\t'.join([filename, '0', '0', '0', '0', '0', '0', '0', '0', '0', '-', '-', "-", "-", "-", "-", "-", "-", "-", "-"])
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
	suml=0
	#print lengths[0]
	
	while suml<fifty:
		N50=lengths[count]
		suml+=lengths[count]
		count+=1
	try:
		print '\t'.join(map(str,[filename, numpy.sum(lengths), len(lengths), numpy.mean(lengths), numpy.std(lengths), numpy.median(lengths), numpy.max(lengths), numpy.min(lengths), scipy.stats.skew(lengths), scipy.stats.kurtosis(lengths), N50, count, numpy.mean(GCs), numpy.std(GCs), numpy.median(GCs), numpy.max(GCs), numpy.min(GCs), sum(Ns), sum(gaps), sum(lenacgt)]))
	except StandardError:
		print filename+"\tfailed"
#	print "N50 =", N50
#	print "N50n =", count
#	print "Mean GC% =", numpy.mean(GCs)
#	print "GC% standard deviation =", numpy.mean(GCs)
#	print "Median GC% =", numpy.mean(GCs)
#	print "Maximum GC% =", numpy.max(GCs)
#	print "Minimum GC% =", numpy.min(GCs)

	
	
	
