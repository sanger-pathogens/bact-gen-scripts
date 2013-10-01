#!/usr/bin/env python

#################################
# Import some necessary modules #
#################################
import string, re, math
import os, sys
from optparse import OptionParser, OptionGroup
import pysam
from Bio.Seq import Seq
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pylab


################################
# Get the command line options #
################################


def main():
	usage = "usage: %prog [options] args"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-r", "--ref", action="store", dest="ref", help="Acessory reference fasta file", type="string", metavar="FILE", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="output prefix [default = %default]", type="string", metavar="FILE", default="")
	parser.add_option("-c", "--core", action="store", dest="corebam", help="Core bam file", type="string", metavar="FILE", default="")
	parser.add_option("-a", "--accessory", action="store", dest="accbam", help="Accessory bam file", type="string", metavar="FILE", default="")
	parser.add_option("-q", "--base_qual_filter", action="store", dest="base_qual_filter", help="Base quality filter for bam file mapping plots [default= %default]", default=0, type="float")
	parser.add_option("-Q", "--mapping_qual_filter", action="store", dest="mapping_qual_filter", help="Mapping quality filter for bam file plots [default= %default]", default=0, type="float")
	parser.add_option("-m", "--min_length", action="store", dest="minlen", help="Minimum contig length to keep", type="int", metavar="INT", default=1000)
	parser.add_option("-p", "--percentcutoff", action="store", dest="percent", help="Remove accessorycontigs with coverage < this percent of the core mean coverage [default = %default]", type="float", default=50)
	parser.add_option("-C", "--cliplen", action="store", dest="cliplen", help="Length to ignore at end of each contig (should be ~ one read length or more) [default = %default]", type="int", metavar="INT", default=100)
	
	return parser.parse_args()


################
# Main program #
################		

if __name__ == "__main__":

	#starttime=time.clock()

	#Get command line arguments

	(options, args) = main()
	
	
	if not os.path.isfile(options.corebam):
		print 'Cannot find file', options.corebam
		sys.exit()
	
	if not os.path.isfile(options.ref):
		print 'Cannot find file', options.ref
		sys.exit()
	
	if not os.path.isfile(options.accbam):
		print 'Cannot find file', options.accbam
		sys.exit()
	
	if options.output=="":
		#print 'No output file name given'
		#sys.exit()
		options.output="test"
	
	if options.base_qual_filter<0 or options.base_qual_filter>100:
		print 'Base quality filter must be between 0 and 100'
		sys.exit()
	if options.mapping_qual_filter<0 or options.mapping_qual_filter>100:
		print 'Mapping quality filter must be between 0 and 100'
		sys.exit()
	

	print "Calculating coverage for core file:", options.corebam
	sys.stdout.flush()
	
	try:
		coresamfile = pysam.Samfile( options.corebam, "rb" )
	except StandardError:
		print 'Failed to open '+options.corebam+'. Is it in bam format?'
		sys.exit()
	
	refs=coresamfile.references
	lengths=coresamfile.lengths
	
	contiglocs={}
	totallength=0
	
	for x, ref in enumerate(refs):
		contiglocs[ref]=totallength
		totallength+=lengths[x]
	
			
	
	
	poscount=0
	
	totcov=0.0
	for x, ref in enumerate(refs):
		lastcolumn=-1
		zerocount=0
		for pileupcolumn in coresamfile.pileup(ref):
			
			poscount+=1
			while pileupcolumn.pos!=lastcolumn+1:
				lastcolumn+=1
				zerocount+=1
			
			if pileupcolumn.n>0:
				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
					#print str(pileupcolumn)
					filtered_depth=0
					for pileupread in pileupcolumn.pileups:
						#print pileupread.alignment
						q=ord(pileupread.alignment.qual[pileupread.qpos])-33
						Q=pileupread.alignment.mapq
						#print q, Q, options.base_qual_filter, options.mapping_qual_filter
						if q>=options.base_qual_filter and Q>=options.mapping_qual_filter:
							filtered_depth+=1
					totcov+=filtered_depth
					#sys.exit()
				else:
					totcov+=pileupcolumn.n
			lastcolumn=pileupcolumn.pos

		
	coremean=totcov/totallength
	print coremean
	
	
	
	try:
		contigs=SeqIO.parse(open(options.ref), "fasta")
	except StandardError:
		print "Could not open file", options.ref
		sys.exit()
	contigseqs={}
	for contig in contigs:
		contigseqs[contig.id]=str(contig.seq)
	
	
	
	print "Calculating coverage for accessory file:", options.corebam
	sys.stdout.flush()
	
	try:
		accsamfile = pysam.Samfile( options.accbam, "rb" )
	except StandardError:
		print 'Failed to open '+options.accbam+'. Is it in bam format?'
		sys.exit()
	
	refs=accsamfile.references
	lengths=accsamfile.lengths
	
	contiglocs={}
	totallength=0
	
	for x, ref in enumerate(refs):
		contiglocs[ref]=totallength
		totallength+=lengths[x]
	
			
	
	
	poscount=0
	
	totcov=0.0
	totlen=0
	numkept=0
	
	keptcovs=[]
	filteredcovs=[]
	
	keepout=open(options.output+"_retained.fasta", "w")
	lenout=open(options.output+"_tooshort.fasta", "w")
	covout=open(options.output+"_lowcov.fasta", "w")
	covfile=open(options.output+"_coverage.csv", "w")
	print >> covfile, '\t'.join(["Core", str(coremean), "Core"])
	for x, ref in enumerate(refs):
		lastcolumn=-1
		zerocount=0
		contigcov=0.0
		contiglen=lengths[x]-(options.cliplen*2)
		if lengths[x]<options.minlen:
			print "Removed contig", ref, "on length =", lengths[x]
			if not ref in contigseqs:
				print "Error, contig not in contigseqs"
				sys.exit()
			print >> lenout, ">"+ref
			print >> lenout, contigseqs[ref]
			continue
		poscount=0
		for pileupcolumn in accsamfile.pileup(ref):
			
			poscount+=1
			while pileupcolumn.pos!=lastcolumn+1:
				lastcolumn+=1
				zerocount+=1
			
			if pileupcolumn.n>0 and poscount>options.cliplen and poscount<(lengths[x]-options.cliplen):
				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
					#print str(pileupcolumn)
					filtered_depth=0
					for pileupread in pileupcolumn.pileups:
						#print pileupread.alignment
						q=ord(pileupread.alignment.qual[pileupread.qpos])-33
						Q=pileupread.alignment.mapq
						#print q, Q, options.base_qual_filter, options.mapping_qual_filter
						if q>=options.base_qual_filter and Q>=options.mapping_qual_filter:
							filtered_depth+=1
					contigcov+=filtered_depth
					#sys.exit()
				else:
					contigcov+=pileupcolumn.n
			lastcolumn=pileupcolumn.pos
		if ((contigcov/contiglen)/coremean)*100<options.percent:
			print "Removed contig", ref, "on coverage =", contigcov/contiglen, "=", ((contigcov/contiglen)/coremean)*100, "percent of core"
			print >> covout, ">"+ref
			print >> covout, contigseqs[ref]
			filteredcovs.append(contigcov/contiglen)
			print >> covfile, '\t'.join([ref, str(contigcov/contiglen), 'removed'])
		else:
			print "Kept contig", ref, "length =", contiglen, "coverage =", contigcov/contiglen
			print >> keepout, ">"+ref
			print >> keepout, contigseqs[ref]
			print >> covfile, '\t'.join([ref, str(contigcov/contiglen), 'kept'])
			totcov+=contigcov
			totlen+=contiglen
			if (contigcov/contiglen)>2*coremean:
				keptcovs.append(2*coremean)
			else:
				keptcovs.append(contigcov/contiglen)
			numkept+=1
	
	keepout.close()
	lenout.close()
	covout.close()
	covfile.close()
	
	accmean=totcov/totlen
	print "Core mean =", coremean, "accessory mean =", accmean, "in", numkept, "contigs"
	
	legend_colours=[(pylab.Rectangle((0,0), 1, 1, fc="b")), (pylab.Rectangle((0,0), 1, 1, fc="r"))]
	n, bins, patches = plt.hist([keptcovs, filteredcovs], bins=50, color=["b", "r"], histtype="barstacked")
	
	ln=plt.axvline(color="g", x=coremean, lw=2, ls="--")
	cutoffline=plt.axvline(color="c", x=coremean/2, lw=2, ls="--")
	plt.legend(legend_colours+[ln]+[cutoffline], ["Retained", "Filtered", "Core mean", "Filter cutoff"])
	plt.xlabel("Mean coverage")
	plt.ylabel("Contig frequency")
	plt.savefig(options.output+"_coverage_histogram.pdf", format="pdf")
