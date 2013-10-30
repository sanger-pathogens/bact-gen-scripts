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
	
	for x, ref in enumerate(refs):
		contiglocs[ref]=totallength
		totallength+=lengths[x]
		indelreads[ref]={"insertion":0, "deletion":0, "clipped":0}
	
	
	for read in samfile:
		
		if not read.is_unmapped:
			
			for cig in read.cigar:
				if cig[0]==1:
					indelreads[refs[read.tid]]["insertion"]+=1
				elif cig[0]==2:
					indelreads[refs[read.tid]]["deletion"]+=1
				elif cig[0]==4:
					indelreads[refs[read.tid]]["clipped"]+=1
	
	for ref in indelreads:
		if indelreads[ref]["clipped"]>0 or indelreads[ref]["deletion"]>0 or indelreads[ref]["insertion"]>0:
			print ref,  indelreads[ref]
	sys.exit()
	
	
	
	
	
 	#This code searches for repeated mappings
#	low_qual_reads={}
#	
#	for read in samfile:
#	
#		if read.is_read1:
#			fr="1"
#		else:
#			fr="2"
#		name=read.qname+"/"+fr
#	
#		if not read.is_unmapped and read.mapq<options.mapping_qual_filter:
#			if not name in low_qual_reads:
#				low_qual_reads[name]=[0,0,set([]),[]]
#			low_qual_reads[name][0]+=1
#			low_qual_reads[name][2].add(refs[read.tid])
#		elif not read.is_unmapped and name in low_qual_reads:
#			low_qual_reads[name][1]+=1
#			low_qual_reads[name][3].append(refs[read.tid])
#		
#			
#		
#	shared_ref_low_qual_reads={}
#	
#	for name in low_qual_reads:
#		if len(low_qual_reads[name][2])>1:
#			lqr_list=list(low_qual_reads[name][2])
#			lqr_list.sort()
#			for x in xrange(len(lqr_list)):
#				if not lqr_list[x] in shared_ref_low_qual_reads:
#					shared_ref_low_qual_reads[lqr_list[x]]={}
#				for y in xrange(x+1, len(lqr_list)):
#					if not lqr_list[y] in shared_ref_low_qual_reads[lqr_list[x]]:
#						shared_ref_low_qual_reads[lqr_list[x]][lqr_list[y]]=set([])
#					shared_ref_low_qual_reads[lqr_list[x]][lqr_list[y]].add(name)
#		
#	
#	for ref1 in shared_ref_low_qual_reads:
#		for ref2 in shared_ref_low_qual_reads[ref1]:
#			if ref1!=ref2:
#				print ref1, ref2, shared_ref_low_qual_reads[ref1][ref2]
#		
#	
#	samfile.reset()
#	
#	
#	sys.exit()
	
	
	
	
	poscount=0
	
	totlen=0
	genes_to_keep=[]
	
	for x, ref in enumerate(refs):
	
		totlen+=lengths[x]
		lastcolumn=-1
		zerocount=0
		poscount=-1
		high_quality_coverage=[0]*lengths[x]
		low_quality_coverage=[0]*lengths[x]
		mapped=0
		for pileupcolumn in samfile.pileup(ref):
			
			poscount+=1
			
			while pileupcolumn.pos!=lastcolumn+1:
				
				poscount+=1
#				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
#					print >> output, ' '.join([str(poscount), "0", "0"])
#				else:
#					print >> output, ' '.join([str(poscount), "0"])
				lastcolumn+=1
				zerocount+=1
			
			if pileupcolumn.n>0:
				mapped=1
				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
					filtered_depth=0
					for pileupread in pileupcolumn.pileups:
						
						q=ord(pileupread.alignment.qual[pileupread.qpos])-33
						Q=pileupread.alignment.mapq
						#print q, Q, options.base_qual_filter, options.mapping_qual_filter
						if q>=options.base_qual_filter and Q>=options.mapping_qual_filter:
							filtered_depth+=1
					high_quality_coverage[poscount]=filtered_depth
					low_quality_coverage[poscount]=pileupcolumn.n-filtered_depth				
				else:
					high_quality_coverage[poscount]=pileupcolumn.n
					low_quality_coverage[poscount]=0.0
			lastcolumn=pileupcolumn.pos
		
		
#		while lastcolumn+1<lengths[x]:
#			poscount+=1
#			if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
#				print >> output, ' '.join([str(poscount), "0", "0"])
#			else:
#				print >> output, ' '.join([str(poscount), "0"])
#			lastcolumn+=1
#			zerocount+=1
		#print len(depths)
		while lastcolumn+1<lengths[x]:
			poscount+=1
#			if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
#				print >> output, ' '.join([str(poscount), "0", "0"])
#			else:
#				print >> output, ' '.join([str(poscount), "0"])
			lastcolumn+=1
			zerocount+=1
		
				
		high_percent=0.0
		low_percent=0.0
		both_percent=0.0
		either_percent=0.0
		if mapped==1:
			
			for y, base in enumerate(high_quality_coverage):
				if base>0 and low_quality_coverage[y]>0:
					high_percent+=1
					both_percent+=1
					low_percent+=1
					either_percent+=1
				elif base>0:
					high_percent+=1
					either_percent+=1
				elif low_quality_coverage[y]>0:
					low_percent+=1
					either_percent+=1
			
			high_percent=(high_percent/lengths[x])*100
			low_percent=(low_percent/lengths[x])*100
			both_percent=(both_percent/lengths[x])*100
			either_percent=(either_percent/lengths[x])*100
			
		if high_percent!=0:
			print refs[x], high_percent, low_percent, both_percent, either_percent
			genes_to_keep.append(refs[x])
	
	
	for gene in genes_to_keep:
		print gene
		returnval=os.system("~sh16/scripts/resistome/bam_filter.py -t contigs -c "+gene+" -f pairedfastq -o tmp_"+gene+" -b "+options.bamfile)
		print "bam_filter return value:", returnval
		
		returnval=os.system("velveth vdir 57 -shortPaired -fastq -separate  tmp_"+gene+"_1.fastq tmp_"+gene+"_2.fastq")
		returnval=os.system("velvetg vdir -ins_length 250 -min_contig_lgth 200 -exp_cov auto -cov_cutoff auto")
		returnval=os.system("cat vdir/contigs.fa >> tmp_db.fasta")
		
		
			
		
		
		
			
					
			
			#print high_quality_coverage
			#print low_quality_coverage
		
#		if ref=="5463_3#11_shuffled_141_cov_6":
#			sys.exit()
	
	#print genes_to_keep
	
