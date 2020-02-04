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
	
	
	parser.add_option("-o", "--output", action="store", dest="outputfile", help="output plot file name [default= %default]", type="string", metavar="FILE", default="")
	parser.add_option("-b", "--bam", action="store", dest="bamfile", help="Input bam file", type="string", metavar="FILE", default="")
	
	parser.add_option("-q", "--base_qual_filter", action="store", dest="base_qual_filter", help="Base quality filter for bam file mapping plots [default= %default]", default=0, type="float")
	parser.add_option("-Q", "--mapping_qual_filter", action="store", dest="mapping_qual_filter", help="Mapping quality filter for bam file plots [default= %default]", default=0, type="float")
	parser.add_option("-l", "--read_length", action="store_true", dest="read_length", help="Also include mean read length in plot [default= %default]", default=False)
	
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
	
	if options.outputfile=="":
		print 'No output file name given'
		sys.exit()
	
	if options.base_qual_filter<0 or options.base_qual_filter>100:
		print 'Base quality filter must be between 0 and 100'
		sys.exit()
	if options.mapping_qual_filter<0 or options.mapping_qual_filter>100:
		print 'Mapping quality filter must be between 0 and 100'
		sys.exit()
	
	output=open(options.outputfile,"w")

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
	
	for x, ref in enumerate(refs):
		contiglocs[ref]=totallength
		totallength+=lengths[x]
	
			
	
	
	poscount=0
	if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
		header=["#BASE", "High_Quality_Coverage", "Low_Quality_Coverage"]
		if options.read_length:
			header.append("Mean_read_length")
		print >> output, ' '.join(header)
	else:
		header=["#BASE", "Coverage"]
		if options.read_length:
			header.append("Mean_read_length")
		print >> output, ' '.join(header)
	
	totlen=0
	
	for x, ref in enumerate(refs):
		print ref
		totlen+=lengths[x]
		lastcolumn=-1
		zerocount=0
		for pileupcolumn in samfile.pileup(ref):
			poscount+=1
			#print pileupcolumn.pos, lastcolumn+1, poscount
			while pileupcolumn.pos!=lastcolumn+1:
				#print lastcolumn, pileupcolumn.pos
				poscount+=1
#				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
#					print >> output, ' '.join([str(poscount), "0", "0"])
#				else:
#					print >> output, ' '.join([str(poscount), "0"])
				lastcolumn+=1
				zerocount+=1
			
			if pileupcolumn.n>0:
				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
					#print str(pileupcolumn)
					filtered_depth=0
					readlengths=0.0
					for pileupread in pileupcolumn.pileups:
						readlengths+=pileupread.alignment.alen
						q=ord(pileupread.alignment.qual[pileupread.qpos])-33
						Q=pileupread.alignment.mapq
						#print q, Q, options.base_qual_filter, options.mapping_qual_filter
						if q>=options.base_qual_filter and Q>=options.mapping_qual_filter:
							filtered_depth+=1
					outlist=[str(poscount), str(filtered_depth), str(pileupcolumn.n-filtered_depth)]
					if options.read_length:
						outlist.append(str(readlengths/pileupcolumn.n))
					print >> output, ' '.join(outlist)
					#sys.exit()
				else:
					outlist=[str(poscount), str(pileupcolumn.n)]
					if options.read_length:
						readlengths=0.0
						for pileupread in pileupcolumn.pileups:
							readlengths+=pileupread.alignment.alen
						outlist.append(str(readlengths/pileupcolumn.n))
					print >> output, ' '.join(outlist)
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
#		print poscount, totlen, lengths[x], lastcolumn
		while lastcolumn+1<lengths[x]:
			poscount+=1
#			if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
#				print >> output, ' '.join([str(poscount), "0", "0"])
#			else:
#				print >> output, ' '.join([str(poscount), "0"])
			lastcolumn+=1
			zerocount+=1
#		print poscount
		if lastcolumn+1==lengths[x]:
			if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
				outlist=[str(totlen), "0", "0"]
			else:
				outlist=[str(totlen), "0"]
			if options.read_length:
				outlist.append("0")
			print >> output, ' '.join(outlist)
		
#		if ref=="5463_3#11_shuffled_141_cov_6":
#			sys.exit()
	
	output.close()
