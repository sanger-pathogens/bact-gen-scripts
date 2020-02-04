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
	
	options.separate_RGs=True
	RGs=[]
	if options.separate_RGs and 'RG' in samfile.header:
		for RG in samfile.header['RG']:
			RGs.append(RG['ID'])
	else:
		options.separate_RGs=False

	
	
	contiglocs={}
	totallength=0
	
	for x, ref in enumerate(refs):
		contiglocs[ref]=totallength
		totallength+=lengths[x]
	
			
	
	
	poscount=0
	if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
		print >> output, ' '.join(["#BASE", "High_Quality_Coverage", "Low_Quality_Coverage"])
	else:
		if options.separate_RGs:
			joinlist=["#BASE"]
			for RG in RGs:
				joinlist.append(RG)
			print >> output, ' '.join(joinlist)
		else:
			print >> output, ' '.join(["#BASE", "Coverage"])
	
	totlen=0
	

	for x, ref in enumerate(refs):
		print ref
		totlen+=lengths[x]
		lastcolumn=-1
		zerocount=0
		for pileupcolumn in samfile.pileup(ref):
			
			poscount+=1
			while pileupcolumn.pos!=lastcolumn+1:
				poscount+=1
				lastcolumn+=1
				zerocount+=1
				
			if options.separate_RGs:
				
				if pileupcolumn.n>0:
					if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
						RG_depth={}
						for RG in RGs:
							RG_depth[RG]=0
						for pileupread in pileupcolumn.pileups:
							q=ord(pileupread.alignment.qual[pileupread.qpos])-33
							Q=pileupread.alignment.mapq
							if q>=options.base_qual_filter and Q>=options.mapping_qual_filter:
								print pileupread.tags
								sys.exit()
						print >> output, ' '.join([str(poscount), str(filtered_depth), str(pileupcolumn.n-filtered_depth)])
						
					else:
						RG_depth={}
						for RG in RGs:
							RG_depth[RG]=0
						for pileupread in pileupcolumn.pileups:
							for tag in pileupread.alignment.tags:
								if tag[0]=='RG' and tag[1] in RG_depth:
									RG_depth[tag[1]]+=1
									break
						joinlist=[str(poscount)]
						for RG in RGs:
							joinlist.append(str(RG_depth[RG]))
						print >> output, ' '.join(joinlist)
			else:
				if pileupcolumn.n>0:
					if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
						filtered_depth=0
						for pileupread in pileupcolumn.pileups:
							q=ord(pileupread.alignment.qual[pileupread.qpos])-33
							Q=pileupread.alignment.mapq
							if q>=options.base_qual_filter and Q>=options.mapping_qual_filter:
								filtered_depth+=1
						print >> output, ' '.join([str(poscount), str(filtered_depth), str(pileupcolumn.n-filtered_depth)])
						
					else:
						print >> output, ' '.join([str(poscount), str(pileupcolumn.n)])
			lastcolumn=pileupcolumn.pos
		

		while lastcolumn+1<lengths[x]:
			poscount+=1
			lastcolumn+=1
			zerocount+=1

		if lastcolumn+1==lengths[x]:
			if options.separate_RGs:
				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
					print >> output, ' '.join([str(totlen), "0", "0"])
				else:
					joinlist=[str(totlen)]
					for RG in RGs:
						joinlist.append("0")
					print >> output, ' '.join(joinlist)
			else:
				if options.base_qual_filter!=0 or options.mapping_qual_filter!=0:
					print >> output, ' '.join([str(totlen), "0", "0"])
				else:
					print >> output, ' '.join([str(totlen), "0"])

	
	output.close()
