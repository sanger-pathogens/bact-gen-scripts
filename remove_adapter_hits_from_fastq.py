#!/usr/bin/env python

import os, sys
from optparse import OptionParser
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-f", "--fastq", action="store", dest="fastq", help="fastq file name", default="", metavar="FILE")
	parser.add_option("-d", "--database", action="store", dest="adapterdb", help="location of fasta file continaing adapter (or other)  DNA sequeces to remove from fastq [Default= %default]", default="/nfs/users/nfs_s/sh16/data/adapters.fasta", metavar="FILE")
	parser.add_option("-e", "--evalue", action="store", dest="evalue", help="evalue cutoff for blast [Default= %default]", type="float", default=0.0001, metavar="FLOAT")
	parser.add_option("-p", "--paired", action="store_true", dest="paired", help="use if the fastq file has a matching mate-pair. The two files must have the suffix (_1.fastq and _2.fastq). Either can be given as the -f argument. If the input file is a shuffled fastq containing both pairs, do not use this option.", default=False, metavar="BOOLEAN")
	parser.add_option("-a", "--processors", action="store_true", dest="processors", help="Number of processors to use for blast [Default= %default]", default=1, metavar="int")
	
	return parser.parse_args()
	
	
################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.fastq=='':
		DoError('No fastq file selected!')
	elif not os.path.isfile(options.fastq):
		DoError('Cannot find file '+options.fastq+'!')
	if options.paired:
	
		if options.fastq.replace(".fastq","")[-2]=="_" and options.fastq.replace(".fastq","")[-1]=="1":
			pair=options.fastq.replace("_1.fastq","_2.fastq")
		elif options.fastq.replace(".fastq","")[-2]=="_" and options.fastq.replace(".fastq","")[-1]=="2":
			pair=options.fastq.replace("_2.fastq","_1.fastq")
		else:
			DoError('Fastq file '+options.fastq+' has invalid suffix for a paired file (should be _1.fastq or _2.fastq)!')
		if not os.path.isfile(pair):
			DoError('Cannot find mate file ('+pair+') corresponding to fastq file '+options.fastq+'!')
	
	if options.adapterdb=='':
		DoError('No fastq file selected!')
	elif not os.path.isfile(options.adapterdb):
		DoError('Cannot find file '+options.adapterdb+'!')
	elif options.evalue>1 or options.evalue<0:
		DoError('Evalue (-e) must be between 0 and 1!')
		
	return
	



################
# Main program #
################		

if __name__ == "__main__":
	
	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	fastqs=[options.fastq]
	evalue=options.evalue
	adapterdb=options.adapterdb
	paired=options.paired
	if options.paired:
	
		if options.fastq.replace(".fastq","")[-2]=="_" and options.fastq.replace(".fastq","")[-1]=="1":
			pair=options.fastq.replace("_1.fastq","_2.fastq")
		elif options.fastq.replace(".fastq","")[-2]=="_" and options.fastq.replace(".fastq","")[-1]=="2":
			pair=options.fastq.replace("_2.fastq","_1.fastq")
	
		fastqs.append(pair)
	
	#check if sequences to be searched for have been formatted into a db
	print "Formatting database"
	sys.stdout.flush()
	if not os.path.isfile(adapterdb+".nhr") or not os.path.isfile(adapterdb+".nin") or not os.path.isfile(adapterdb+".nsq"):
		os.system("formatdb -p F -i "+adapterdb)
		os.system("rm formatdb.log")
		
	
	blastlist=set()
	
	for fastq in fastqs:
	
		#Create fasta from fastq
		print "Converting fastq to fasta"
		sys.stdout.flush()
		os.system("/nfs/users/nfs_s/sh16/scripts/fastq2fasta.pl "+fastq+" "+fastq+".fasta")
	
		#Run blast
		print "Running blast on "+fastq+" (this may take a while)"
		sys.stdout.flush()
		os.system("blastall -p blastn -v 1 -b 1 -m 8 -o "+fastq+".blast -e "+str(evalue)+" -a "+str(options.processors)+" -i "+fastq+".fasta -d "+adapterdb+" > /dev/null 2>&1")
	
		os.system("rm "+fastq+".fasta")
		blast=open(fastq+".blast")
		
		print "Parsing blast hits"
		sys.stdout.flush()
		
		for line in blast:
			if line.split()[0][:-2] not in blastlist:#len(blastlist)==0 or 
				blastlist.add(line.split()[0][:-2])
		blast.close()
		os.system("rm "+fastq+".blast")
	
	for fastq in fastqs:
		print "Filtering blast hits from "+fastq+" into "+fastq.replace(".fastq", "_filtered.fastq")
		sys.stdout.flush()
		fastqfile=open(fastq)
		output=open(fastq.replace(".fastq", "_filtered.fastq"),"w")
		
		count=0
		for line in fastqfile:
			if line.strip().replace("@","")[:-2] not in blastlist:
				print >> output, line.strip()
				for x in range(0,3):
					print >> output, fastqfile.next().strip()
				count+=1
			else:
				for x in range(0,3):
					output, fastqfile.next().strip()
		fastqfile.close()
		output.close()
		print count, "reads remaining"
		sys.stdout.flush()
