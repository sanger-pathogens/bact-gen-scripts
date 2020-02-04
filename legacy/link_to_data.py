#!/usr/bin/env python


##################
# Import modules #
##################

import string, re, os, sys
from optparse import OptionParser

##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options] <list of fastq/bam files to be mapped>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2011"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	parser.add_option("-s", "--studyid", action="store", dest="study", help="study id", default=False, metavar="INT", type="int")
	parser.add_option("-t", "--type", action="store", dest="type", choices=["bam","fastq"], help="type of file to return (choice of fastq or bam) [default= %default]", metavar="STRING", type="choice", default="fastq")
	parser.add_option("-o", "--output", action="store", dest="output", help="output file name for isolate data list (leave blank for no file)", default="", metavar="FILE", type="str")

	return parser.parse_args()

################################
# Check command line arguments #
################################

def check_input_validity(options, args):
	
	force=False
	
	if options.output!='':
		while force==False and os.path.isfile(options.output):
			outopt=""
			outopt=raw_input('\nAn output file with the chosen name already exists.\n\nWould you like to overwrite (o), choose a new output file name (n) or quit (Q): ')
			if outopt=='Q':
				sys.exit()
			elif outopt=="o":
				force=True
			elif outopt=="n":
				options.output=raw_input('Enter a new output file name: ')
			
	
	if options.study=='':
		DoError('No study (-s) selected!')
	
	return
	
	
########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)
	
	if options.output!="":
		output=open(options.output, "w")
		print >> output, "\t".join(["Name", "Run", "Lane", "Tag"])
	
	lanes=set()
	
	#run plexfind
	
	print "Finding lanes for study id", options.study
	sys.stdout.flush()
	
	for line in os.popen("plexfind -s "+str(options.study)):
	
		words=line.strip().replace(" ","").split(",")
		if not '_'.join(words[1:2]) in lanes:
			lanes.add('_'.join(words[1:2]))
		
		if options.output!="":
			print >> output, "\t".join([words[1]+"_"+words[2], words[0], words[1].split("_")[0], words[1].split("_")[1], words[2]])
	
	if options.output!="":
		output.close()
	
	print "Found", len(lanes), "lanes"
	sys.stdout.flush()
	
	if len(lanes)==0:
		sys.exit()
	
	print "Creating links to", options.type, "files..."
	sys.stdout.flush()
	
	for lane in lanes:
		print lane
		sys.stdout.flush()
		for location in os.popen("pathfind -l "+lane+" -t "+options.type):
			sys.stdout.flush()
			location=location.strip()
			filename=location.split("/")[-1].replace("#", "_")
			
			os.system("cp -s "+location+" "+filename)
			
	
	
