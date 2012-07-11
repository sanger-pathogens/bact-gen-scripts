#!/usr/bin/env python


##################
# Import modules #
##################

import string, re
import os, sys
from optparse import OptionParser, OptionGroup



##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options] <script/program to bsub>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2011"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	parser.add_option("-q", "--queue", action="store", dest="LSFQ", help="LSF queue to submit to. [Default= %default]", default="normal", type="choice", choices=["normal","long", "basement", "hugemem"])
	parser.add_option("-m", "--memory", action="store", dest="mem", help="Amount of memory required for analysis (Gb). [Default= None]", default=False, type="int")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file name. [Default= None]", default=False)
	parser.add_option("-e", "--error", action="store", dest="error", help="Error file name (Gb). [Default= None]", default=False)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.LSFQ!="hugemem" and (options.mem>30 or options.mem<0):
		DoError('Memory requirement (-M) must be between 0 and 30Gb')
	elif options.LSFQ=="hugemem" and (options.mem>250 or options.mem<30):
		DoError('Memory requirement (-M) for hugemem queue must be between 30 and 250Gb')

	
	return


		
########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)
	
	rstring=""
	ostring=""
	estring=""
	memstring=""
	memrstring=""
	
	if options.mem>0:
		memkb=str(options.mem*1000000)
		memmb=str(options.mem*1000)
		memrstring='\'select[mem>'+memmb+'] rusage[mem='+memmb+']\''
		memstring="-M "+memkb
	
	if memstring:
		rstring=' '.join(["-R", memrstring])
	
	if options.LSFQ:
		qstring='-q '+options.LSFQ
	if options.output:
		ostring="-o "+options.output
	if options.error:
		estring="-e "+options.error
	
	
	bsubstring= ' '.join(["bsub", memstring, rstring, qstring, ostring, estring, '"', ' '.join(args), '"'])
	print bsubstring
	os.system(bsubstring)
	
	