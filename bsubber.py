#!/usr/bin/env python


##################
# Import modules #
##################

import string, re
import os, sys
from optparse import OptionParser, OptionGroup
from socket import gethostname
import shlex
import subprocess

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

def get_user_options(args=[]):
	usage = "usage: %prog [options] <script/program to bsub>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	parser.add_option("-q", "--queue", action="store", dest="LSFQ", help="LSF queue to submit to. [Default= %default]", default="normal", type="choice", choices=["normal","long", "basement", "hugemem", "yesterday"])
	parser.add_option("-m", "--memory", action="store", dest="mem", help="Amount of memory required for analysis (Gb). [Default= None]", default=False, type="float")
	parser.add_option("-p", "--processors", action="store", dest="CPUs", help="Number of processors to use on each node. [Default= %default]", default=1, type="int")
	parser.add_option("-c", "--c", action="store_true", dest="checkpoint", help="Checkpoint your job (so it can be restarted if it dies due to hardware failure). [Default= False]", default=False)
	parser.add_option("-d", "--checkpoint_directory", action="store", dest="checkpoint_directory", help="Directory to store checkpoints in. [Default= pwd]", default="")
	parser.add_option("-f", "--checkpoint_frequency", action="store", dest="checkpoint_period", help="Frequency to save checkpoints in minutes. [Default= %default]", default=120, type="int")
	parser.add_option("-r", "--restart", action="store_true", dest="restart", help="Restart a checkpointed job. To do this you will need to provide the path to the checkpoint file including the jobid (you will normally probably want to use the last one in the directory), e.g. /my/checkpoint/directory/<JOBID>. Important: Check that your job has really died; restarting a job from a checkpoint whilst the original is running will lead to unpredictable results! If your jobs needs special resources, you should add them to the command. (LSF does not remember them from the original bsub). ", default=False)
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file name. [Default= None]", default=False)
	parser.add_option("-e", "--error", action="store", dest="error", help="Error file name. [Default= None]", default=False)
	parser.add_option("-E", "--exclude", action="store", dest="exclude", help="Comma separated list of nodes to exclude. [Default= None]", default="")
	
	if args==[]:
		return parser.parse_args()
	else:
		return parser.parse_args(args)


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.LSFQ!="hugemem" and (options.mem>30 or options.mem<0):
		DoError('Memory requirement (-M) must be between 0 and 30Gb')
	elif options.LSFQ=="hugemem" and (options.mem>250 or options.mem<30):
		DoError('Memory requirement (-M) for hugemem queue must be between 30 and 250Gb')
	if options.CPUs>10 or options.CPUs<1:
		DoError('Number of CPUs must be between 1 and 10')
	if options.checkpoint and not options.checkpoint_directory:
		options.checkpoint_directory=os.getcwd()
	elif options.checkpoint and not os.path.isdir(options.checkpoint_directory):
		DoError('Checkpoint directory '+options.checkpoint_directory+' does not exist')
	elif options.restart and not os.path.isdir(options.checkpoint_directory):
		DoError('Checkpoint directory '+options.checkpoint_directory+' does not exist')
#	if options.checkpoint and host[:4]=="pcs4":
#		DoError('Checkpointing does not work on pcs4. You will need to long on to the farm to use this option.')
	if options.checkpoint and (options.checkpoint_period<1):
		DoError('Checkpoint frequency must be greater than 1')
#	if options.restart and not os.path.isfile(args[0]):
#		DoError('Checkpoint file '+args[0]+' does not exist')
		
	return

####################
# Get cluster name #
####################

def getclustername():
	mycluster="unknown"
	try:
		lsid_output=subprocess.check_output(["lsid"])
		
		for line in lsid_output.split("\n"):
			words=line.strip().split()
			if len(words)>0:
				if words[1]=="cluster":
					mycluster=words[4]
	
		
	except StandardError:
		return mycluster
	
	return mycluster

		
########
# Main #
########


if __name__ == "__main__":

	host=getclustername()
	print "Running on", host
	(options, args)=get_user_options()
	
	if len(args)==1:
		if len(args[0])==0:
			print "No command found to bsub"
			sys.exit()
		args[0].strip()
		usage = "usage: %prog [options] <script/program to bsub>"
		version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
		parser = OptionParser(usage=usage, version=version)
		(options, args)=get_user_options(args=sys.argv[1:-1]+shlex.split(args[0]))
		
	#print options, args
	check_input_validity(options, args)
	
	
	options.checkpoint_directory=os.path.abspath(options.checkpoint_directory)
	
	rlist=[]
	rstring=""
	ostring=""
	estring=""
	cstring=""
	crstring=""
	memstring=""
	cpustring=""
	excludestring=""
	
	if options.exclude!="":
		toexclude=options.exclude.split(",")
		excludestring='-R "select['
		for x, exc in enumerate(toexclude):
			excludestring=excludestring+"hname!='"+exc+"'"
			if x<len(toexclude)-1:
				excludestring=excludestring+" && "
		excludestring=excludestring+']"'
	
	if options.mem>0 or options.CPUs>1:
		rlist.append("-R")
		rlist.append("'")
	
	if options.CPUs>1:
		rlist.append('span[hosts=1]')
		cpustring="-n "+str(options.CPUs)
	
	if options.mem>0:
		if host=="farm3" or host=="pcs5":
			memlimit=str(int(options.mem*1000))
		else:
			memlimit=str(int(options.mem*1000000))
		memresource=str(int(options.mem*1000))
		if not options.restart:
			rlist.append('select[mem>'+memresource+'] rusage[mem='+memresource+']')
		memstring="-M "+memlimit
	
	if options.checkpoint:
		cstring=' '.join(["-k '"+options.checkpoint_directory, "method=blcr", str(options.checkpoint_period)+"'"])
		crstring="cr_run"
	
	if len(rlist)>2:
		rlist.append("'")
		rstring=' '.join(rlist)
	
	
	if options.LSFQ:
		qstring='-q '+options.LSFQ
	if options.output:
		ostring="-o "+options.output
	if options.error:
		estring="-e "+options.error
	
	for x in xrange(len(args)):
		args[x]=args[x].replace("(","\(")
		args[x]=args[x].replace(")","\)")
		args[x]=args[x].replace(" ","\ ")
		args[x]=args[x].replace('\\"','"')
	
	if options.restart:
		bsubstring="brestart"
		if options.checkpoint_directory[-1]=="/":
			options.checkpoint_directory=options.checkpoint_directory[:-1]
		arguments='/'.join(options.checkpoint_directory.split("/")[:-1])+" "+options.checkpoint_directory.split("/")[-1]
	else:
		bsubstring="bsub"
		if options.checkpoint:
			arguments="'"+' '.join([crstring]+args)+"'"
		else:
			arguments="'"+' '.join(args)+"'"
	
	
	#print arguments
	
	if len(arguments)==0:
		print "Nothing to bsub"
		sys.exit()		
	
	submitstring= ' '.join(' '.join([bsubstring, cstring, memstring, cpustring, rstring, qstring, ostring, estring, excludestring, arguments]).split())
	print submitstring
	#sys.exit()
	os.system(submitstring)
	
	
