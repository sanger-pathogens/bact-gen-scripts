#!/usr/bin/env python

import os, sys
from optparse import OptionParser, OptionGroup
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
import subprocess
import shlex
import gzip
import shutil
from random import *
import math
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
import farm

##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-q", "--query", action="store", dest="query", help="Fasta file to use as the BLAST database", default="", metavar="FILE")
	group.add_option("-s", "--subject", action="store", dest="subject", help="Fasta file to use as the BLAST query", default="", metavar="FILE")
	group.add_option("-o", "--output", action="store", dest="prefix", help="Prefix for output file(s). Note: Currently the script only produces a gzipped crunch file. Please ask for additional files to be added. [Default= %default]", default="", metavar="FILE")
	group.add_option("-Q", "--queue", action="store", dest="queue", help="Queue to bsub to (choose from normal, long or basement). [Default= %default]", default="normal", type="choice", choices=['normal', 'long', 'basement'])
	group.add_option("-p", "--program", action="store", dest="blastprog", help="BLAST program to use (choose from blastn or tblastx) This script doesn't support protein BLASTs, and won't unless people ask very nicely. [Default= %default]", default="blastn", type="choice", choices=['blastn', 'tblastx'])
	group.add_option("-e", "--evalue", action="store", dest="e", help="evalue cutoff for BLAST. [Default= %default]", default=0.00001, type="float")
	group.add_option("-E", "--extras", action="store", dest="extras", help="Extra BLAST options to use. Note: These will not be sanity checked.", default="")
	group.add_option("-d", "--tmpdir", action="store", dest="tmpdir", help="Temporary directory name. [Default= %default]", default="better_blast_tmp_dir")
	
	parser.add_option_group(group)


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.query=='':
		DoError('No query file selected')
	elif not os.path.isfile(options.query):
		DoError('Cannot find file '+options.query)
	if options.subject=='':
		DoError('No subject file selected')
	elif not os.path.isfile(options.subject):
		DoError('Cannot find file '+options.subject)
	if options.tmpdir=='':
		options.tmpdir="better_blast_tmp_dir"
	
	
	return



################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	if options.prefix=="":
		options.prefix=os.path.basename(options.query)+"_vs_"+os.path.basename(options.subject)
	
	#Create a temporary name
	
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	
	if not os.path.exists(options.tmpdir):
		try:
			os.makedirs(options.tmpdir)
		except StandardError:
			DoError("failed to make better blast directory")
		rem_dir=True
	else:
		rem_dir=False
	
	
	
	query_info_file_name=options.tmpdir+"/"+tmpname+".queryinfo"
	query_info_file=open(query_info_file_name, "w")
	
	print "Preprocessing query fasta file"
	sys.stdout.flush()
	
	
	try:
		query_seqs=read_seq_file(options.query)
	except StandardError:
		DoError("Cannot read query file")
	
	query_tot=0
	query_tmp_file_name=options.tmpdir+"/"+tmpname+"_query.fasta"
	query_fragment_count=0
	query_file_count=1
	curr_file_len=0
	query_tmp_file=open(query_tmp_file_name+str(query_file_count), "w")
	for query in query_seqs:
		seqlength=len(str(query.seq))
		if seqlength>50000:
			num_fragments=int(math.ceil(float(seqlength)/50000))
			fragment_length=int(math.ceil(float(seqlength)/num_fragments))
			x=0
			for x in xrange(num_fragments):
				i=x*fragment_length
				if i>seqlength:
					break
				j=(x+1)*fragment_length
				if j>seqlength:	
					j=seqlength
				query_tmp_file.close()
				query_file_count+=1
				curr_file_len=0
				query_tmp_file=open(query_tmp_file_name+str(query_file_count), "w")
				print >> query_info_file, '\t'.join(["q"+str(query_fragment_count), query.id, str(query_tot)])
				print >> query_tmp_file, ">q"+str(query_fragment_count)
				print >> query_tmp_file, str(query.seq)[i:j]
				query_fragment_count+=1
				curr_file_len+=j-i
				query_tot+=j-i
		else:
			if (curr_file_len+seqlength)>50000:
				query_tmp_file.close()
				query_file_count+=1
				curr_file_len=0
				query_tmp_file=open(query_tmp_file_name+str(query_file_count), "w")
			print >> query_info_file, '\t'.join(["q"+str(query_fragment_count), query.id, str(query_tot)])
			print >> query_tmp_file, ">q"+str(query_fragment_count)
			print >> query_tmp_file, str(query.seq)
			query_fragment_count+=1
			curr_file_len+=seqlength
			query_tot+=seqlength
			
	query_tmp_file.close()
	query_info_file.close()
	
	
	print "Preprocessing subject fasta file"
	sys.stdout.flush()
	
	subject_info_file_name=options.tmpdir+"/"+tmpname+".subjectinfo"
	subject_info_file=open(subject_info_file_name, "w")
	
	try:
		subject_seqs=read_seq_file(options.subject)
	except StandardError:
		DoError("Cannot read subject file")
	
	subject_tot=0
	subject_tmp_file_name=options.tmpdir+"/"+tmpname+"_subject.fasta"
	subject_tmp_file=open(subject_tmp_file_name, "w")
	subject_fragment_count=0
	for subject in subject_seqs:
		seqlength=len(str(subject.seq))
		print >> subject_info_file, '\t'.join(["s"+str(subject_fragment_count), subject.id, str(subject_tot)])
		print >> subject_tmp_file, ">s"+str(subject_fragment_count)
		print >> subject_tmp_file, str(subject.seq)
		subject_tot+=seqlength
		subject_fragment_count+=1
	subject_tmp_file.close()
	
	
	print "formatting subject database"
	sys.stdout.flush()
	
	formatdb_command="formatdb -p F -i "+subject_tmp_file_name
	formatdb_args = shlex.split(formatdb_command)
	formatdb_returnval=subprocess.call(formatdb_args)

	if formatdb_returnval!=0:
		DoError("formatdb command "+formatdb_command+" failed with exit value "+str(formatdb_returnval))



	print "Running blast jobs over lsf"
	sys.stdout.flush()
	
	blast_command="blastall -p "+options.blastprog+" -i "+query_tmp_file_name+"INDEX -m 8 -e "+str(options.e)+" -o "+options.tmpdir+"/"+tmpname+".blastINDEX -d "+subject_tmp_file_name+" "+options.extras
	
	job1 = farm.Bsub(options.tmpdir+"/"+tmpname+"_bb_bsub.out", options.tmpdir+"/"+tmpname+"_bb_bsub.err", tmpname+"_blast", "normal", 0.5, blast_command, start=1, end=query_file_count)
	job1_id = job1.run()
	
	print "Job ID =", job1_id
	sys.stdout.flush()
	
	
	
	job2 = farm.Bsub(options.prefix+"_bb_bsub.out", options.prefix+"_bb_bsub.err", tmpname+"_postprocess", "normal", 0.5, "~/scripts/post_process_better_blast.py "+options.tmpdir+" "+tmpname+" "+options.prefix+" "+str(query_file_count))
	job2.add_dependency(job1_id) 
	job2_id = job2.run()
	print "Job ID =", job2_id
	sys.stdout.flush()
	
	print "Once these jobs are finished, you can run a comparison in act using:"
	print "act", options.query, options.prefix+".crunch.gz", options.subject
	
	job3 = farm.Bsub(options.prefix+"_bb_bsub.out", options.prefix+"_bb_bsub.err", tmpname+"_cleanup", "normal", 0.5, "rm -rf formatdb.log "+options.tmpdir)
	job3.add_dependency(job2_id) 
	job3_id = job3.run()
		
	
	
