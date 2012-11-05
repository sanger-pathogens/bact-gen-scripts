#!/usr/bin/env python
import os, sys, string
from random import *
from optparse import OptionParser, OptionGroup
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *




##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	group = OptionGroup(parser, "Input Options")

	group.add_option("-f", "--forward", action="store", dest="forward", help="forward fastq file (for single end, just use this)", default=False)
	group.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fastq file", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Output Options")
	group.add_option("-o", "--output", action="store", dest="outputname", help="output file name (will be appended withg _1.fastq and _2.fastq for forward and reverse reads, respectively)", default="")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Read Filtering Options")
	group.add_option("-d", "--database", action="store", dest="database", help="database of contaminants to be removed (optional)", default="")
	group.add_option("-q", "--quality", action="store", type="int", dest="quality_cutoff", help="Quality cutoff to trim reads to [default=%default]", default=15)
	group.add_option("-l", "--length", action="store", type="int", dest="length_cutoff", help="Length cutoff: remove reads if at least one of the pair is less than this length [default=%default]", default=36)
	group.add_option("-g", "--gc", action="store", dest="gcfilter", help="Filter reads by GC (format must be a letter a (above) or b (below) followed by a percentage) [default=%default]", default=False)
	parser.add_option_group(group)

	
	return parser.parse_args()



def check_input_options(options, args):

	single=False
	if not options.forward or not os.path.isfile(options.forward):
		print "Cannot find file", options.forward
		sys.exit()
	if not options.reverse or not os.path.isfile(options.reverse):
		print "Cannot find file",options.reverse
		single=True

	if options.database and not os.path.isfile(options.database):
		print "Cannot find file", options.database
		sys.exit()

	if options.gcfilter and not options.gcfilter[0] in ["a", "b"]:
		DoError("Invalid read GC cutoff. Must start with a (above) or b (below) followed by a float")
	if options.gcfilter:
		try:
			value=float(options.gcfilter[1:])
		except StandarError:
			DoError("Invalid read GC cutoff. Must start with a (above) or b (below) followed by a float")
		if value <0 or value >100:
			DoError("Invalid read GC cutoff. GC cutoff must be between 0 and 100")
		if options.gcfilter[1] not in ["a", "b"]:
			DoError("Invalid read GC cutoff. Must start with a (above) or b (below) followed by a float")


	return single


    
#######################################
# Function remove reads mapping to db #
#######################################


def remove_unmapping_reads(samfilename):
	
	
	MAX_INSERT=1000
	MIN_INSERT=0
	
	#CIGAR operators: 0=match, 1=insertion, 2=deletion, 3=skipped region from reference, 4=soft clip on read (i.e. sequence still in SEQ, 5=hard clip on read (sequence not in SEQ), 6=Padding (silent deletion from padded reference - for multiple alignment)
	def calculate_cigar_length(cigar_sequence):
		length=0
		
		for f in cigar_sequence:
			if f[0] in [0,1,3,4,5]:
				length+=f[1]
		
		return length
		
		
	def bamline2fastq(bamline, direction, handle=""):
		if handle=="":
			print "@"+bamline.qname+"/"+direction
			print bamline.seq
			print "+"
			print bamline.qual
		else:
			print >> handle, "@"+bamline.qname+"/"+direction
			print >> handle, bamline.seq
			print >> handle, "+"
			print >> handle, bamline.qual
			
		
	
	
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	
	
	if samfilename.split(".")[-1]=="bam":
	
		samfile = pysam.Samfile( samfilename, "rb" )
	
	
	elif samfilename.split(".")[-1]=="sam":
		samfile = pysam.Samfile( samfilename, "r" )
	else:
		print "Not a sam or bam file"
		sys.exit()	
	
	
	outputf=open(samfilename.split(".")[0]+"_1_unmapped.fastq", "w")
	if not single:
		outputr=open(samfilename.split(".")[0]+"_2_unmapped.fastq", "w")
			
		
	reads_iter=samfile
	
	addedcount=0
	
	for count, read in enumerate(reads_iter):
		
		mate = reads_iter.next()
		
		if read.is_reverse:
			readtmp=mate
			mate=read
			read=readtmp
		
		#print read, mate
		
		#if  (not read.is_proper_pair or read.opt("XT")!=85) or (not mate.is_proper_pair or  mate.opt("XT")!=85) or read.isize>MAX_INSERT or read.isize<MIN_INSERT:
		
		if  read.is_unmapped and mate.is_unmapped:
	
			readqual=0
			toprint=True
			for base, basequal in enumerate(read.qual):
				if ord(basequal)-33<15 or ord(mate.qual[base])<15:
					toprint=False
					#print read, mate
					break	
			
			if toprint:
				
				bamline2fastq(read, "1", outputf)
				if not single:
					bamline2fastq(mate, "2", outputr)
					
				addedcount+=1
		
	
	outputf.close()
	outputr.close()






if __name__ == "__main__":


	if sys.version.startswith("3"):
		import io
		io_method = io.BytesIO
	else:
		import cStringIO
		io_method = cStringIO.StringIO


	(options, args) = main()
	
	single=check_input_options(options, args)
	
	
	forward=options.forward
	reverse=options.reverse
	outputname=options.outputname
	db=options.database

	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	currforwardfilename=forward
	currreversefilename=reverse
	
	todelete=""
	
	
	if currforwardfilename.split(".")[-1]=="gz":
		print "Unzipping", currforwardfilename, "to", tmpname+"_1.fastq"
		os.system("zcat "+currforwardfilename+" > "+tmpname+"_1.fastq")
		currforwardfilename=tmpname+"_1.fastq"
		todelete=todelete+currforwardfilename+" "
	if not single:
		if currreversefilename.split(".")[-1]=="gz":
			print "Unzipping", currreversefilename, "to", tmpname+"_2.fastq"
			os.system("zcat "+currreversefilename+" > "+tmpname+"_2.fastq")
			currreversefilename=tmpname+"_2.fastq"
			todelete=todelete+currreversefilename+" "
	
	
	if not single and  db!="" and os.path.isfile(db):
		print "Mapping reads vs", db
		os.system("cp "+db+" "+tmpname+"_db.fasta")
		os.system("bwa index "+tmpname+"_db.fasta")
		os.system(BWA_DIR+"bwa aln -q 15 "+tmpname+"_db.fasta "+currforwardfilename+" > "+tmpname+".F.sai")
		os.system(BWA_DIR+"bwa aln -q 15 "+tmpname+"_db.fasta "+currreversefilename+" > "+tmpname+".R.sai")
		#Join both aligments
		os.system(BWA_DIR+"bwa sampe "+tmpname+"_db.fasta "+tmpname+".F.sai "+tmpname+".R.sai "+currforwardfilename+" "+currreversefilename+" > "+tmpname+".sam")
		print "Finding reads that do not map against", db
		remove_unmapping_reads(tmpname+".sam")
		
		os.system("rm "+tmpname+".sam")
		currforwardfilename=tmpname+"_1_unmapped.fastq"
		currreversefilename=tmpname+"_2_unmapped.fastq"
		if todelete!="":
			os.system("rm -f "+todelete)
		todelete=currforwardfilename+" "+currreversefilename
		

	
	
	#Filter reads for length and GC
	
	
	filter_fastq(currforwardfilename, reversefile=currreversefilename, shuffled=False, quality_cutoff=options.quality_cutoff, length_cutoff=options.length_cutoff, GCcutoff=options.gcfilter)
	
	currforwardfilename='.'.join(currforwardfilename.split(".")[:-1])+"_filtered."+currforwardfilename.split(".")[-1]
	if not single:
		currreversefilename='.'.join(currreversefilename.split(".")[:-1])+"_filtered."+currreversefilename.split(".")[-1]
	
	if todelete!="":
		os.system("rm -f "+todelete)
	
	if single:
		os.system("mv "+currforwardfilename+" "+options.outputname+".fastq")
	else:
		os.system("mv "+currforwardfilename+" "+options.outputname+"_1.fastq")
		os.system("mv "+currreversefilename+" "+options.outputname+"_2.fastq")

	print "Done"
