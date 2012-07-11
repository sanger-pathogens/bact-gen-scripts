#!/usr/bin/env python
import os, sys, string
from random import *
from optparse import OptionParser, OptionGroup
import pysam
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import subprocess
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *


########################
# Define some globals #
#######################

SAMTOOLS_DIR=""
BWA_DIR=""
SOAPdir="/software/pathogen/external/apps/usr/bin/"
velvet_dir="/nfs/users/nfs_s/sh16/velvet-sc/"

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	group = OptionGroup(parser, "Input Options")

	group.add_option("-f", "--forward", action="store", dest="forward", help="forward fastq file", default=False)
	group.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fastq file", default=False)
	group.add_option("-s", "--shuffled", action="store", dest="shuffled", help="shuffled fastq file (cannot be used with contaminant removal step or SOAPdenovo)", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Output Options")
	group.add_option("-o", "--output", action="store", dest="outputname", help="output file name [default=%default]", default="merged.fasta")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Assembly Options")
	group.add_option("-n", "--number", action="store", type="int", dest="readspersplit", help="number of reads to put in each partition (0=use one partition) [default=%default]", default=0)
	group.add_option("-m", "--maxreads", action="store", type="int", dest="maxreads", help="number of reads from the fastq file to use for assembly (these will be randomly selected). This will override the -n option", default=0)
#	group.add_option("-S", "--scaffold", action="store_true", dest="scaffold", help="Create scaffolds rather than contigs [default=%default]", default=False)
	parser.add_option("-l", "--readlength", action="store", dest="readlength", help="Read length (required for multi-line fasta files)", default=0, type="int", metavar="length")
	
	parser.add_option("-M", "--max", action="store", dest="max", help="maximum genome length", default=7000000, type="int")
	
	parser.add_option_group(group)
	
	return parser.parse_args()



def check_input_options(options, args):

	if options.shuffled and not os.path.isfile(options.shuffled):
		print "Cannot find file", options.shuffled
		sys.exit()
	elif not options.shuffled:
		if not options.forward or not os.path.isfile(options.forward):
			print "Cannot find file", options.forward
			sys.exit()
		if not options.reverse or not os.path.isfile(options.reverse):
			print "Cannot find file",options.reverse
			sys.exit()


	
	if options.maxreads<0:
		options.maxreads=0



#################################
# Function count line in a file #
#################################

def bufcount(filename):
    f = open(filename)                  
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)

    return lines
    





	

def get_read_length(fastqfile):
	if options.readlength>0:
		return options.readlength
	if fastqfile.split(".")[-1]=="gz":
		p = subprocess.Popen(["zcat", fastqfile], stdout = subprocess.PIPE)
		fh = io_method(p.communicate()[0])
		assert p.returncode == 0
		readlength=len(fh.readlines()[-1])
	else:
		readlength=len(os.popen("head -n 2 "+fastqfile).readlines()[-1])
	return readlength
	


    
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
		#pysam.sort( "-n", sys.argv[1], tmpname )
	#	os.system("samtools sort -n "+sys.argv[1]+" "+tmpname )
	#	os.system("cp "+sys.argv[1]+".bai "+tmpname+".bam.bai" )
	#	samfile = pysam.Samfile( tmpname+".bam", "rb" )
	#	os.system("rm "+tmpname+".bam*")
	
	
		samfile = pysam.Samfile( samfilename, "rb" )
	
	
	elif samfilename.split(".")[-1]=="sam":
		samfile = pysam.Samfile( samfilename, "r" )
	else:
		print "Not a sam or bam file"
		sys.exit()	
	
	if shuffled:
		output=open(samfilename.split(".")[0]+"_unmapped.fastq", "w")
	else:
		outputf=open(samfilename.split(".")[0]+"_1_unmapped.fastq", "w")
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
				if shuffled:
					bamline2fastq(read, "1", output)
					bamline2fastq(mate, "2", output)
				else:
					bamline2fastq(read, "1", outputf)
					bamline2fastq(mate, "2", outputr)
					
				addedcount+=1
		
	
	print addedcount, "pairs added to unmapped reads fastq file"
	if shuffled:
		output.close()
	else:
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
	
	check_input_options(options, args)
	
	
	forward=options.forward
	reverse=options.reverse
	shuffled=options.shuffled
	readspersplit=options.readspersplit*2
	maxreads=options.maxreads*2
	outputname=options.outputname
	
	if not shuffled:
		readlength=get_read_length(forward)
		nameprefix=forward.replace("_1.fastq", "")
	else:
		nameprefix=shuffled.replace(".fastq", "")
		readlength=get_read_length(shuffled)
	
	if readspersplit<0:
		print "number of reads to put in each partition must be 0 or greater"
		sys.exit()
	
	
	
	
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	
	if not shuffled:
		print "Shuffling sequences"
		
		os.system("~sh16/scripts/shufflefastqSequences.pl "+forward+" "+reverse+" "+tmpname+".fastq")
	else:
		if shuffled.split(".")[-1]=="gz":
			os.system("zcat "+shuffled+" > "+tmpname+".fastq")
		else:
			os.system("cp "+shuffled+" "+tmpname+".fastq")
		

	filetype="fastq"
#	if options.filter:
	os.system("~sh16/scripts/fastq2fasta.pl "+tmpname+".fastq "+tmpname+".fasta")
#	os.system("export EUSRC=/nfs/users/nfs_s/sh16/euler-sr-ec-2.0.2")
#	os.putenv("EUSRC", "/nfs/users/nfs_s/sh16/euler-sr-ec-2.0.2")
#	os.system("~/euler-sr-ec-2.0.2/assembly/EulerEC.pl "+tmpname+".fasta 55 -minMult 10")
#	os.system("mv fixed/"+tmpname+".fasta "+tmpname+".fasta")
#	
#	os.system("rm -rf fixed/"+tmpname+".fasta")
	filetype="fasta"
		
	
	
	
	
	#Count the number of reads in the files
	
	if filetype=="fastq":
		linecount=bufcount(tmpname+".fastq")
		linesperread=4
		filename=tmpname+".fastq"
	else:
		linecount=bufcount(tmpname+".fasta")
		if readlength>60:
			linesperread=3
		elif readlength>120:
			linesperread=4
		elif readlength<=60:
			linesperread=2
		filename=tmpname+".fasta"
			
	linesperpair=linesperread*2
	
	 
	if readspersplit==0 or maxreads>0:
		numfiles=1
		if maxreads>0 and maxreads<=float(linecount)/linesperpair:
			readsperfile=maxreads/2
		else:
			readsperfile=int(float(linecount)/linesperpair)
	else:
		numfiles=(linecount/(linesperpair/2))/readspersplit
		readsperfile=int((float(linecount)/linesperpair)/(numfiles))
	if numfiles<1:
		numfiles=1
	
	
	
	print "There are", linecount/linesperpair, "paired reads in your fastq file. Splitting fastq into", numfiles, "files containing approximately", readsperfile, "paired reads each"
	sys.stdout.flush()
	
	outfiles=[]
	print maxreads, readsperfile
	if maxreads>0 and readsperfile<float(linecount)/linesperpair:
		outfile = open(tmpname+"_1."+filetype,"w")
		readstoadd=sample(xrange(linecount/linesperpair), readsperfile)
		
		readstoadd.sort()
		readsset=set(readstoadd)
		fastqfile=open(filename,"rU")
		linenum=0
		count=0
		for x in fastqfile:
			lines=[x.strip()]
			linenum+=1
			for y in range(0,linesperpair-1):
				lines.append(fastqfile.next().strip())
				linenum+=1
			
			if linenum/linesperpair in readsset:
				for line in lines:
					count+=1
					print >> outfile, line
				readsset.remove(linenum/linesperpair)
				
		outfile.close()
		print count
					
	else:
		fastqfile=open(filename,"rU")
			
		for x in range(0,numfiles):
			outfiles.append(open(tmpname+"_"+str(x+1)+"."+filetype,"w"))
			
		for x in fastqfile:
		
			lines=[x.strip()]
				
			for y in range(0,linesperpair-1):
				lines.append(fastqfile.next().strip())
			
			rannum=randint(0, len(outfiles)-1)
			
			for line in lines:
				print >> outfiles[rannum], line
			
		
		for output in outfiles:
			output.close()
		
			

	
	expectedcoverage=(float(readlength*(readsperfile))/options.max)/numfiles

	for x in range(0,numfiles):
		if os.path.isfile(tmpname+"contigs"):
			os.system(velvet_dir+"velveth "+tmpname+"_velvet 21,61,10 -"+filetype+" -shortPaired "+tmpname+"_"+str(x+1)+"."+filetype+" -long "+tmpname+"contigs")
			os.system(velvet_dir+"velvetg "+tmpname+"_velvet -ins_length 200 -min_contig_lgth 100 -unused_reads yes -long_mult_cutoff 0 -exp_cov "+str(expectedcoverage)+" -cov_cutoff "+str(0.5*expectedcoverage))
		else:
			os.system(velvet_dir+"velveth "+tmpname+"_velvet 21,61,10 -"+filetype+" -shortPaired "+tmpname+"_"+str(x+1)+"."+filetype)
			os.system(velvet_dir+"velvetg "+tmpname+"_velvet -ins_length 200 -min_contig_lgth 100 -unused_reads yes -exp_cov "+str(expectedcoverage)+" -cov_cutoff "+str(0.5*expectedcoverage))
		
		if os.path.isfile(tmpname+"_velvet/contigs.fa"):
			os.system("cp "+tmpname+"_velvet/contigs.fa "+options.outputname+str(x+1))
			os.system("mv "+tmpname+"_velvet/contigs.fa "+tmpname+"contigs")
			if x<numfiles-1:
				os.system("cat "+tmpname+"_velvet/UnmappableReads.fa >> "+tmpname+"_"+str(x+2)+"."+filetype)
		else:
			print "Velvetg failed"
			sys.exit()
		os.system("rm -rf "+tmpname+"_velvet")

	os.system("mv "+tmpname+"contigs "+options.outputname)
	os.system("rm -rf "+tmpname+"*")
	print tmpname

	
