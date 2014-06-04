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
import gzip
import mimetypes


def cleanup():
	if forward!=options.forward:
		os.system("rm -f "+forward)
	if reverse!=options.reverse:
		os.system("rm -f "+reverse)

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)


	parser.add_option("-f", "--forward", action="store", dest="forward", help="forward fastq file", default=False)
	parser.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fastq file", default=False)
	parser.add_option("-o", "--output", action="store", dest="output", help="output file prefix", default=False)
	parser.add_option("-p", "--proportion", action="store", dest="proportion", help="proportion of reads to keep [default=%default]", type="float", default=0.0)
	parser.add_option("-n", "--number", action="store", dest="number", help="number of reads to keep [default=%default]", type="int", default=2000000)
	parser.add_option("-c", "--coverage", action="store", dest="coverage", help="coverage to aim for. This will override the number of reads option. Note, a reference length must be specified. [default=%default]", type="int", default=0)
	parser.add_option("-l", "--length", action="store", dest="length", help="expected length of genome (required if coverage is being used for read count) [default=%default]", type="int", default=0)
	parser.add_option("-z", "--zip", action="store_true", dest="zip", help="Gzip output files", default=False)
	
	
	return parser.parse_args()



def check_input_options(options, args):

	if not options.forward or not os.path.isfile(options.forward):
		print "Cannot find file", options.forward
		sys.exit()
	if not options.reverse or not os.path.isfile(options.reverse):
		print "Cannot find file",options.reverse
		sys.exit()
	
	if not options.output:
		print "No output file prefix specified"
		sys.exit()
	
	if options.number<=0:
		print "Number of reads must be >0"

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

if __name__ == "__main__":

	if sys.version.startswith("3"):
		import io
		io_method = io.BytesIO
	else:
		import cStringIO
		io_method = cStringIO.StringIO


	(options, args) = main()
	
	check_input_options(options, args)
	
#	if options.forward.split('.')[-1]=="gz":
#		print "Unzipping forward fastq file"
#		os.system("zcat "+options.forward+" > "+'.'.join(options.forward.split("/")[-1].split(".")[:-1]))
#		forward='.'.join(options.forward.split("/")[-1].split(".")[:-1])
#	else:
#		forward=options.forward
#	if options.reverse.split('.')[-1]=="gz":
#		print "Unzipping reverse fastq file"
#		os.system("zcat "+options.reverse+" > "+'.'.join(options.reverse.split("/")[-1].split(".")[:-1]))
#		reverse='.'.join(options.reverse.split("/")[-1].split(".")[:-1])
#	else:
#		reverse=options.reverse
	
	forward=options.forward
	reverse=options.reverse
	
	if mimetypes.guess_type(forward)[1]=="gzip":
		linecount1=0
		for line in gzip.open(forward,"r"):
			linecount1+=1
	else:
		linecount1=bufcount(forward)
	if mimetypes.guess_type(reverse)[1]=="gzip":
		linecount2=0
		for line in gzip.open(reverse,"r"):
			linecount2+=1
	else:
		linecount2=bufcount(reverse)
	if linecount1!=linecount2:
		cleanup()
		DoError("Forward and reverse line counts are not the same")
	linecount=float(linecount1)
	
	
	if options.coverage>0:
		if options.length<=0:
			DoError("When using coverage to subsample, the reference length must be specified.")
		if mimetypes.guess_type(forward)[1]=="gzip":
			fastqfile=gzip.open(forward,"r")
		else:
			fastqfile=open(forward,"rU")
		lines=[]
		for y in range(0,2):
			lines.append(fastqfile.next().strip())
		readlen=len(lines[1])
		subsample=int(float(options.coverage*options.length)/(readlen*2))
		fastqfile.close()
	elif options.proportion>0:
		subsample=int(options.proportion*(linecount/4))
	else:
		subsample=options.number
	
	
	if int(subsample)>=linecount/4:
		cleanup()
		print "Subsample is >= number of reads in the input files. No need to subsample. Exiting..."
		if mimetypes.guess_type(forward)[1]=="gzip":
			os.system("cp "+forward+" "+options.output+"_1.fastq.gz")
		else:
			os.system("cp "+forward+" "+options.output+"_1.fastq")
		if mimetypes.guess_type(reverse)[1]=="gzip":
			os.system("cp "+reverse+" "+options.output+"_2.fastq.gz")
		else:
			os.system("cp "+reverse+" "+options.output+"_2.fastq")
		sys.exit()
	
	if options.zip:
		outfilef=gzip.open(options.output+"_1.fastq.gz","w")
		outfiler=gzip.open(options.output+"_2.fastq.gz","w")
	else:
		outfilef=open(options.output+"_1.fastq","w")
		outfiler=open(options.output+"_2.fastq","w")
	readstoadd=sample(xrange(1, int(linecount/4)+1), subsample)
	readstoadd.sort()
	readstoadd.reverse()
	#readsset=set(readstoadd)
	
	
	print "Sampling", subsample, "paired reads from each fastq file"
	
	linenum=0
	if mimetypes.guess_type(forward)[1]=="gzip":
		fastqfilef=gzip.open(forward,"r")
	else:
		fastqfilef=open(forward,"rU")
	if mimetypes.guess_type(reverse)[1]=="gzip":
		fastqfiler=gzip.open(reverse,"r")
	else:
		fastqfiler=open(reverse,"rU")
	
	for x in fastqfilef:
		linesr=[]
		linesr.append(fastqfiler.next().strip())
		linesf=[x.strip()]
		
		linenum+=1
		for y in range(0,3):
			linesf.append(fastqfilef.next().strip())
			linesr.append(fastqfiler.next().strip())
			#linenum+=1
		
		
		if linenum==readstoadd[-1]:
			for line in linesf:
				print >> outfilef, line
			for line in linesr:
				print >> outfiler, line
			readstoadd.pop()
			
		if len(readstoadd)==0:
			break
	fastqfilef.close()
	fastqfiler.close()
	outfilef.close()
	outfiler.close()
#	if options.zip:
#		print "Zipping forward output file"
#		os.system("gzip "+options.output+"_1.fastq")
#		print "Zipping reverse output file"
#		os.system("gzip "+options.output+"_2.fastq")
	cleanup()
	
	
	
