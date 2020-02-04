#!/usr/bin/env python
import os, sys
from optparse import OptionParser


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

#	parser.add_option("-f", "--forward", action="store", dest="forward", help="forward fastq file", default="", metavar="FILE")
#	parser.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fastq file", default="", metavar="FILE")
	parser.add_option("-s", "--shuffled", action="store", dest="shuffled", help="shuffled fastq/fasta file", default="", metavar="FILE")
	parser.add_option("-f", "--filetype", action="store", dest="filetype", help="file type (fastq or fastq) [default=%default]", type="choice", choices=["fastq", "fasta"], default="fastq", metavar="FILE")
	parser.add_option("-l", "--readlength", action="store", dest="readlength", help="Read length (required for multi-line fasta files)", default=0, type="int", metavar="length")
	parser.add_option("-k", "--kmer", action="store", dest="kmer", help="kmer (0=auto)", default=0, type="int", metavar="kmer")
	parser.add_option("-m", "--max", action="store", dest="max", help="minimum genome length", default=7000000, type="int")
	parser.add_option("-M", "--memory", action="store", dest="memory", help="memory to use in Gb", default=2, type="int")
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.shuffled=="":
		DoError("No tree provided")
	elif not os.path.isfile(options.shuffled):
		DoError("Fastq file does not exist")
		
	return


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


###############################
# Function to get read length #
###############################

def get_read_length(fastqfile):
	if fastqfile.split(".")[-1]=="gz":
		p = subprocess.Popen(["zcat", fastqfile], stdout = subprocess.PIPE)
		fh = io_method(p.communicate()[0])
		assert p.returncode == 0
		readlength=len(fh.readlines()[-1])
	else:
		readlength=len(os.popen("head -n 2 "+fastqfile).readlines()[-1])
	return readlength

velvet_dir="/nfs/users/nfs_s/sh16/velvet-sc/"

################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	if options.readlength<1:
		readlength=get_read_length(options.shuffled)
	else:
		readlength=options.readlength
	
	if options.filetype=="fastq":
		filelength=bufcount(options.shuffled)
	else:
		filelength=0
		for line in open(options.shuffled):
			if len(line)>0 and line[0]==">":
				filelength+=4
	
	expectedcoverage=float(readlength*(filelength/4))/options.max
	
	if options.kmer==0:
		options.kmer=(2*readlength)/3
	if options.kmer>61:
		options.kmer=61
	
	print options.kmer, expectedcoverage
	
	membig=options.memory*1000000
	memsmall=options.memory*1000
	
	
	print velvet_dir+"velveth "+'.'.join(options.shuffled.split("/")[-1].split(".")[:-1])+"_velvet_"+str(options.kmer)+" "+str(options.kmer)+" -"+options.filetype+" -shortPaired "+options.shuffled
	
	os.system(velvet_dir+"velveth "+'.'.join(options.shuffled.split("/")[-1].split(".")[:-1])+"_velvet_"+str(options.kmer)+" "+str(options.kmer)+" -"+options.filetype+" -shortPaired "+options.shuffled)
	
	x=0.5
	for x in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]:
		
		os.system("mkdir "+'.'.join(options.shuffled.split("/")[-1].split(".")[:-1])+"_velvet_"+str(options.kmer)+"_"+str(x))
		os.system("cp -r "+'.'.join(options.shuffled.split("/")[-1].split(".")[:-1])+"_velvet_"+str(options.kmer)+"/* "+'.'.join(options.shuffled.split(".")[:-1])+"_velvet_"+str(options.kmer)+"_"+str(x))
		print "bsub -M "+str(membig)+" -R 'select[mem>"+str(memsmall)+"] rusage[mem="+str(memsmall)+"]' -o "+str(x)+".out -e "+str(x)+".err "+velvet_dir+"velvetg "+".".join(options.shuffled.split("/")[-1].split(".")[:-1])+"_velvet_"+str(options.kmer)+"_"+str(x)+" -ins_length 200 -exp_cov "+str(expectedcoverage)+" -cov_cutoff "+str(x*expectedcoverage)
		os.system("bsub -M "+str(membig)+" -R 'select[mem>"+str(memsmall)+"] rusage[mem="+str(memsmall)+"]' -o "+str(x)+".out -e "+str(x)+".err "+velvet_dir+"velvetg "+".".join(options.shuffled.split("/")[-1].split(".")[:-1])+"_velvet_"+str(options.kmer)+"_"+str(x)+" -ins_length 200 -exp_cov "+str(expectedcoverage)+" -cov_cutoff "+str(x*expectedcoverage))
		
	
	
	