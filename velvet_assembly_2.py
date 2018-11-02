#!/usr/bin/env python
import string, re, numpy
import os, sys
from random import *
from optparse import OptionParser, OptionGroup
import subprocess
import time
from Bio.Seq import Seq
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_SeqIO import *

########################
# Define some globals #
#######################

SAMTOOLS_DIR=""
BWA_DIR=""
velvet_dir=""


def get_read_length(fastqfile):
	if fastqfile.split(".")[-1]=="gz":
		p = subprocess.Popen(["zcat", fastqfile], stdout = subprocess.PIPE)
		fh = io_method(p.communicate()[0])
		assert p.returncode == 0
		readlength=len(fh.readlines()[-1])
	else:
		readlength=len(os.popen("head -n 2 "+fastqfile).readlines()[-1])
	return readlength


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	group = OptionGroup(parser, "Input Options")

	group.add_option("-f", "--forward", action="store", dest="forward", help="forward fastq file", default=False)
	group.add_option("-r", "--reverse", action="store", dest="reverse", help="reverse fastq file", default=False)
	group.add_option("-s", "--shuffled", action="store", dest="shuffled", help="shuffled fastq file", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Output Options")
	group.add_option("-o", "--output", action="store", dest="outputname", help="output file name [default=%default]", default="best_assembly.mfa")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Assembly Options")
	group.add_option("-a", "--optimiser", action="store", type="choice", dest="optimiser", choices=["max", "n50", "total", "coverage"], help="Parameter to use for optimisation (max, n50, total, coverage) [default= %default]", default="max")
	group.add_option("-k", "--min_kmer", action="store", type="int", dest="minkmer", help="Minimum kmer to use [default= %default]", default=19)
	group.add_option("-K", "--max_kmer", action="store", type="int", dest="maxkmer", help="Maximum kmer to use [default= %default]", default=41)
	group.add_option("-e", "--exp_cov", action="store", type="float", dest="ecov", help="Expected coverage to aim for (if coverage chosen as optimiser) [default= %default]", default=20.0)
	parser.add_option_group(group)
	
	
	
	group = OptionGroup(parser, "Contig Filtering Options")
	group.add_option("-L", "--contiglength", action="store", type="int", dest="contig_length_cutoff", help="Contig length cutoff: remove contigs if they are less than this length (0 for no filtering) [default=%default]", default=100)
	group.add_option("-G", "--contiggc", action="store", dest="contiggcfilter", help="Filter contigs by GC (format must be a letter a (above) or b (below) followed by a percentage) [default=%default]", default=False)
	parser.add_option_group(group)
	
	
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
			print "Cannot find file", options.reverse
			sys.exit()
	
	if options.maxkmer%2==0:
		options.maxkmer+=1
		print "Maximum kmer cannot be odd. Increasing by 1"
	
	if options.minkmer%2==0:
		options.minkmer-=1
		print "Minimum kmer cannot be odd. Reducing by 1"
	
	if options.maxkmer>61:
		options.maxkmer=61
		print "Maximum kmer cannot be larger than 61. Resetting to 61"
	elif options.maxkmer<19:
		options.maxkmer=19
		print "Maximum kmer cannot be smaller than 19. Resetting to 19"
	
	if options.minkmer>61:
		options.minkmer=61
		print "Minimum kmer cannot be larger than 61. Resetting to 61"
	elif options.minkmer<19:
		options.minkmer=19
		print "Minimum kmer cannot be larger than 19. Resetting to 19"
	
	if options.minkmer>options.maxkmer:
		print "Minimum kmer must be smaller than maximum kmer"
		sys.exit()
	
	
	




if __name__ == "__main__":

	if sys.version.startswith("3"):
		import io
		io_method = io.BytesIO
	else:
		import cStringIO
		io_method = cStringIO.StringIO
		
	(options, args) = main()
	
	check_input_options(options, args)

	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))

	forward=options.forward
	reverse=options.reverse
	shuffled=options.shuffled
	outputname=options.outputname
	optimiser=options.optimiser
	ecov=options.ecov
	minkmer=options.minkmer
	maxkmer=options.maxkmer
	
	if not shuffled:
		readlength=get_read_length(forward)
		nameprefix=forward.replace("_1.fastq", "")
		print "Shuffling sequences"
		os.system("/nfs/pathogen/sh16_scripts/shufflefastqSequences.pl "+forward+" "+reverse+" "+tmpname+".fastq")
		filename=tmpname+".fastq"
	else:
		nameprefix=shuffled.replace(".fastq", "")
		readlength=get_read_length(shuffled)
		filename=shuffled

	if maxkmer>=readlength-2:
		maxkmer=readlength-2
	
	
	insert=250

	#make the velvet input string - can add to this later
	velvethstring = " ".join(map(str, [velvet_dir+"velveth ", tmpname+"_${LSB_JOBINDEX} ${LSB_JOBINDEX} -fastq -shortPaired", filename, ]))#"> /dev/null 2>&1"]))
	velvetgstring = " ".join(map(str, [velvet_dir+"velvetg ", tmpname+"_${LSB_JOBINDEX} -ins_length", insert, ]))#"> /dev/null 2>&1"]))

	
	os.system('echo \''+velvethstring+'\' | bsub -o '+tmpname+'.out -e '+tmpname+'.err -J "'+tmpname+'_velveth['+str(minkmer)+'-'+str(maxkmer)+':2]" > '+tmpname+'jobstring')
	
	os.system('echo \''+velvetgstring+'\' | bsub -o '+tmpname+'.out -e '+tmpname+'.err -w \'ended('+tmpname+'_velveth)\' -J "'+tmpname+'_velvetg['+str(minkmer)+'-'+str(maxkmer)+':2]" > '+tmpname+'jobstring2' )
	
		
	jobnum=open(tmpname+'jobstring', "rU").read().split(">")[0].split("<")[1]
	jobbnum=open(tmpname+'jobstring2', "rU").read().split(">")[0].split("<")[1]
	
	os.system("rm -f "+tmpname+'jobstring '+tmpname+'jobstring2')

	todo=1
	print "JOBID    ARRAY_SPEC  OWNER  NJOBS PEND DONE  RUN EXIT SSUSP USUSP PSUSP"
	#while os.path.isfile(tmpname+'waitfile'):
	while todo!=0:
		time.sleep(5)
		bjoblines=os.popen("bjobs -A "+jobnum).readlines()
		if len(bjoblines)<2:
			continue
		bjobsstring=bjoblines[1]
		print bjobsstring.strip()+"\r",
		sys.stdout.flush()
		pend=int(bjobsstring.split()[4])
		done=int(bjobsstring.split()[5])
		run=int(bjobsstring.split()[6])
		exited=int(bjobsstring.split()[7])
		todo=run+pend
	
	print "\nFinished velveth."
	
	todo=1
	print "JOBID    ARRAY_SPEC  OWNER  NJOBS PEND DONE  RUN EXIT SSUSP USUSP PSUSP"
	#while os.path.isfile(tmpname+'waitfile'):
	while todo!=0:
		time.sleep(5)
		bjoblines=os.popen("bjobs -A "+jobbnum).readlines()
		if len(bjoblines)<2:
			continue
		bjobsstring=bjoblines[1]
		print bjobsstring.strip()+"\r",
		sys.stdout.flush()
		pend=int(bjobsstring.split()[4])
		done=int(bjobsstring.split()[5])
		run=int(bjobsstring.split()[6])
		exited=int(bjobsstring.split()[7])
		todo=run+pend
	
	print "\nFinished velvetg."
	
	
	
	
	optimised=[]
	
	for kmer in range(minkmer, maxkmer+1,2):
	
		filename=tmpname+"_"+str(kmer)+"/contigs.fa"

		try:
			contigs=SeqIO.parse(open(filename), "fasta")
		except StandardError:
			print "Could not open file", filename
			continue
		
		lengths=[]
		GCs=[]

		for contig in contigs:

			
			seq=str(contig.seq)
			length=len(seq)
			
			lengths.append(length)
			
			#lengthnons=len(seq.upper().replace("N",""))

		
		if options.optimiser=="max":
			optimised.append([numpy.max(lengths),kmer])
		elif options.optimiser=="max":
			optimised.append([numpy.sum(lengths),kmer])
		elif options.optimiser=="n50":
			lengths.sort()
			lengths.reverse()
			fifty=float(numpy.sum(lengths))/2
			count=0
			sum=0
			while sum<fifty:
				N50=lengths[count]
				sum+=lengths[count]
				count+=1
			optimised.append([N50,kmer])
		else:
			mincov=float(options.ecov)/2
			if mincov<3:
				mincov=3
			mediancov=int(os.popen("/nfs/pathogen/sh16_scripts/velvet_stats_2_av_cov.py "+tmpname+"_"+str(kmer)+"/stats.txt"+" "+str(mincov)).read().strip())
			mediff=options.ecov-mediancov
			if mediff<0:
				mediff=1-mediff
			optimised.append([mediff,kmer])
			
	
	optimised.sort()
	optimised.reverse()
	
	print "\nBest Kmer is", optimised[0][1], ":\n"
	
	print '\t'.join(["RANK", "KMER", options.optimiser.upper()])
	for x, kmer in enumerate(optimised):
		print '\t'.join(map(str,[x+1, kmer[1], kmer[0]]))
	
	mincov=float(options.ecov)/2
	if mincov<3:
		mincov=3
	mediancov=int(os.popen("/nfs/pathogen/sh16_scripts/velvet_stats_2_av_cov.py "+tmpname+"_"+str(optimised[0][1])+"/stats.txt"+" "+str(mincov)).read().strip())
	
	print mediancov
	if mediancov<6:
		cutoff=float(mediancov)/2
	else:
		cutoff=3
		
	
	
	os.system(" ".join(map(str, [velvet_dir+"velvetg ", tmpname+"_"+optimised[0][1], "-ins_length", insert, "-exp_cov", mediancoverage, "-cov_cutoff", cutoff])))
	
	os.system("mv "+tmpname+"_"+optimised[0][1]+"/contigs.fa "+options.outputname) 
	
	#os.system("rm -rf "+tmpname+"_*/Graph "+tmpname+"_*/LastGraph "+tmpname+"_*/Log "+tmpname+"_*/PreGraph "+tmpname+"_*/Roadmaps "+tmpname+"_*/Sequences")
	os.system("rm -rf "+tmpname+"*")

	
		
		
	
		