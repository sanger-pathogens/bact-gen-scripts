#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re
import os, sys, getopt, random, math, time
from random import *



#########
# Usage #
#########

def Usage():
	print '\nrun_multiple_mappings.py Usage:'
	print '\nrun_multiple_mappings.py [OPTIONS] <list of fastq files>'
	print "\nINPUT OPTIONS:"
	print "-r Reference dna sequence:\t<fasta file> (REQUIRED)"
	print "\nMAPPING OPTIONS:"
	print "-P mapping program to use:\t<maq/ssaha> (REQUIRED)"
	print "\nMAQ MAPPING OPTIONS:"
	print "-n max no. of mismatches for mapping:\t\t<integer between 1 and 3>"
	print "-m max no. of snps per read for consensus:\t<integer between 1 and 100>"
	print "-p use paired-end reads"
	print "-s reads are multiplexed with stupid naming convention e.g. 2610_4_1_10 and 2610_4_2_10 for pairs of strain 10 in pool 4"
	print "-i max insert size:\t\t\t\t<integer between 100 and 10,000> [default=300]"
	print "-R read length:\t\t\t\t\t<integer between 36 and 1000> [default=54]"
	print "-q minimum mapping quality:\t\t\t<integer between 1 and 60> [default=30]"	
	print "-d minimum mapping depth:\t\t\t<integer between 1 and 100,000> [default=5]"
	print "\nSsaha MAPPING OPTIONS:"
	print "-t data type:\t\t\t\t\t<454/solexa>"
	print "-p use paired-end reads"
	print "-i max insert size:\t\t\t\t<integer between 10 and 10,000> [default=300]"
	print "-j min insert size:\t\t\t\t<integer between 10 and 10,000> [default=100]"
	print "-R read length:\t\t\t\t\t<integer between 36 and 1000> [default=54]"
	print "-q minimum mapping score:\t\t\t<integer between 1 and 100> [default=30]"
	print "\nUSAGE OPTIONS:"
	print "-h show this help"
	print "-I Turn off interactive phylip-style menu"
	print '\nCopyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'

	

##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hd:r:i:j:pq:vIm:n:P:t:R:e:s", ["depth=", "help", "ref=", "align", "quality=", "program=", "velvet=", "interactive", "maxsnps=", "pairedend", "maxinsert=", "mininsert=", "mismatches=", "rtype=", "readlength=", "embl=", "stupid"])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	ref=''
	inputdirs=[]
	quality=30
	mindepth=5
	program=''
	velvet='n'
	interactive='y'
	maxsnps=2
	pairedend='n'
	maxinsertsize=300
	mininsertsize=100
	mismatches=2
	rtype='solexa'
	readlength=54
	embl=''
	multiplexnamed='n'

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
		elif opt in ("-e", "--embl"):
			embl=arg
		elif opt in ("-q", "--quality"):
			quality=int(arg)
		elif opt in ("-d", "--depth"):
			mindepth=int(arg)
		elif opt in ("-m", "--maxsnps"):
			maxsnps=int(arg)
		elif opt in ("-P", "--program"):
			program=arg.lower()
		elif opt in ("-v", "--velvet"):
			velvet='y'
		elif opt in ("-I", "--interactive"):
			interactive='n'
		elif opt in ("-p", "--pairedend"):
			pairedend='y'
		elif opt in ("-s", "--stupid"):
			multiplexnamed='y'
		elif opt in ("-i", "--maxinsert"):
			maxinsertsize=int(arg)
		elif opt in ("-j", "--mininsert"):
			mininsertsize=int(arg)
		elif opt in ("-n", "--mismatches"):
			mismatches=int(arg)
		elif opt in ("-t", "--rtype"):
			rtype=arg.lower()
		elif opt in ("-R", "--readlength"):
			readlength=int(arg)
		
	inputdirs=args

	if ref=='':
		print 'Error: No reference dna file (-r) selected!'
		Usage()
		sys.exit()	
	elif inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	elif rtype not in ['454', 'solexa']:
		print 'Error: Data type (-t) must be 454 or solexa!'
		Usage()
		sys.exit()
	elif program not in ['ssaha', 'maq']:
		print 'Error: Program (-P) must be maq or ssaha!'
		Usage()
		sys.exit()
	elif maxinsertsize<=mininsertsize:
		print 'Error: Minimum insert size (-j) must be smaller than maximum (-i)!'
		Usage()
		sys.exit()
	elif mininsertsize>10000 or mininsertsize<10:
		print 'Error: Minimum insert size (-j) must be between 10 and 10,000!'
		Usage()
		sys.exit()
	elif quality>100 or quality<1:
		print 'Error: Mapping quality score (-q) must be between 1 and 100!'
		Usage()
		sys.exit()
	elif program=='maq' and (mindepth>100000 or mindepth<1):
		print 'Error: Depth value (-d) must be between 1 and 100000!'
		Usage()
		sys.exit()
	elif program=='maq' and (maxsnps>100 or maxsnps<1):
		print 'Error: Maximum number of snps per read (-m) for consensus must be between 1 and 100!'
		Usage()
		sys.exit()
	elif program=='maq' and (mismatches>3 or mismatches<1):
		print 'Error: Maximum number of mismatches for mapping (-n) 1 and 3!'
		Usage()
		sys.exit()
	elif readlength>1000 or readlength<36:
		print 'Error: Read length must be between (-R) 36 and 1000!'
		Usage()
		sys.exit()
	elif program=='maq' and rtype=='454':
		print 'Error: Only ssaha can use 454 reads!'
		Usage()
		sys.exit()

	snpquality=quality
	snpdepth=mindepth
	raxml='n'
	bootstrap=0
	tabfile='n'
	align='n'
	graphs='n'
	
	if interactive=='y':
		ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed=menusystem(ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed)
	
	
	outfile=ref.split('.')[0]+"_q"+str(snpquality)+"_d"+str(snpdepth)
	model="GTRGAMMA"
	

	#summarymenu(outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, snpquality, snpdepth)
	
	
	return ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed



####################
# Interactive Menu #
####################

def menusystem(ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed):
	
	os.system('clear')
	
	
	print "\nrun_multiple_mappings.py: Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009"
	
	print "\nINPUT OPTIONS:"
	
	if ref=='':
		print "r: Reference dna sequence:\t\tNone selected (required)"
	else:
		print "r: Reference dna sequence:\t\t"+ref
		if len(inputdirs)==1:
			print len(inputdirs), "file to be mapped"
		else:
			print len(inputdirs), "files to be mapped"
		print "\nMAPPING OPTIONS:"
		print "P: Program:\t\t\t\t"+program
		if program=='ssaha':
			print "t: Read type:\t\t\t\t"+rtype
		if rtype=='solexa':
			print "R: Read length:\t\t\t\t"+str(readlength)
			
		if pairedend=='n':
			print "p: Use paired-end reads:\t\tno"
		else:
			print "p: Use paired-end reads:\t\tyes"
			if program=='maq':
				print "i: Maximum insert size:\t\t\t"+str(maxinsertsize)
				
			elif program=='ssaha':
				print "i: Maximum insert size:\t\t\t"+str(maxinsertsize)
				print "j: Minimum insert size:\t\t\t"+str(mininsertsize)
				
		if program=='maq':
			print "n: Maximum mismatches:\t\t\t"+str(mismatches)
			print "m: Maximum SNPs per read:\t\t"+str(maxsnps)
			print "d: Minimum mapping depth:\t\t"+str(mindepth)
			print "q: Minimum mapping quality:\t\t"+str(quality)
		elif program=='ssaha':
			print "q: Minimum mapping score:\t\t"+str(quality)
		
		if velvet=='n':
			print "v: Assemble non-mapping reads:\t\tno"
		else:
			print "v: Assemble non-mapping reads:\t\tyes"
		print "D: Reset to default mapping parameters"
		
	print "\nQ: QUIT"
	
	if ref=="":
		message="\nPlease select an option:"
		inputlist=['r', 'Q']
	else:
		message="\nPlease select an option or type y to run:"
		if program=='ssaha':
			inputlist=['r','q','v','y','t','Q','D','p','R','P']
		else:
			inputlist=['r','q','d','m','n','v','y','Q','D','p','P', 'R']
		if pairedend=='y':
			inputlist=inputlist+['i']
			if program=='ssaha':
				inputlist=inputlist+['j']
			
	ui=''
	while ui not in inputlist:
		ui=raw_input(message+' ')

	if ui=='y':
		os.system('clear')
		return ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed
		
	elif ui=='r':
		oldref=ref
		ref=''
		while not os.path.isfile(ref):
			ref=raw_input('Enter reference file name including path or Q to go back to the menu: ')
			if ref=='Q':
				ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed=menusystem(oldref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed)
			elif not os.path.isfile(ref):
				print "File not found"
				
	
	elif ui=='n' and program=='maq':
		mismatches=0
		while mismatches > 3 or mismatches < 1:
			mismatches=int(raw_input('Enter maximum mismatches (1-3): '))
		
	elif ui=='m' and program=='maq':
		maxsnps=0
		while maxsnps > 100 or maxsnps < 1:
			maxsnps=int(raw_input('Enter maximum number of SNPs per read (1-100): '))
	
	elif ui=='q':
		quality=0
		while quality > 100 or quality < 1:
			quality=int(raw_input('Enter minimum mapping quality (1-100): '))
	
	elif ui=='d' and program=='maq':
		mindepth=0
		while mindepth > 100000 or mindepth < 1:
			mindepth=int(raw_input('Enter minimum read depth for mapping and SNP calling (1-100,000): '))
	
	elif ui=='v':
		if velvet=='n':
			velvet='y'
		else:
			velvet='n'
	
#	elif ui=='o':
#		outfile=''
#		while outfile=='':
#			outfile=raw_input('Enter prefix for output file names, or D to use the default: ')
#			if outfile!='D' and outfile!='':
#				outorig='n'
#			elif outfile=='D':
#				outorig='y'
	
	elif ui=='P':
		if program=='ssaha':
			program='maq'
		else:
			program='ssaha'
	
	elif ui=='t':
		if rtype=='solexa':
			rtype='454'
		else:
			rtype='solexa'
	
	elif ui=='D':
		quality=30
		mindepth=5
		maxsnps=2
	
	elif ui=='p':
		if pairedend=='n':
			pairedend='y'
		else:
			pairedend='n'
	
	elif ui=='i':
		maxinsertsize=0
		while maxinsertsize > 10000 or maxinsertsize < 10 or maxinsertsize<=mininsertsize:
			maxinsertsize=int(raw_input('Enter maximum insert size (10-10,000). Must be more than min: '))
	
	elif ui=='j' and program=='ssaha':
		mininsertsize=0
		while mininsertsize > 10000 or mininsertsize < 10 or maxinsertsize<=mininsertsize:
			mininsertsize=int(raw_input('Enter minimum insert size (10-10,000). Must be less than max: '))
	
	elif ui=='R' and rtype=='solexa':
		readlength=0
		while readlength > 1000 or readlength < 10:
			readlength=int(raw_input('Enter read length (10-1000): '))
			
	elif ui=='Q':
		sys.exit()
	
#	if outorig=='y':
#		outfile=ref.split('.')[0]+"_q"+str(quality)+"_d"+str(mindepth)
		
	ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed=menusystem(ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed)

	return ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed






#####################
# SNPanalysis class #
#####################

class SNPanalysis:
	def __init__(self, fastq='', name='', mapped={}, runssaha='n', CDSseq='', number=0):
		self.fastq=fastq
		self.name=name
		self.runname=''
		self.quality=30
		self.rtype=rtype
		self.mindepth=5
		self.maxsnps=2
		self.mismatches=2
		self.mininsertsize=300
		self.maxinsertsize=300
		self.pairedend='n'
		self.readlength=54
		self.fastqdir=''
		self.number=number
	def runMaq(self, ref, tmpname):
		print "\nRunning Maq on "+self.name+'...',
		sys.stdout.flush()
		bashfile=open(self.number+tmpname+'_mbs.sh','w')
		
		print >> bashfile, "maq fasta2bfa "+ref+" "+self.runname+"/ref.bfa"
		
		maxinsert=(self.readlength*2)+self.maxinsertsize
		
		if self.pairedend=='n':
			print "maq fastq2bfq "+self.fastqdir+self.name+".fastq "+self.runname+"/reads_1.bfq"
			print >> bashfile, "maq fastq2bfq "+self.fastqdir+self.name+".fastq "+self.runname+"/reads_1.bfq"
			print >> bashfile,  "maq map -n "+str(self.mismatches)+" -u "+self.runname+"/unmap_1.txt "+self.runname+"/reads_1.map "+self.runname+"/ref.bfa  "+self.runname+"/reads_1.bfq"
			print >> bashfile, "maq mapcheck -s -m "+str(self.maxsnps)+" -q "+str(self.quality)+" "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map >"+self.runname+"/mapcheck.txt"
			print >> bashfile, "maq assemble -s -m "+str(self.maxsnps)+" -q "+str(self.quality)+" "+self.runname+"/consensus.cns "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map"
			print >> bashfile, "maq cns2fq "+self.runname+"/consensus.cns >"+self.runname+"/cns.fq"
			print >> bashfile, "maq cns2snp "+self.runname+"/consensus.cns >"+self.runname+"/cns.snp"
			print >> bashfile, "maq.pl SNPfilter -a -d "+str(self.mindepth)+" -q 20 "+self.runname+"/cns.snp >"+self.runname+"/cns.final.snp"
			print >> bashfile, "maq pileup -s -m "+str(self.maxsnps)+" -q "+str(self.quality)+" -v "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map > "+self.runname+"/all.pileup"
		else:
			
			if multiplexnamed=='y':
				print >> bashfile, "maq fastq2bfq "+self.fastqdir+'_'.join(self.name.split('_')[:2])+"_1_"+self.name.split('_')[-1]+".fastq "+self.runname+"/reads_1.bfq"
				print >> bashfile, "maq fastq2bfq "+self.fastqdir+'_'.join(self.name.split('_')[:2])+"_2_"+self.name.split('_')[-1]+".fastq "+self.runname+"/reads_2.bfq"
			else:
				print >> bashfile, "maq fastq2bfq "+self.fastqdir+self.name+"_1.fastq "+self.runname+"/reads_1.bfq"
				print >> bashfile, "maq fastq2bfq "+self.fastqdir+self.name+"_2.fastq "+self.runname+"/reads_2.bfq"
			print >> bashfile, "maq map -n "+str(self.mismatches)+" -a "+str(maxinsert)+" -u "+self.runname+"/unmap_1.txt "+self.runname+"/reads_1.map "+self.runname+"/ref.bfa  "+self.runname+"/reads_1.bfq "+self.runname+"/reads_2.bfq"
			print >> bashfile, "maq mapcheck -m "+str(self.maxsnps)+" -q "+str(self.quality)+" "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map > "+self.runname+"/mapcheck.txt"
			print >> bashfile, "maq assemble -m "+str(self.maxsnps)+" -q "+str(self.quality)+" "+self.runname+"/consensus.cns "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map"
			print >> bashfile, "maq cns2fq "+self.runname+"/consensus.cns > "+self.runname+"/cns.fq"
			print >> bashfile, "maq cns2snp "+self.runname+"/consensus.cns > "+self.runname+"/all.snp"
			print >> bashfile, "maq.pl SNPfilter -d "+str(self.mindepth)+" -q 20 "+self.runname+"/all.snp > "+self.runname+"/cns.final.snp"
			print >> bashfile, "maq pileup -m "+str(self.maxsnps)+" -v "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map > "+self.runname+"/all_cov.pileup"
			print >> bashfile, "awk '{print $4}' "+self.runname+"/all_cov.pileup > "+self.runname+"/random_exact_repeats_coverage.plot"
			print >> bashfile, "maq pileup -m "+str(self.maxsnps)+" -q "+str(self.quality)+" -v "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map > "+self.runname+"/all.pileup"
		print >> bashfile, "rm "+self.runname+"/*.bfq "+self.runname+"/*.bfa "+self.runname+"/all_cov.pileup"
		print >> bashfile, "gzip -f "+self.runname+"/*"
		bashfile.close()
	
	
	def runSsaha(self, ref, tmpname):	
		print "\nRunning Ssaha on "+self.name+'...',
		sys.stdout.flush()
		meaninsert=((self.maxinsertsize-self.mininsertsize)/2)+self.mininsertsize
		bashfile=open(self.number+tmpname+'_sbs.sh','w')
		#Ssaha commands.
		
		#single end mapping
		if pairedend=='n':
			if rtype=='454':
				#extract fasta and fasta.qual from sff using linker
				print >> bashfile, "~sh16/mira_2.9.37_dev_linux-gnu_x86_64/3rdparty/sff_extract -c -o "+self.runname+"/extract.tmp "+self.fastqdir+self.name+".sff"

				print >> bashfile, "PERL5LIB=$PERL5LIB:~tdo/bin/oldPerl"

				print >> bashfile, "export PERL5LIB"
				#turn fasta+fasta.qual to fastq
				print >> bashfile, "perl -w -e \"use AssemblyTools_unstable;AssemblyTools_unstable::fasta2fastq(  '"+self.runname+"/extract.tmp.fasta', '"+self.runname+"/shuffled.tmp' );\""		
			
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -seeds 5 -score "+str(self.quality)+" -kmer 13 -skip 4 -diff 0 -output cigar "+ref+" "+self.runname+"/shuffled.tmp > "+self.runname+"/cigar1.tmp"
			else:
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -score "+str(self.quality)+" -skip 2 -diff 0 -output cigar "+ref+" "+self.fastqdir+self.name+".fastq > "+self.runname+"/cigar1.tmp"
			print >> bashfile, 'grep "^cigar" '+self.runname+"/cigar1.tmp > "+self.runname+"/cigar2.tmp"
			print >> bashfile, "mv "+self.runname+"/cigar2.tmp "+self.runname+"/cigar1.tmp"
			
			print >> bashfile, "perl /nfs/users/nfs_s/sh16/scripts/cigar00_2plot.pl "+self.runname+"/cigar1.tmp "+self.runname+"/allcoverage.plot "+ str(self.readlength)
			print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_cigar "+self.runname+"/cigar1.tmp "+self.runname+"/cigar2.tmp"
			print >> bashfile, "awk '{print $2}' "+self.runname+"/cigar2.tmp > "+self.runname+"/readnames.tmp"
			
			print >> bashfile, "/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/other_codes/get_seqreads/get_seqreads "+self.runname+"/readnames.tmp "+self.fastqdir+self.name+".fastq "+self.runname+"/fastq.tmp"
			
			
			
			if rtype=='454':
				print >> bashfile, "/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/other_codes/get_seqreads/get_seqreads "+self.runname+"/readnames.tmp "+self.runname+"/shuffled.tmp "+self.runname+"/fastq.tmp"
				print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/exclude/get_excreads "+self.runname+"/readnames.tmp "+self.runname+"/shuffled.tmp "+self.name+"unmap.fastq"
				print >> bashfile, "mv "+self.name+"unmap.fastq "+self.runname+"/unmap.fastq"
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 0 -trans 0 "+self.runname+"/cigar2.tmp "+ref+" "+self.runname+"/fastq.tmp > "+self.runname+"/snp.tmp"
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 0 -trans 0 "+self.runname+"/cigar2.tmp "+ref+" "+self.runname+"/fastq.tmp > "+self.runname+"/pileup.tmp"
			elif rtype=='solexa':
				print >> bashfile, "/nfs/users/nfs_t/tdo/bin/pileup_v0.4/ssaha_pileup/other_codes/get_seqreads/get_seqreads "+self.runname+"/readnames.tmp "+self.fastqdir+self.name+".fastq "+self.runname+"/fastq.tmp"
				if velvet=='y':
					print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/exclude/get_excreads "+self.runname+"/readnames.tmp "+self.fastqdir+self.name+".fastq "+self.name+"unmap.fastq"
					print >> bashfile, "mv "+self.name+"unmap.fastq "+self.runname+"/unmap.fastq"
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 1 -trans 0 "+self.runname+"/cigar2.tmp "+ref+" "+self.runname+"/fastq.tmp > "+self.runname+"/snp.tmp"
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 1 -trans 0 "+self.runname+"/cigar2.tmp "+ref+" "+self.runname+"/fastq.tmp > "+self.runname+"/pileup.tmp"
			print >> bashfile, "egrep ^SNP "+self.runname+"/snp.tmp > "+self.runname+"/all.snp"
			print >> bashfile, "egrep ^cons "+self.runname+"/pileup.tmp > "+self.runname+"/all.pileup"
			
			
		#paired end
		else:
			if rtype=='454':
			
				#extract fasta and fasta.qual from sff using linker
				#titanium linker
				print >> bashfile, "~sh16/mira_2.9.37_dev_linux-gnu_x86_64/3rdparty/sff_extract -c -o "+self.runname+"/extract.tmp -l ~tdo/work/linker.fasta -i \"insert_size:"+str(meaninsert)+",insert_stdev:"+str(meaninsert-mininsertsize)+"\" "+self.fastqdir+self.name+".sff"
				#older linker
				#print >> bashfile, "~sh16/mira_2.9.37_dev_linux-gnu_x86_64/3rdparty/sff_extract -c -o "+self.runname+"/extract.tmp -l /nfs/pathdata/Salmonella/typhimurium/DT56/linker.txt -i \"insert_size:"+str(meaninsert)+",insert_stdev:"+str(meaninsert-mininsertsize)+"\" "+self.fastqdir+self.name+".sff"
				
				print >> bashfile, "PERL5LIB=$PERL5LIB:~tdo/bin/oldPerl"

				print >> bashfile, "export PERL5LIB"
				#turn fasta+fasta.qual to fastq
				print >> bashfile, "perl -w -e \"use AssemblyTools_unstable;AssemblyTools_unstable::fasta2fastq(  '"+self.runname+"/extract.tmp.fasta', '"+self.runname+"/shuffled.tmp' );\""
				
				#
				
				print >> bashfile, "echo Extracting reads"
				print >> bashfile, "~sh16/scripts/sff_extract_fastq_splitter.py "+self.runname+"/shuffled.tmp"
				print >> bashfile, "cat "+self.runname+"/shuffled.forward.tmp "+self.runname+"/shuffled.reverse.tmp > "+self.runname+"/shuffled.tmp"

				#do I need to revcomp forward reads? yes...done
				
				#map single reads
				#print >> bashfile, "echo Mapping single end reads"
				#print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -seeds 5 -score "+str(self.quality)+" -kmer 13 -skip 4 -diff 0 -output cigar "+ref+" "+self.runname+"/shuffled.single.tmp > "+self.runname+"/cigarsingle.tmp"
				
				#map paired reads
				print >> bashfile, "echo Mapping paired reads"
				#print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -pair "+str(self.mininsertsize)+','+str(self.maxinsertsize)+" -seeds 5 -score "+str(self.quality)+" -kmer 13 -skip 4 -diff 0 -output cigar "+ref+" "+self.runname+"/shuffled.forward.tmp "+self.runname+"/shuffled.reverse.tmp > "+self.runname+"/cigarpair.tmp"
				#print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -score "+str(self.quality)+" -kmer 13 -skip 4 -diff 0 -output cigar "+ref+" "+self.runname+"/shuffled.paired.tmp > "+self.runname+"/cigarpair.tmp"
				#print >> bashfile, "cat "+self.runname+"/cigarpair.tmp "+self.runname+"/cigarsingle.tmp > "+self.runname+"/cigar1.tmp"
				#print >> bashfile, "rm "+self.runname+"/extract.tmp* "+self.runname+"/cigarpair.tmp "+self.runname+"/cigarsingle.tmp"
				#print >> bashfile, "cat "+self.runname+"/shuffled.paired.tmp "+self.runname+"/shuffled.single.tmp > "+self.runname+"/shuffled.tmp"
				#print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -score "+str(self.quality)+" -kmer 13 -skip 4 -diff 0 -output cigar "+ref+" "+self.runname+"/shuffled.tmp > "+self.runname+"/cigar1.tmp"
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -score "+str(self.quality)+" -skip 2 -diff 0 -kmer 13 -outfile "+self.runname+"/abnormal.tmp -pair "+str(self.mininsertsize)+","+str(self.maxinsertsize)+" -output cigar "+ref+" "+self.runname+"/shuffled.forward.tmp "+self.runname+"/shuffled.reverse.tmp > "+self.runname+"/cigar1.tmp"
				print >> bashfile, "rm "+self.runname+"/extract.tmp*"
				
				
			else:
				
				if multiplexnamed=='y':
					
					print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -score "+str(self.quality)+" -skip 2 -diff 0 -kmer 13 -outfile "+self.runname+"/abnormal.tmp -pair "+str(self.mininsertsize)+","+str(self.maxinsertsize)+" -output cigar "+ref+" "+self.fastqdir+'_'.join(self.name.split('_')[:2])+"_1_"+self.name.split('_')[-1]+".fastq "+self.fastqdir+'_'.join(self.name.split('_')[:2])+"_2_"+self.name.split('_')[-1]+".fastq > "+self.runname+"/cigar1.tmp"
					
				else:
				
					print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -score "+str(self.quality)+" -skip 2 -diff 0 -kmer 13 -outfile "+self.runname+"/abnormal.tmp -pair "+str(self.mininsertsize)+","+str(self.maxinsertsize)+" -output cigar "+ref+" "+self.fastqdir+self.name+"_1.fastq "+self.fastqdir+self.name+"_2.fastq > "+self.runname+"/cigar1.tmp"

			print >> bashfile, 'grep "^cigar" '+self.runname+"/abnormal.tmp > "+self.runname+"/abnormal.cigar"
			print >> bashfile, 'grep "^cigar" '+self.runname+"/cigar1.tmp > "+self.runname+"/cigar2.tmp"
			print >> bashfile, "mv "+self.runname+"/cigar2.tmp "+self.runname+"/cigar1.tmp"
			
			#print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/cigar2plot.py "+ref+" "+self.runname+"/cigar1.tmp allcoverage.plot"
			

			print >> bashfile, "awk '{print $2}' "+self.runname+"/cigar1.tmp > "+self.runname+"/readnames.tmp"
			
			if rtype=='solexa' and multiplexnamed=='y':
				print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/shufflefastqSequences.pl "+self.fastqdir+'_'.join(self.name.split('_')[:2])+"_1_"+self.name.split('_')[-1]+".fastq "+self.fastqdir+'_'.join(self.name.split('_')[:2])+"_2_"+self.name.split('_')[-1]+".fastq "+self.runname+"/shuffled.tmp"
			elif rtype=='solexa':
				print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/shufflefastqSequences.pl "+self.fastqdir+self.name+"_1.fastq "+self.fastqdir+self.name+"_2.fastq "+self.runname+"/shuffled.tmp"
			
			
			print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.6/ssaha_pileup/other_codes/get_seqreads/get_seqreads "+self.runname+"/readnames.tmp "+self.runname+"/shuffled.tmp "+self.runname+"/fastq.tmp"
			
			if velvet=='y' or pairedend=='y':
				print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/exclude/get_excreads "+self.runname+"/readnames.tmp "+self.runname+"/shuffled.tmp "+self.name+"unmap.fastq"
				print >> bashfile, "mv "+self.name+"unmap.fastq "+self.runname+"/unmap.fastq"
			
			if pairedend=='y':
				if rtype=='454':
					print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -seeds 5 -score "+str(self.quality)+" -kmer 13 -skip 4 -diff 0 -output cigar "+ref+" "+self.runname+"/unmap.fastq > "+self.runname+"/unmapcigar1.tmp"
				else:
					print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/ssaha2 -"+self.rtype+" -score "+str(self.quality)+" -skip 2 -diff 0 -output cigar "+ref+" "+self.runname+"/unmap.fastq > "+self.runname+"/unmapcigar1.tmp"
				print >> bashfile, 'grep "^cigar" '+self.runname+"/unmapcigar1.tmp > "+self.runname+"/unmapcigar2.tmp"
				print >> bashfile, "mv "+self.runname+"/unmapcigar2.tmp "+self.runname+"/unmapcigar1.tmp"
			
			
			if rtype=='454':
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 0 -trans 0 "+self.runname+"/cigar1.tmp "+ref+" "+self.runname+"/fastq.tmp > "+self.runname+"/snp.tmp"
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 0 -trans 0 "+self.runname+"/cigar1.tmp "+ref+" "+self.runname+"/fastq.tmp > "+self.runname+"/pileup.tmp"
			elif rtype=='solexa':
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 1 -trans 0 "+self.runname+"/cigar1.tmp "+ref+" "+self.runname+"/fastq.tmp > "+self.runname+"/snp.tmp"
				print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 1 -trans 0 "+self.runname+"/cigar1.tmp "+ref+" "+self.runname+"/fastq.tmp > "+self.runname+"/pileup.tmp"
			print >> bashfile, "egrep ^SNP "+self.runname+"/snp.tmp > "+self.runname+"/all.snp"
			print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_indel -insertion 1 "+self.runname+"/cigar1.tmp "+ref+" "+self.runname+"/pileup.tmp > "+self.runname+"/insertions.txt"
			print >> bashfile, "/nfs/users/nfs_s/sh16/ssaha_2.2/pileup_v0.5/ssaha_pileup/ssaha_pileup/ssaha_indel -deletion 1 "+self.runname+"/cigar1.tmp "+ref+" "+self.runname+"/pileup.tmp > "+self.runname+"/deletions.txt"
			print >> bashfile, "egrep ^cons "+self.runname+"/pileup.tmp > "+self.runname+"/all.pileup"
		
		print >> bashfile, "mv "+self.runname+"/cigar1.tmp "+self.runname+"/all.cigar"	
		#print >> bashfile, "mv "+self.runname+"/cigar2.tmp "+self.runname+"/final.cigar"
		print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/cigar2plot.py "+ref+" "+self.runname+"/all.cigar "+self.runname+"/"+self.name+"_mappedcoverage.plot"
		
		if pairedend=='y':
			print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/cigar2_odd_plot.py "+ref+" "+self.runname+"/abnormal.cigar "+str(self.maxinsertsize)+" "+str(self.mininsertsize)+" 50000 "+self.runname+"/"+self.name
			print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/cigar2_fri_plot.py "+ref+" "+self.runname+"/all.cigar "+str(self.maxinsertsize)+" "+str(self.mininsertsize)+" "+self.runname+"/"+self.name+"_fricoverage.plot"
		else:
			print >> bashfile, "/nfs/users/nfs_s/sh16/scripts/cigar2_fr_plot.py "+ref+" "+self.runname+"/all.cigar "+str(self.maxinsertsize)+" "+str(self.mininsertsize)+" "+self.runname+"/"+self.name+"_frcoverage.plot"
		print >> bashfile, "rm "+self.runname+"/*.tmp"
		print >> bashfile, "gzip -f "+self.runname+"/*"
		if velvet=='y':
			print >> bashfile, "gunzip "+self.runname+"/unmap.fastq.gz"
			
		bashfile.close()
		
		
	def maq_unmap_to_fastq(self):
		
		bashfile=open(self.number+ref+'_velvetscript.sh','w')
		
		print "\nRunning velvet on unassembled reads from "+self.name+'...'
		#Creating single unmap file where there are multiple
		
		print >> bashfile, "test"
		
		if not os.path.isfile(self.runname+"/unmap.fastq"):
			print >> bashfile, "cat "+self.runname+"/unmap*_*.txt > "+self.runname+"/unmap.txt"
		
			#Changing unmap txt file to fastq format
			unmapout=open(self.runname+"/unmap.fastq", "w")
			lines=open(self.runname+"/unmap.txt", "rU").readlines()
			lastread=''
			for line in lines:
				words=line.split()
				if words[0]==lastread and program=='maq':
					print >> unmapout, '@'+words[0].replace('/1','/2')
				else:
					print >> unmapout, '@'+words[0]
				print >> unmapout, words[2]
				print >> unmapout, '+'
				print >> unmapout, words[3]
				lastread=words[0]
				unmapout.flush()
			unmapout.close()
		bashfile.close()
		
		
########
# Main #
########


if __name__ == "__main__":
        argv=sys.argv[1:]
        ref, inputdirs, quality, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mininsertsize, mismatches, program, rtype, readlength, multiplexnamed=getOptions(argv)	
	

	print '\nChecking input files...'
	sys.stdout.flush()
	
	chars = string.ascii_letters + string.digits

	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	pools=[]
	count=0
	poolsort=[]
	for pool in inputdirs:
		
		if not os.path.isfile(pool):
			print "File "+pool+" not found! Skipping..."
			continue
		fastqdir=''
		if pool[-1]=='/':
			pool=pool[:-1]
		if len(pool.split('/'))>1:
			fastqdir='/'.join(pool.split('/')[:-1])+'/'
		filetype='.'+pool.split('.')[-1]
		
		if rtype=='solexa' and filetype!='.fastq':
			print "WARNING: Input file name is not fastq."
		elif rtype=='454' and filetype!='.sff':
			print "WARNING: Input file name is not sff."
		
		pool='.'.join(pool.split('/')[-1].split('.')[:-1])

		
		if multiplexnamed=='y':
			if rtype=='solexa' and pairedend=='y':
				if pool.split('_')[-2]=='1':
					pool='_'.join(pool.split('_')[:2])+'_'+pool.split('_')[-1]
				else:
					continue
		
		else:
			if pairedend=='y':
				if rtype=="solexa" and pool[-2:]=='_1':
					pool=pool[:-2]
					if not os.path.isfile(fastqdir+pool+"_2.fastq"):
						print "File "+pool+"_2.fastq not found! Skipping..."
						continue
				elif rtype=="solexa" and pool[-2:]=='_2':
					pool=pool[:-2]
					if not os.path.isfile(fastqdir+pool+"_1.fastq"):
						print "File "+pool+"_1.fastq not found! Skipping..."
						continue
				elif rtype!='454':
					continue		
		
		

		name=pool
		
		
		if program=='maq':
			pool=pool+'_maq'
		elif program=='ssaha':
			pool=pool+'_ssaha'
		if pool in poolsort:
			continue
		print pool+'...',
		sys.stdout.flush()
		if not os.path.isdir(pool):
		 	print "pool "+pool+" not found! Creating...",
			os.system("mkdir "+pool)
						
		
		pools.append(SNPanalysis())
		pools[count].number=str(count+1)
		pools[count].runname=pool
		pools[count].name=name
		pools[count].fastqdir=fastqdir
		pools[count].filetype=filetype
		#the rest can be globals
		pools[count].quality=quality
		pools[count].rtype=rtype
		pools[count].mindepth=mindepth
		pools[count].maxsnps=maxsnps
		pools[count].mismatches=mismatches
		pools[count].mininsertsize=mininsertsize
		pools[count].maxinsertsize=maxinsertsize
		pools[count].pairedend=pairedend
		pools[count].readlength=readlength
		poolsort.append(pool)
		print 'ok'
		sys.stdout.flush()
		
		count=count+1
		
		
	if len(pools)==0:
		print "\nError: No valid input files!"
		sys.exit()

	#Running Maq, Ssaha and Velvet where required
			
	if program=='ssaha':
		for pool in pools:
			pool.runSsaha(ref, tmpname)
			
		#os.system('echo \'bash ${LSB_JOBINDEX}'+tmpname+'_sbs.sh\' | bsub -R \'select[mem>24000] rusage[mem=24000]\'  -J'+tmpname+'_ssaha"[1-'+str(count)+']%16" -o '+tmpname+'ssaha-%I.out -e '+tmpname+'ssaha-%I.err')# run all ssaha jobs in job array
		os.system('echo \'bash ${LSB_JOBINDEX}'+tmpname+'_sbs.sh\' | bsub -J'+tmpname+'_ssaha"[1-'+str(count)+']%20" -o '+tmpname+'ssaha-%I.out -e '+tmpname+'ssaha-%I.err')# run all ssaha jobs in job array
		
		os.system('bsub -w \'ended('+tmpname+'_ssaha)\' rm *'+tmpname+'_sbs.sh '+tmpname+'*.out '+tmpname+'*.err')#when job array is all done delete the ssaha bash scripts
		
		if velvet=='y':
			for pool in pools:#For each pool run velvet if required, but wait for mapping to finish first
				meaninsert=((pool.maxinsertsize-pool.mininsertsize)/2)+pool.mininsertsize
				if pool.pairedend=='y':
					os.system('bsub -w \'ended('+tmpname+'_ssaha)\' /nfs/users/nfs_s/sh16/scripts/velvet_assembly.sh -p -i '+str(meaninsert)+' -f '+pool.runname+'/unmap.fastq -n')
				else:
					os.system('bsub -w \'ended('+tmpname+'_ssaha)\' /nfs/users/nfs_s/sh16/scripts/velvet_assembly.sh -f '+pool.runname+'/unmap.fastq -n')
		#testing calling my snp caller. Needs options
		#os.system('bsub -w \'ended('+ref+'_ssaha)\' /nfs/users/nfs_s/sh16/scripts/summarise_snps_from_ssahatest.py -o test -a -g -t -r '+ref+' -e '+ref.replace('.dna','.embl')+' -q 60 -d 5 '+' '.join(inputdirs))

	elif program=='maq':
		for pool in pools:
			pool.runMaq(ref, tmpname)
		
		os.system('echo \'bash ${LSB_JOBINDEX}'+tmpname+'_mbs.sh\' | bsub -R \'select[mem>24000] rusage[mem=24000]\' -J'+tmpname+'_maq"[1-'+str(count)+']%20" -o '+tmpname+'maq-%I.out -e '+tmpname+'maq-%I.err')# run all maq jobs in job array
		os.system('bsub -w \'ended('+tmpname+'_maq)\' rm *'+tmpname+'_mbs.sh')#when job array is all done delete the ssaha bash scripts
#		if velvet=='y':
#			for pool in pools:
#				pool.maq_unmap_to_fastq()# not sure how best to do this
#			os.system('echo \'bash ${LSB_JOBINDEX}'+ref+'_velvetscript.sh\' | bsub -w \'ended('+ref+'_ssaha)\' -J'+ref+'_velvet"[1-'+str(count)+']%16" -o '+ref+'velvet-%I.out -e '+ref+'velvet-%I.err')
#			for pool in pools:
#				if pool.pairedend=='y':
#					os.system('bsub -w \'ended('+ref+'_velvet)\' /nfs/users/nfs_s/sh16/scripts/velvet_assembly.sh -p -i '+str(meaninsert)+' -f '+pool.runname+'/unmap.fastq -n')
#				else:
#					os.system('bsub -w \'ended('+ref+'_velvet)\' /nfs/users/nfs_s/sh16/scripts/velvet_assembly.sh -f '+pool.runname+'/unmap.fastq -n')
		#os.system('bsub -w \'ended('+tmpname+'_maq)\' /nfs/users/nfs_s/sh16/scripts/summarise_snps_from_maqtest2.py -o test -a -g -t -r '+tmpname+' -e '+ref.replace('.dna','.embl')+' -q 60 -d 5 '+' '.join(inputdirs))
			
			

			
			
			
			
