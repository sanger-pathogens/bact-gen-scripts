#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math
from scipy.stats import chi2



def Usage():
	print 'moving_trees.py Usage:'
	print 'moving_trees.py [options] <input alignment file>'
	print 'Options:'
	print '-w\t\twindow size'
	print '-s\t\tstep size'
	print '-r\t\treference tree'
	print '-f\n\nforce (overwrites previous moving_tree files)'


#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hw:s:fr:", ["window=", "step=", "force=", "ref="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	window=1000
	step=500
	force='n'
	ref=''
	
	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-w", "--window"):
			window=int(arg)
		elif opt in ("-s", "--step"):
			step=int(arg)
		elif opt in ("-r", "--ref"):
			ref=arg
		elif opt in ("-f", "--force"):
			force='y'
	
	inputfile=args[0]
	
	
	if inputfile==[]:
		print 'Error: No input file selected!'
		Usage()
		sys.exit()
	if ref!='' and not os.path.isfile(ref):
		print 'Error: Cannot find reference file!'
		Usage()
		sys.exit()

	return window, step, inputfile, force, ref
	


	
if __name__ == "__main__":
	argv=sys.argv[1:]
	window, step, inputfile, force, ref=getOptions(argv)
	
	if force=='n' and os.path.isdir("moving_trees_dir"):
		userinput=''
		print "\nmoving_trees_dir directory already exists. If you continue some contents may be replaced."
		while userinput not in ['y','n', 'Y', 'N', 'Q']:
			userinput=raw_input("Continue? (y/n/Q): ")
		userinput=userinput.lower()
		if userinput in ['n', 'q']:
			sys.exit()
		else:
			os.system("rm -f moving_trees_dir/RAxML*")
	elif not os.path.isdir("moving_trees_dir"):
		os.mkdir("./moving_trees_dir/")
	
	print '\nReading input alignment...',
	sys.stdout.flush()
	
	sequences={}
	#for x in hassnps.keys():
		#line=lines[x-1]
	
	
	lines=open(inputfile, "rU").read().split('>')[1:]
	
	sequences={}
	
	for line in lines:
		words=line.strip().split('\n')
		sequences[words[0]]=''.join(words[1:])
	
		
	seqlen=len(sequences[sequences.keys()[0]])
	
	if ref!='':
		os.system("cp "+ref+" moving_trees_dir/whole.tre")
	
	os.chdir("moving_trees_dir/")
	
	for sequence in sequences.keys():
		if len(sequences[sequence])!=seqlen:
			print "\nERROR!: sequences are not all of the same length!!!\n"
			sys.exit()
		sequences[sequence]=sequences[sequence].replace('N','-')
			
	
	print "Found", len(sequences.keys()), "sequences of length", seqlen
	sys.stdout.flush()	
	
	print "\nIdentifying variable sites...",
	sys.stdout.flush()	
	
	snplocations=[]
	snpbases=[]
	basecounts=[]
	snpsequences={}
	for key in sequences.keys():
		snpsequences[key]=''
	
	for x in range(seqlen):
		numbases=0
		foundbases={}	
		for key in sequences.keys():
			base=sequences[key][x].upper()
			if base not in foundbases.keys() and base!='-':
				foundbases[base]=1
				numbases=numbases+1
			elif base!='-':
				foundbases[base]=foundbases[base]+1
		
		if numbases>1:
			snplocations.append(x)
			snpbases.append(foundbases)
			basecounts.append(numbases)
			for key in sequences.keys():
				snpsequences[key]=snpsequences[key]+sequences[key][x]
				
	
	snpseqlen=len(snpsequences[snpsequences.keys()[0]])
	
	for sequence in snpsequences.keys():
		if len(snpsequences[sequence])!=snpseqlen:
			print "\nERROR!: snp sequences are not all of the same length!!!\n"
			sys.exit()
	
	print "Found", snpseqlen, "variable sites"
	sys.stdout.flush()	
	
	
	if ref=='':
		
		allaln=open("whole.aln", 'w')
		print >> allaln, len(snpsequences.keys()), snpseqlen
		for sequence in sequences.keys():
			print >> allaln, sequence, snpsequences[sequence]
		allaln.close()
		
		print "\nRunning tree of whole alignment...",
		sys.stdout.flush()	
		os.system("RAxML -f d -m GTRGAMMA -s whole.aln -n whole > RAxML.tmp")
		os.system("mv RAxML_result.whole whole.tre")
		os.system("rm -f whole.aln RAxML*")
		print "Done"
		
	
	count=0
	snprfs=[0.0]*snpseqlen
	snpgds=[0.0]*snpseqlen
	numcalcs=[0]*snpseqlen
	for base in range(0,snpseqlen,int(step)):
		alignfile=open("window.aln", "w")
		start=base-(window/2)
		end=base+(window/2)
		print >> alignfile, len(snpsequences.keys()), (base+(window/2))-(base-(window/2))
		for sequence in snpsequences.keys():
		
			if start>=0 and end<snpseqlen:	
				print >> alignfile, sequence, snpsequences[sequence][start:end]
			elif start<0 and end<snpseqlen:
				print >> alignfile, sequence, snpsequences[sequence][snpseqlen+start:]+snpsequences[sequence][:end]
			elif start>=0 and end>=snpseqlen:
				print >> alignfile, sequence, snpsequences[sequence][start:]+snpsequences[sequence][:end-snpseqlen]
			else:
				print "Error! Your window size is larger than your sequence!"
				sys.exit()
			
		alignfile.close()
		count=count+1
		os.system("RAxML -f d -m GTRGAMMA -s window.aln -n "+str(count)+" > RAxML.tmp")
		if os.path.isfile("RAxML_result."+str(count)):
			

			treeout=open("trees.tre","w")
			treein=open("RAxML_result."+str(count), 'rU').readlines()
			print >> treeout, treein[0].strip().replace(':0.0;',';')

			treein=open("whole.tre", 'rU').readlines()
                        print >> treeout, treein[0].strip().replace(':0.0;',';')
                        print >> treeout, treein[0].strip().replace(':0.0;',';')

			treeout.close()
			sys.exit()
			os.system('~sh16/hashRF/hashrf trees.tre 2 > rfout.tmp')
			rfout=open('rfout.tmp', 'rU').readlines()
			rf=float(rfout[11].split()[1])
			os.system('java -jar ~sh16/geodemaps_v0.2/geodeMAPS.jar -v trees.tre > gdout.tmp')
			#os.system("cat gdout.tmp")
			gdout=open('output.txt', 'rU').readlines()
			gd=float(gdout[0].split()[2])
		else:
			rf=-1
			gd=-1
		print rf, gd
		
		if start>=0 and end<snpseqlen and rf>-1 and gd>-1:
			for x in range(start,end):
				snprfs[x]=snprfs[x]+rf
				snpgds[x]=snpgds[x]+gd
				numcalcs[x]=numcalcs[x]+1

		elif start<0:
			for x in range(0,end):
                                snprfs[x]=snprfs[x]+rf
                                snpgds[x]=snpgds[x]+gd
                                numcalcs[x]=numcalcs[x]+1
			for x in range(snpseqlen+start,snpseqlen):
                                snprfs[x]=snprfs[x]+rf
                                snpgds[x]=snpgds[x]+gd
                                numcalcs[x]=numcalcs[x]+1

		elif end>snpseqlen:
                        for x in range(0,end-snpseqlen):
                                snprfs[x]=snprfs[x]+rf
                                snpgds[x]=snpgds[x]+gd
                                numcalcs[x]=numcalcs[x]+1
                        for x in range(start,snpseqlen):
                                snprfs[x]=snprfs[x]+rf
                                snpgds[x]=snpgds[x]+gd
                                numcalcs[x]=numcalcs[x]+1


		os.system("rm -f trees.tre rfout.tmp gdout.tmp RAxML*")
	

	gdoutput=open("gd_w"+str(window)+"_w"+str(step)+".plot",'w')
	rfoutput=open("rf_w"+str(window)+"_w"+str(step)+".plot",'w')
	
	print >> rfoutput, "# BASE Robinson_Foulds_Distance"
	print >> gdoutput, "# BASE Geodseic_Distance"

	for x, snp in enumerate(snplocations):
		if numcalcs[x]>0:
			print >> gdoutput, snp, float(snpgds[x])/numcalcs[x]
			print >> rfoutput, snp, float(snprfs[x])/numcalcs[x]
		else:
			print >> gdoutput, snp, "0"
                        print >> rfoutput, snp, "0"

	os.system("rm -f window.aln*")
	os.chdir("../")
	
	
	
	
	
	
	
	
	
	
	
	
	

