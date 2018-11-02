#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from Bio import Seq
from Bio import AlignIO
from optparse import OptionParser, OptionGroup
import shlex, subprocess

sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from Si_SNPs_temp import *
#from  multiprocessing import cpu_count



import time



####################
# Set some globals #
####################




##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-c", "--converter", action="store", dest="converter", help="backconverter file created by Picky_oligo_designer_in.py", default="", metavar="FILE")
	group.add_option("-p", "--pickyfile", action="store", dest="picky", help="Picky output file containing candidate oligos", default="")
	group.add_option("-o", "--output", action="store", dest="output", help="Prefix for output files", default="")
	group.add_option("-d", "--contaminant_database", action="store", dest="contaminants", help="Name file containing contaminant accession numbers", default=False, metavar="FILE")
	group.add_option("-H", "--human", action="store_true", dest="human", help="Blast primers against human genome", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "specificity options")
	group.add_option("-M", "--mismatch", action="store", dest="mismatch", help="Minimum mismatch score (=length-mispriming score calculated as +1 for a matching base, -1 for a mismatching base and -2 for a gap) [default= %default]", default=30, type="int", metavar="FLOAT")
	group.add_option("-3", "--3prime", action="store", dest="threeprime", help="Minimum 3' mismatch score [default= %default]", default=2, type="float", metavar="int")
	group.add_option("-I", "--ignore", action="store", dest="ignore", help="Ignore all hits with mismatch score greater than this value (whether or not they meet the 3' mismatch score) [default= %default]", default=60, type="float", metavar="int")
#	group.add_option("-a", "--self_any", action="store", dest="self_any", help="Maximum self complementarity [default= %default]", default=8.0, type="float", metavar="FLOAT")
#	group.add_option("-e", "--self_end", action="store", dest="self_end", help="Maximum 3' self complementarity [default= %default]", default=3.0, type="float", metavar="FLOAT")
	parser.add_option_group(group)
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.converter=='':
		DoError('No converter file selected')
	elif not os.path.isfile(options.converter):
		DoError('Cannot find file '+options.converter)
	
	if options.picky=='':
		DoError('No picky file selected')
	elif not os.path.isfile(options.picky):
		DoError('Cannot find file '+options.picky)

	
	if options.output=='':
		DoError('No output file prefix selected')
	elif options.output[-1]!="_":
		options.output=options.output+"_"
	return





def print_primers_to_tabfile(handle, primers):
	print >> handle, "ID   primers"
	for primer in primers.keys():
#		if primers[primer]["strand"]=="f":
#			print >> handle, "FT   primer          "+str(primers[primer]["location"][0])+".."+str(primers[primer]["location"][1])
#		elif primers[primer]["strand"]=="r":
#			print >> handle, "FT   primer          complement("+str(primers[primer]["location"][0])+".."+str(primers[primer]["location"][1])+")"
		print >> handle, "FT   primer          "+str(primers[primer]["location"][0])+".."+str(primers[primer]["location"][1])
		print >> handle, "FT                   /name="+primer
		keys=primers[primer].keys()
		keys.sort()
		for qualifier in keys:
			if qualifier!="location":
				print >> handle, "FT                   /"+str(qualifier)+"="+str(primers[primer][qualifier])
		
		#print >> handle, "FT                   /colour="+str(primers[primer]["hit"])



################
# Main program #
################		

if __name__ == "__main__":
	
	starttime=time.clock()

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	#make random name for temporary files
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	
	#Read the backconverter file
	genomelength=0
	converter_dict={}
	foundlength=False
	for line in open(options.converter, "rU"):
		words=line.strip().split()
		if not foundlength:
			if len(words)!=1:
				print "Expecting first line of backconverter file to be genome length. Instead found:"
				print line.strip()
				sys.exit()
			else:
				genomelength=int(words[0])
				foundlength=True
				continue
	
		if len(words)!=3:
			print "Found invalid line in backconverter file:"
			print line.strip()
			sys.exit()
		else:
			if words[0] not in converter_dict:
				try:
					converter_dict[words[0]]=[int(words[1]), int(words[2])]
				except StandardError:
					print "Found invalid line in backconverter file:"
					print line.strip()
					sys.exit()
			else:
				print "Found duplicate entry in backconverter file:"
				print line.strip()
				print words[0], converter_dict[words[0]]
	
	
	pickyin=open(options.picky, "rU").read()
	pickyblocks=pickyin.split("\n\n")
	
	primers={}
	
	for block in pickyblocks:
		lines=block.split("\n")
		try:
			key=lines[1].split()[6]
			primers[key]={}
			primers[key]["sequence"]=lines[0].split()[0]
			primers[key]["location"]=[int(lines[1].split()[4])+converter_dict[key][0],int(lines[1].split()[5])+converter_dict[key][0]]
			primers[key]["length"]=int(lines[0].split()[1])
			primers[key]["melting_temp"]=float(lines[1].split()[1])
		except StandardError:
			if block=="":
				continue
			else:
				print "Found non-standard Picky block:"
				print block
				sys.exit()
		
	
#	print primers
#	sys.exit()
#	
#	
#	position=1
#	primers={}
#	
#	count=0
#	total=0.0
#	hundredth=float(len(consensus_chunks))/100
#	
#	print "Identifying potential primers..."
#	sys.stdout.flush()
#	fcount=0
#	rcount=0
#	
#	
#	
#	for chunk in consensus_chunks:
#	
#		count=count+1
#		if count>=hundredth:
#			total=total+count
#			count=0
#			print "%.0f%% complete\r" % (100*(total/len(consensus_chunks))),
#			sys.stdout.flush()
#	
#		if len(chunk)>options.max_size:
#			
#			prev=0
#			for x in range(0,len(chunk),options.max_size):
#				
#				find potential forward primers
#				
#				handle=open(tmpname+"primer_3_input", "w")
#				
#				if x+(options.max_size*2)<len(chunk):
#					make_primer3_input_file(handle, chunk[x:x+(options.max_size*2)], position+x, opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=options.self_any, self_end=options.self_end)
#					end=x+40
#				elif len(chunk)>x+options.max_size:
#					make_primer3_input_file(handle, chunk[x:], position+x, opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=options.self_any, self_end=options.self_end)
#					end=len(chunk)
#				else:
#					continue
#				handle.close()
#				
#				primerdata=run_primer3(primer_3_input=tmpname+"primer_3_input")		
#
#				keys=primerdata.keys()
#				keys.sort()
#				hit=1
#				for key in keys:
#					primerdata[key]["location"][0]=primerdata[key]["location"][0]+position+x
#					primerdata[key]["location"][1]=primerdata[key]["location"][1]+position+x-1
#					fcount+=1
#					primername="f"+str(primerdata[key]["location"][0])
#					
#					if not primers.has_key(primername):
#						primers[primername]={"hit":hit, "strand":"f"}
#						hit+=1
#						for value in primerdata[key].keys():
#							primers[primername][value]=primerdata[key][value]
#							
#				find potential reverse primers
#				
#				handle=open(tmpname+"primer_3_input", "w")
#				
#				if x+(options.max_size*2)<len(chunk):
#					make_primer3_input_file(handle, revcomp(chunk[x:x+(options.max_size*2)]), position+x, opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=options.self_any, self_end=options.self_end)
#				elif len(chunk)>x+options.max_size:
#					make_primer3_input_file(handle, revcomp(chunk[x:]), position+x, opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=options.self_any, self_end=options.self_end)
#				handle.close()
#				
#				primerdata=run_primer3(primer_3_input=tmpname+"primer_3_input")		
#				
#				keys=primerdata.keys()
#				keys.sort()
#				hit=1
#				
#				for key in keys:
#					primerdata[key]["location"][0]=position+end-primerdata[key]["location"][0]-1
#					primerdata[key]["location"][1]=position+end-primerdata[key]["location"][1]
#					
#					primername="r"+str(primerdata[key]["location"][0])
#					rcount+=1
#					if not primers.has_key(primername):
#						primers[primername]={"hit":hit, "strand":"r"}
#						hit+=1
#						for value in primerdata[key].keys():
#							primers[primername][value]=primerdata[key][value]
#				
#		
#		position+=len(chunk)+1
#		if position>10000:
#			break
#	print "100% complete"
#	sys.stdout.flush()
#
#	
#	print fcount, rcount
	
	
	
	primerlist=primers.keys()
	numseqs=len(primers.keys())
	
	if numseqs==0:
		DoError("No primers identified")
	else:
		print numseqs, "candidate primers identified"
		sys.stdout.flush()
	
	
#		
#	
#	
	
	
	if options.contaminants and os.path.isfile(options.contaminants):
	
		filenum=1
		output=open(tmpname+"primerseqs.fasta"+str(filenum),"w")
		count=0
		
		
		for x, primername in enumerate(primers.keys()):
			count+=1
			if count==100:
				filenum+=1
				output.close()
				output=open(tmpname+"primerseqs.fasta"+str(filenum),"w")
				count=0
			print >> output, ">"+primername
			print >> output, primers[primername]["sequence"]
			primers[primername]["CONTAMINANT_MIN_MISPRIMING_PASS"]=True
			primers[primername]["CONTAMINANT_MIN_MISPRIMING"]=1000
#			if x>90:
#				break
		output.close()
	
		numprimers=x+1
		
		print "Creating database of specified contaminant sequences"
		sys.stdout.flush()
		
		if os.path.isfile(tmpname+"contaminantdb.fasta"):
			os.remove(tmpname+"contaminantdb.fasta")
	
		seqlines=open(options.contaminants, "rU").readlines()
		
		for x, line in enumerate(seqlines):
			accession=line.strip()
			
			print "Downloading sequence for", accession
			sys.stdout.flush()
			os.system("pfetch "+accession+" >> "+tmpname+"contaminantdb.fasta")
#			if x>0:
#				break
		os.system("/software/pubseq/bin/ncbi_blast+/makeblastdb -in "+tmpname+"contaminantdb.fasta -dbtype nucl")
		
		print "Running blast search of "+str(numprimers)+" primers vs "+str(x+1)+" contaminant sequences..."
		sys.stdout.flush()
	
		
		blastallstring='echo \'/software/pubseq/bin/ncbi_blast+/blastn -evalue 100 -db '+tmpname+'contaminantdb.fasta -query '+tmpname+'primerseqs.fasta${LSB_JOBINDEX} -outfmt \"\'10 qseqid sseqid length qseq sseq\'\" -out '+tmpname+'blast.tempout${LSB_JOBINDEX} -word_size 7 -dust no\''
		
		os.system(blastallstring+' | bsub -o '+tmpname+'_contaminantbsubout -e '+tmpname+'_contaminantbsuberr -J'+tmpname+'primerseqs.fasta"[1-'+str(filenum)+']" > '+tmpname+'jobstring')

		
		jobnum=open(tmpname+'jobstring', "rU").read().split(">")[0].split("<")[1]
		todo=numprimers
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
		
		if exited>0:
			DoError("Some of your blasts did not complete")
		
		os.system('cat '+tmpname+'blast.tempout* >>'+tmpname+'blast.out; rm '+tmpname+'blast.tempout* '+tmpname+'jobstring')
		
		#os.system("blastall -p blastn -a "+str(options.processors)+" -d contaminantdb.fasta -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -e 100 -m 3 -W 7 -F F")
		print "\nParsing blast results"
		sys.stdout.flush()
		
		readfile=open(tmpname+"blast.out","rU")
		for line in readfile:
		
			words=line.strip().split(',')
			query=words[0]
			subject=words[1]
			matchlen=int(words[2])
			qalignment=words[3]
			salignment=words[4]
			if query==subject:
				continue
				
			matches=0
			score=0.0
			
			#print matches, cmatches, score, matchpos, started, matchlen
			for matchpos, qhit in enumerate(qalignment):
				shit=salignment[matchpos]
				if qhit==shit:
					matches+=1
					score+=1.0
				elif qhit=="-" or shit=="-":
					score-=2
				elif qhit=="N" or shit=="N":
					score-=0.25
				else:
					score-=1.0
			#print matches, cmatches, score, matchpos, started, matchlen
			if (primers[query]["length"]-score)<options.mismatch:
				primers[query]["CONTAMINANT_MIN_MISPRIMING_PASS"]=False

			if (primers[query]["length"]-score)<primers[query]["CONTAMINANT_MIN_MISPRIMING"]:	
				primers[query]["CONTAMINANT_MIN_MISPRIMING"]=(primers[query]["length"]-score)
				
				
		for primername in primers.keys():
			if primers[primername]["CONTAMINANT_MIN_MISPRIMING"]==1000:
				primers[primername]["CONTAMINANT_MIN_MISPRIMING"]=0
			
			if not primers[primername]["CONTAMINANT_MIN_MISPRIMING_PASS"]:
				del primers[primername]
		primerlist=primers.keys()
		numseqs=len(primers.keys())
		os.system('rm '+tmpname+'blast.out')

		
	
				
	elif options.contaminants:
		print "Cannot find contaminants file", options.contaminants, "skipping..."
	
	
#	handle=open(options.output+"testl_primers.tab", "w")
#	print_primers_to_tabfile(handle, primers)
#	handle.close()
	
	
	if options.human:
	
	
		filenum=1
		output=open(tmpname+"primerseqs.fasta"+str(filenum),"w")
		count=0
		
		
		for x, primername in enumerate(primers.keys()):
			count+=1
			if count==10:
				filenum+=1
				output.close()
				output=open(tmpname+"primerseqs.fasta"+str(filenum),"w")
				count=0
		
			print >> output, ">"+primername
			print >> output, primers[primername]["sequence"]
			primers[primername]["HUMAN_MIN_MISPRIMING_PASS"]=True
			primers[primername]["HUMAN_MIN_MISPRIMING"]=1000
	#		if x>99:
	#			break
		output.close()
	
	
		print "Running blast search of "+str(str(x+1))+" primers vs Human... may be slow!"
		sys.stdout.flush()
		
		
		blastallstring='echo \'/software/pubseq/bin/ncbi_blast+/blastn -evalue 100 -db /lustre/scratch101/blastdb/Supported/hs.fna -query '+tmpname+'primerseqs.fasta${LSB_JOBINDEX} -outfmt \"\'10 qseqid sseqid length qseq sseq\'\" -out '+tmpname+'blast.tempout${LSB_JOBINDEX} -word_size 7 -dust no\''
		
		os.system(blastallstring+' | bsub -o '+tmpname+'_humanbsubout -e '+tmpname+'_humanbsuberr -R \'select[mem>10000] rusage[mem=10000]\'  -J'+tmpname+'primerseqs.fasta"[1-'+str(filenum)+']"  -M 10000000 > '+tmpname+'jobstring')

		
		jobnum=open(tmpname+'jobstring', "rU").read().split(">")[0].split("<")[1]

		todo=x+1
		print "JOBID    ARRAY_SPEC  OWNER  NJOBS PEND DONE  RUN EXIT SSUSP USUSP PSUSP"
		#while os.path.isfile(tmpname+'waitfile'):
		while todo!=0:
			time.sleep(10)
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
		
		if exited>0:
			DoError("Some of your blasts did not complete")
		
		os.system('cat '+tmpname+'blast.tempout* >>'+tmpname+'blast.out; rm '+tmpname+'blast.tempout*')

		
		#os.system("blastall -e 100 -p blastn -a "+str(options.processors)+" -d /lustre/scratch103/sanger/sh16/Chlamydia/primerdesigns/softmasked_dusted.fa -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -m 3 -W 7 -F F")
		print "\nParsing blast results"
		sys.stdout.flush()
		
		readfile=open(tmpname+"blast.out","rU")
		for line in readfile:
		
			words=line.strip().split(',')
			query=words[0]
			subject=words[1]
			matchlen=int(words[2])
			qalignment=words[3]
			salignment=words[4]
			if query==subject:
				continue
				
			matches=0
			score=0.0
			
			#print matches, cmatches, score, matchpos, started, matchlen
			for matchpos, qhit in enumerate(qalignment):
				shit=salignment[matchpos]
				if qhit==shit:
					matches+=1
					score+=1.0
				elif qhit=="-" or shit=="-":
					score-=2
				elif qhit=="N" or shit=="N":
					score-=0.25
				else:
					score-=1.0
			#print matches, cmatches, score, matchpos, started, matchlen
			if (primers[query]["length"]-score)<options.mismatch:
				primers[query]["HUMAN_MIN_MISPRIMING_PASS"]=False

			if (primers[query]["length"]-score)<primers[query]["HUMAN_MIN_MISPRIMING"]:	
				primers[query]["HUMAN_MIN_MISPRIMING"]=(primers[query]["length"]-score)
				
				
		for primername in primers.keys():
			if primers[primername]["HUMAN_MIN_MISPRIMING"]==1000:
				primers[primername]["HUMAN_MIN_MISPRIMING"]=0
			
			if not primers[primername]["HUMAN_MIN_MISPRIMING_PASS"]:
				del primers[primername]
		primerlist=primers.keys()
		numseqs=len(primers.keys())
		os.system('rm '+tmpname+'blast.out')
		
	
	
	
		
	
	
	print len(primers), "primers match specificity requirements\nChecking for primer-dimers..."
	
	
	
	filenum=1
	output=open(tmpname+"primerseqs.fasta"+str(filenum),"w")
	outputb=open(tmpname+"allprimerseqs.fasta","w")
	count=0
	
	
	for x, primername in enumerate(primers):
		count+=1
		if count==100:
			filenum+=1
			output.close()
			output=open(tmpname+"primerseqs.fasta"+str(filenum),"w")
			count=0
	
		print >> output, ">"+primername
		print >> output, primers[primername]["sequence"]
		print >> outputb, ">"+primername
		print >> outputb, primers[primername]["sequence"]
		primers[primername]["PAIRWISE_FAILS"]=set([])
#		if x>99:
#			break
	output.close()
	outputb.close()

	os.system("/software/pubseq/bin/ncbi_blast+/makeblastdb -in "+tmpname+"allprimerseqs.fasta -dbtype nucl")
	print "Running blast search of "+str(str(x+1))+" primers vs themselves"
	sys.stdout.flush()
	
	
	blastallstring='echo \'/software/pubseq/bin/ncbi_blast+/blastn -evalue 100 -db '+tmpname+'allprimerseqs.fasta -query '+tmpname+'primerseqs.fasta${LSB_JOBINDEX} -outfmt \"\'10 qseqid sseqid length qseq sseq\'\" -out '+tmpname+'blast.tempout${LSB_JOBINDEX} -word_size 7 -dust no\''
	
	print blastallstring
	os.system(blastallstring+' | bsub -R \'select[mem>5000] rusage[mem=5000]\'  -J'+tmpname+'primerseqs.fasta"[1-'+str(filenum)+']"  -M 5000000 > '+tmpname+'jobstring')

	
	jobnum=open(tmpname+'jobstring', "rU").read().split(">")[0].split("<")[1]

	todo=x+1
	print "JOBID    ARRAY_SPEC  OWNER  NJOBS PEND DONE  RUN EXIT SSUSP USUSP PSUSP"
	#while os.path.isfile(tmpname+'waitfile'):
	while todo!=0:
		time.sleep(10)
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
	print
	if exited>0:
		DoError("Some of your blasts did not complete")
	
	os.system('cat '+tmpname+'blast.tempout* >>'+tmpname+'blast.out; rm '+tmpname+'blast.tempout*')
	
	#os.system("blastall -e 100 -p blastn -a "+str(options.processors)+" -d /lustre/scratch103/sanger/sh16/Chlamydia/primerdesigns/softmasked_dusted.fa -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -m 3 -W 7 -F F")
	sys.stdout.flush()
	
	readfile=open(tmpname+"blast.out","rU")
	for line in readfile:
	
		words=line.strip().split(',')
		query=words[0]
		subject=words[1]
		matchlen=int(words[2])
		qalignment=words[3]
		salignment=words[4]
		if query==subject:
			continue
			
		matches=0
		score=0.0
		
		#print matches, cmatches, score, matchpos, started, matchlen
		for matchpos, qhit in enumerate(qalignment):
			shit=salignment[matchpos]
			if qhit==shit:
				matches+=1
				score+=1.0
			elif qhit=="-" or shit=="-":
				score-=2
			elif qhit=="N" or shit=="N":
				score-=0.25
			else:
				score-=1.0
		#print matches, cmatches, score, matchpos, started, matchlen
		if (primers[query]["length"]-score)<options.mismatch:
			primers[query]["PAIRWISE_FAILS"].add(subject)
			
			
	os.system('rm '+tmpname+'blast.out')

	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	genomeblocks={}
	blockstarts=[]
	
	for x in range(0,genomelength,5000):
		genomeblocks[x]=0
		blockstarts.append(x)
	
	for primername in primers.keys():
		x=0
		#print x, blockstarts[x],genomelength,  seqFeature.location.nofuzzy_start
		while x<len(blockstarts) and int(primers[primername]["location"][0])>blockstarts[x]:
			#print x, blockstarts[x],genomelength,  seqFeature.location.nofuzzy_start
			x+=1
		genomeblocks[blockstarts[x-1]]+=1
		primers[primername]["block"]=blockstarts[x-1]
		
	
	print "Calculating free energies..."
	count=0.0
	total=0.0
	onepercent=float(len(primers))/100
	for x, primer in enumerate(primers):
		count+=1
		if count>=onepercent:
			total+=count
			print str(int((total/len(primers))*100))+"%\r",
			sys.stdout.flush()
			count=0.0
			break
		handle=open(tmpname+"vrnaeval_in", "w")
		print >> handle, ">tmp"
		print >> handle, primers[primer]["sequence"]
		handle.close()
		os.system("vrnafold -sequence "+tmpname+"vrnaeval_in -outfile "+tmpname+"vrna.fold > /dev/null 2>&1")
		
		primers[primer]["free_energy"]=float(open(tmpname+"vrna.fold","rU").readlines()[1].split()[-1].replace(")","").replace("(",""))
	print "100%"
	
	
	handle=open(options.output+"candidate_primers.tab", "w")
	print_primers_to_tabfile(handle, primers)
	handle.close()
	

	
	discarded=0
	kept=0
	finalprimers={}
	
	for i, primer1 in enumerate(primerlist):
	
		#count=count+1
		#if count>=hundredth:
			#total=total+count
			#count=0
			#print "%.0f%% complete\r" % (100*(total/len(primerlist))),
		print kept, "kept,", discarded, "discarded\r",
		sys.stdout.flush()
	

	
	
	
	
	

		toremove=[]
		keep=True
		for primer2 in primerlist[i+1:]:

			
			if primer2 in primers[primer1]["PAIRWISE_FAILS"]:
			
				
				winner=None

				if genomeblocks[primers[primer1]["block"]]>=genomeblocks[primers[primer2]["block"]]:
					winner=primer2
				elif genomeblocks[primers[primer1]["block"]]<genomeblocks[primers[primer2]["block"]]:
					winner=primer1
	#			elif primers[primer1]["PENALTY"]<primers[primer2]["PENALTY"]:
	#				winner=primer1
	#			elif primers[primer2]["PENALTY"]<primers[primer1]["PENALTY"]:
	#				winner=primer2
				
	#				print "Winner =", winner, primer1score, primer2score, primer1mindiffs, primer2mindiffs, primers[primer1]["PENALTY"], primers[primer2]["PENALTY"]
				
				if winner==primer2:
					keep=False
					genomeblocks[primers[primer1]["block"]]-=1
					break
				elif winner==primer1:
					toremove.append(primer2)
					genomeblocks[primers[primer2]["block"]]-=1
					
		
		
				
		if keep:
			finalprimers[primer1]=primers[primer1]
			kept+=1
			for loser in toremove:
				primerlist.remove(loser)
				discarded+=1
		else:
			for loser in toremove:
				genomeblocks[primers[loser]["block"]]+=1
			discarded+=1
	
#	os.system("blastall -p blastn -S 1 -e 100 -a 4 -d "+tmpname+"primerseqs.fasta -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -m 3 -W 7 -F F")
#	print "Parsing blast results"
#	sys.stdout.flush()
	
	
	
	
	
	print kept, "kept,", discarded, "discarded"
	print len(finalprimers), "primers found that do not form primer-dimers"
	
	
#	rcount=0
#	fcount=0
#	for f in finalprimers.keys():
#		if finalprimers[f]["strand"]=="f":
#			fcount+=1
#		elif finalprimers[f]["strand"]=="r":
#			rcount+=1
#	
#	print fcount, rcount
	
	

	
	handle=open(options.output+"best_primers.tab", "w")
	print_primers_to_tabfile(handle, finalprimers)
	handle.close()
	
	print "Cleaning up..."
	os.system("rm "+tmpname+"*")
	print "Results can be found in "+options.output+"all_primers.tab and "+options.output+"best_primers.tab"
	
	print time.clock()-starttime
	sys.exit()
	
	
