#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from Bio import Seq
from Bio import AlignIO
from optparse import OptionParser, OptionGroup
import tre

sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from Si_SNPs_temp import *
#from  multiprocessing import cpu_count



import time



####################
# Set some globals #
####################


PRIMER_3_LOCATION="/nfs/users/nfs_s/sh16/primer3-2.2.2-beta/src/"






##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-i", "--input", action="store", dest="alignment", help="Input alignment file name", default="", metavar="FILE")
	group.add_option("-p", "--prefix", action="store", dest="prefix", help="Prefix for output files", default="")
	group.add_option("-d", "--contaminant_database", action="store", dest="contaminants", help="Name file containing contaminant accession numbers", default=False, metavar="FILE")
	group.add_option("-H", "--human", action="store_true", dest="human", help="Blast primers against human genome", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "primer3 options")
	group.add_option("-o", "--opt_size", action="store", dest="opt_size", help="Optimum primer size [default= %default]", default=17, type="int", metavar="INT")
	group.add_option("-s", "--min_size", action="store", dest="min_size", help="Minimum primer size [default= %default]", default=15, type="int", metavar="INT")
	group.add_option("-S", "--max_size", action="store", dest="max_size", help="Maximum primer size [default= %default]", default=24, type="int", metavar="INT")
	
	group.add_option("-m", "--opt_tm", action="store", dest="opt_tm", help="Optimum melting temperature (in Celcius) [default= %default]", default=42.5, type="float", metavar="FLOAT")
	group.add_option("-t", "--min_tm", action="store", dest="min_tm", help="Minimum melting temperature (in Celcius) [default= %default]", default=40, type="float", metavar="FLOAT")
	group.add_option("-T", "--max_tm", action="store", dest="max_tm", help="Maximum melting temperature (in Celcius) [default= %default]", default=45, type="float", metavar="FLOAT")

	group.add_option("-c", "--opt_gc", action="store", dest="opt_gc", help="Optimum GC content [default= %default]", default=50, type="float", metavar="FLOAT")
	group.add_option("-g", "--min_gc", action="store", dest="min_gc", help="Minimum GC content [default= %default]", default=20, type="float", metavar="FLOAT")
	group.add_option("-G", "--max_gc", action="store", dest="max_gc", help="Maximum GC content [default= %default]", default=80, type="float", metavar="FLOAT")
	
	group.add_option("-l", "--gc_clamp", action="store", dest="gc_clamp", help="Length of 3' GC clamp [default= %default]", default=1, type="int", metavar="INT")
	group.add_option("-x", "--poly_x", action="store", dest="poly_x", help="Maximum poly-x length allowed in primer [default= %default]", default=5, type="int", metavar="INT")
	
	group.add_option("-a", "--self_any", action="store", dest="self_any", help="Maximum self complementarity [default= %default]", default=8.0, type="float", metavar="FLOAT")
	group.add_option("-e", "--self_end", action="store", dest="self_end", help="Maximum 3' self complementarity [default= %default]", default=3.0, type="float", metavar="FLOAT")
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "specificity options")
	group.add_option("-M", "--mismatch", action="store", dest="mismatch", help="Minimum mismatch score (=length-mispriming score calculated as +1 for a matching base, -1 for a mismatching base and -2 for a gap) [default= %default]", default=4, type="int", metavar="FLOAT")
	group.add_option("-3", "--3prime", action="store", dest="threeprime", help="Minimum 3' mismatch score [default= %default]", default=2, type="float", metavar="int")
	group.add_option("-I", "--ignore", action="store", dest="ignore", help="Ignore all hits with mismatch score greater than this value (whether or not they meet the 3' mismatch score) [default= %default]", default=6, type="float", metavar="int")
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
	if options.alignment=='':
		DoError('No alignment file selected')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment)
	
	elif options.opt_size<1 or options.opt_size>36:
		DoError('Primer optimum size must be between 1 and 36 (This upper limit is governed by maximum oligo size for which melting-temperature calculations are valid)')
	elif options.opt_size>options.max_size:
		DoError('Primer optimum size must be the same or smaller than maximum size')
	elif options.opt_size<options.min_size:
		DoError('Primer optimum size must be the same or larger than minimum size')
		
	
	elif options.opt_tm<1:
		DoError('Primer optimum melting temperature must be greater than 1')
	elif options.opt_tm>options.max_tm:
		DoError('Primer optimum melting temperature must be the same or lower than maximum melting temperature')
	elif options.opt_tm<options.min_tm:
		DoError('Primer optimum melting temperature must be the same or higher than minimum melting temperature')
	
	
	elif options.opt_gc<0 or options.opt_gc>100:
		DoError('Primer optimum gc must be between 0 and 100')
	elif options.opt_gc>options.max_gc:
		DoError('Primer optimum gc must be the same or smaller than maximum size')
	elif options.opt_gc<options.min_gc:
		DoError('Primer optimum gc must be the same or larger than minimum size')
	
	
	elif options.gc_clamp<0 or options.gc_clamp>options.min_size:
		DoError('GC-clamp length must be between 0 and '+str(options.min_size))
	
	if options.prefix!="":
		options.prefix=options.prefix+"_"
	return



def make_primer3_input_file(handle, sequence1, starting_pos, sequence2="", opt_size=17, min_size=15, max_size=24, opt_tm=42.5, min_tm=40, max_tm=45, min_gc=20, opt_gc=50, max_gc=80, gc_clamp=1, poly_x=5, self_any=6.00, self_end=4.00):
	
	
	
	print >> handle, "SEQUENCE_ID="+str(starting_pos)
	
	if sequence2=="":
		print >> handle, "PRIMER_TASK=pick_sequencing_primers"
		print >> handle, "SEQUENCE_TEMPLATE="+sequence1
		print >> handle, "PRIMER_PICK_LEFT_PRIMER=1"
		print >> handle, "PRIMER_PICK_INTERNAL_OLIGO=0"
		print >> handle, "PRIMER_PICK_RIGHT_PRIMER=0"
		print >> handle, "PRIMER_LEFT_NUM_RETURN=1"
	else:
		print >> handle, "PRIMER_TASK=check_primers"
		print >> handle, "SEQUENCE_PRIMER="+sequence1
		print >> handle, "SEQUENCE_PRIMER_REVCOMP="+sequence2
		print >> handle, "PRIMER_PICK_LEFT_PRIMER=1"
		print >> handle, "PRIMER_PICK_INTERNAL_OLIGO=0"
		print >> handle, "PRIMER_PICK_RIGHT_PRIMER=1"

	print >> handle, "PRIMER_OPT_SIZE="+str(opt_size)
	print >> handle, "PRIMER_MIN_SIZE="+str(min_size)
	print >> handle, "PRIMER_MAX_SIZE="+str(max_size)
	print >> handle, "PRIMER_OPT_TM="+str(opt_tm)
	print >> handle, "PRIMER_MIN_TM="+str(min_tm)
	print >> handle, "PRIMER_MAX_TM="+str(max_tm)
	print >> handle, "PRIMER_MAX_NS_ACCEPTED=0"
	print >> handle, "PRIMER_MIN_GC="+str(min_gc)
	print >> handle, "PRIMER_GC_CLAMP="+str(gc_clamp)
	print >> handle, "PRIMER_OPT_GC_PERCENT="+str(opt_gc)
	print >> handle, "PRIMER_MAX_GC="+str(max_gc)
	print >> handle, "PRIMER_MAX_POLY_X="+str(poly_x)
	print >> handle, "PRIMER_SELF_ANY="+str(self_any)
	print >> handle, "PRIMER_SELF_END="+str(self_end)
	print >> handle, "PRIMER_TM_FORMULA=1"
	print >> handle, "PRIMER_SALT_CORRECTIONS=1"
	#print >> handle, "PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION="+str(len(sequence))+",20 16,20"
	print >> handle, "="




def print_primers_to_tabfile(handle, primers):
	print >> handle, "ID   primers"
	for primer in primers.keys():
		if primers[primer]["strand"]=="f":
			print >> handle, "FT   primer          "+str(primers[primer]["location"][0])+".."+str(primers[primer]["location"][1])
		elif primers[primer]["strand"]=="r":
			print >> handle, "FT   primer          complement("+str(primers[primer]["location"][0])+".."+str(primers[primer]["location"][1])+")"
		print >> handle, "FT                   /name="+primer
		print >> handle, "FT                   /length="+str(len(primers[primer]["SEQUENCE"]))
		keys=primers[primer].keys()
		keys.sort()
		for qualifier in keys:
			if qualifier!="location":
				print >> handle, "FT                   /"+str(qualifier)+"="+str(primers[primer][qualifier])
		
		print >> handle, "FT                   /colour="+str(primers[primer]["hit"])



def run_primer3(primer_3_input="primer_3_input"):

	primer3out=os.popen(PRIMER_3_LOCATION+"primer3_core < "+primer_3_input)
	primerdata={}
	
	for line in primer3out:
		
		if line[:11]=="PRIMER_LEFT":
		
			variable=line.split("=")[0]
			variablebits=variable.split("_")
			value=line.split("=")[1].strip()
			try:
				hit=int(variablebits[2])+1
				prefix="_".join(variablebits[:2])
				prefixlen=3
			except StandardError:
				hit=1
				prefix=""
				prefixlen=2
			
			if "_".join(variablebits[prefixlen:])=="NUM_RETURNED":
				if int(value)==0:
					return {}
				continue
			
			if not primerdata.has_key(hit):
				primerdata[hit]={}
			
			if "_".join(variablebits[prefixlen:])=="":
				primerdata[hit]["location"]=[int(value.split(',')[0]), int(value.split(',')[1])+int(value.split(',')[0])]
			else:
				primerdata[hit]["_".join(variablebits[prefixlen:])]=value
	return primerdata
	
def run_primer3_dimer_check(primer_3_input, sequence1, sequence2):

	primer3out=os.popen(PRIMER_3_LOCATION+"primer3_core < "+primer_3_input)
	primerdata={}
	
	COMP_ANY=""
	COMP_END=""
	phit=-1
	lhit=-2
	rhit=-3
	lseq=""
	rseq=""
	
	for line in primer3out:
		
		if line[:11]=="PRIMER_LEFT":
			
			lvariable=line.split("=")[0]
			lvariablebits=lvariable.split("_")
			lvalue=line.split("=")[1].strip()
			try:
				lhit=int(lvariablebits[2])+1
				lprefix="_".join(lvariablebits[:2])
				lprefixlen=3
			except StandardError:
				lhit=1
				lprefix=""
				lprefixlen=2
			
			if "_".join(lvariablebits[lprefixlen:])=="SEQUENCE":
				lseq=lvalue
		
		elif line[:12]=="PRIMER_RIGHT":
		
			rvariable=line.split("=")[0]
			rvariablebits=rvariable.split("_")
			rvalue=line.split("=")[1].strip()
			try:
				rhit=int(rvariablebits[2])+1
				rprefix="_".join(rvariablebits[:2])
				rprefixlen=3
			except StandardError:
				rhit=1
				rprefix=""
				rprefixlen=2

			
			if "_".join(rvariablebits[rprefixlen:])=="SEQUENCE":
				rseq=rvalue
		
		elif line[:11]=="PRIMER_PAIR":
		
			pvariable=line.split("=")[0]
			pvariablebits=pvariable.split("_")
			pvalue=line.split("=")[1].strip()
			try:
				phit=int(pvariablebits[2])+1
				pprefix="_".join(pvariablebits[:2])
				pprefixlen=3
			except StandardError:
				phit=1
				pprefix=""
				pprefixlen=2

			
			if "_".join(pvariablebits[pprefixlen:])=="COMPL_ANY":
				COMP_ANY=pvalue
			elif "_".join(pvariablebits[pprefixlen:])=="COMPL_END":
				COMP_END=pvalue					
			
		if phit==rhit and phit==lhit and sequence1==lseq and sequence2==rseq and COMP_ANY!="" and COMP_END!="":
			
			return float(COMP_ANY), float(COMP_END)
				
	
	if COMP_ANY=="":
		COMP_ANY=float(len(sequence1))
	if COMP_END=="":
		COMP_END=float(len(sequence1))
	return float(COMP_ANY), float(COMP_END)

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
	
	
	#Read the alignment file

	try:
		alignment=read_alignment(options.alignment)
	except StandardError:
		DoError("Cannot open alignment file")

	consensus=consensus_from_alignment(alignment, unknowns=["N"]).replace("-", "N")
	
	consensus_chunks=consensus.split("N")
	
	
	pickyout=open("Picky_input.seq","w")
	conversion=open("Picky_backconverter.txt","w")
	count=0
	position=1
	for chunk in consensus_chunks:
		if len(chunk)>options.max_size:
			for x in range(0,len(chunk),options.max_size):
				if x+(options.max_size*2)<len(chunk):
					print >> pickyout, ">Picky"+str(count)
					print >> pickyout, chunk[x:x+(options.max_size*2)]
					print >> conversion, str(count), str(position+x), str(position+x+(options.max_size*2))
					count+=1
				elif len(chunk)>x+options.max_size:
					print >> pickyout, ">Picky"+str(count)
					print >> pickyout, chunk[x:]
					print >> conversion, str(count), str(position+x), str(position+len(chunk))
					count+=1
				else:
					continue
				
				
			count+=1
		position+=len(chunk)+1
	pickyout.close()
	sys.exit()
	
	position=1
	primers={}
	
	count=0
	total=0.0
	hundredth=float(len(consensus_chunks))/100
	
	print "Identifying potential primers..."
	sys.stdout.flush()
	fcount=0
	rcount=0
	
	
	
	for chunk in consensus_chunks:
	
		count=count+1
		if count>=hundredth:
			total=total+count
			count=0
			print "%.0f%% complete\r" % (100*(total/len(consensus_chunks))),
			sys.stdout.flush()
	
		if len(chunk)>options.max_size:
			
			prev=0
			for x in range(0,len(chunk),options.max_size):
				
				#find potential forward primers
				
				handle=open(tmpname+"primer_3_input", "w")
				
				if x+(options.max_size*2)<len(chunk):
					make_primer3_input_file(handle, chunk[x:x+(options.max_size*2)], position+x, opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=options.self_any, self_end=options.self_end)
					end=x+40
				elif len(chunk)>x+options.max_size:
					make_primer3_input_file(handle, chunk[x:], position+x, opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=options.self_any, self_end=options.self_end)
					end=len(chunk)
				else:
					continue
				handle.close()
				
				primerdata=run_primer3(primer_3_input=tmpname+"primer_3_input")		

				keys=primerdata.keys()
				keys.sort()
				hit=1
				for key in keys:
					primerdata[key]["location"][0]=primerdata[key]["location"][0]+position+x
					primerdata[key]["location"][1]=primerdata[key]["location"][1]+position+x-1
					fcount+=1
					primername="f"+str(primerdata[key]["location"][0])
					
					if not primers.has_key(primername):
						primers[primername]={"hit":hit, "strand":"f"}
						hit+=1
						for value in primerdata[key].keys():
							primers[primername][value]=primerdata[key][value]
							
				#find potential reverse primers
				
				handle=open(tmpname+"primer_3_input", "w")
				
				if x+(options.max_size*2)<len(chunk):
					make_primer3_input_file(handle, revcomp(chunk[x:x+(options.max_size*2)]), position+x, opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=options.self_any, self_end=options.self_end)
				elif len(chunk)>x+options.max_size:
					make_primer3_input_file(handle, revcomp(chunk[x:]), position+x, opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=options.self_any, self_end=options.self_end)
				handle.close()
				
				primerdata=run_primer3(primer_3_input=tmpname+"primer_3_input")		
				
				keys=primerdata.keys()
				keys.sort()
				hit=1
				
				for key in keys:
					primerdata[key]["location"][0]=position+end-primerdata[key]["location"][0]-1
					primerdata[key]["location"][1]=position+end-primerdata[key]["location"][1]
					
					primername="r"+str(primerdata[key]["location"][0])
					rcount+=1
					if not primers.has_key(primername):
						primers[primername]={"hit":hit, "strand":"r"}
						hit+=1
						for value in primerdata[key].keys():
							primers[primername][value]=primerdata[key][value]
				
		
		position+=len(chunk)+1
#		if position>10000:
#			break
	print "100% complete"
	sys.stdout.flush()

	
	print fcount, rcount
	
	
	
	primerlist=primers.keys()
	numseqs=len(primers.keys())
	
	if numseqs==0:
		DoError("No primers found by primer3")
	else:
		print numseqs, "candidate primers identified by primer3"
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
			print >> output, primers[primername]["SEQUENCE"]
			primers[primername]["CONTAMINANT_MIN_MISPRIMING_PASS"]=True
			primers[primername]["CONTAMINANT_MIN_MISPRIMING"]=0
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
		os.system("formatdb -i "+tmpname+"contaminantdb.fasta -p F")

		
		print "Running blast search of "+str(numprimers)+" primers vs "+str(x+1)+" contaminant sequences..."
		sys.stdout.flush()
	
		
		blastallstring='echo \'blastall -e 100 -p blastn -d '+tmpname+'contaminantdb.fasta -i '+tmpname+'primerseqs.fasta${LSB_JOBINDEX} -o '+tmpname+'blast.tempout${LSB_JOBINDEX} -m 3 -W 7 -F F\''
		
		os.system(blastallstring+' | bsub  -J'+tmpname+'primerseqs.fasta"[1-'+str(filenum)+']" > '+tmpname+'jobstring')

		
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
			
		
		os.system('cat '+tmpname+'blast.tempout* >>'+tmpname+'blast.out; rm '+tmpname+'blast.tempout* '+tmpname+'jobstring')
		
		#os.system("blastall -p blastn -a "+str(options.processors)+" -d contaminantdb.fasta -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -e 100 -m 3 -W 7 -F F")
		print "\nParsing blast results"
		sys.stdout.flush()
		
		currseq=""
		readfile=open(tmpname+"blast.out","rU")
		for line in readfile:
			line.strip()
			stop=False
			if len(line.split())==2 and line.split()[0]=="Query=":
				currseq=line.split()[1]
				currseqlen=len(primers[currseq]["SEQUENCE"])
				primers[currseq]["CONTAMINANT_MIN_MISPRIMING"]=1000
				
				while not ( len(line.split())==4 and len(line.split()[0])>2 and line.split()[0][-2:]=="_0") and not (len(line.split())>0 and line.split()[0]=="BLASTN"):
					try:
						line=readfile.next().strip()
					except StopIteration:
						stop=True
						break
						
				if stop or (len(line.split())>0 and line.split()[0]=="BLASTN"):
					continue
				
				start=string.find(line, line.split()[2])
				end=start+len(line.split()[2])
				missingcbit=currseqlen-int(line.strip().split()[3])
				bestmatch=4
				line=readfile.next().strip()
				
				
				while len(line.split())==4 and line.split()[0]!="BLASTN":
					matches=0
					cmatches=0
					score=0.0
					matchlen=len(line[end:start:-1].strip())
					matchpos=0
					started=False
					#print matches, cmatches, score, matchpos, started, matchlen
					for x, hit in enumerate(line[end:start:-1]):
						if started:
							matchpos+=1
						elif hit!=" ":
							started=True
						if hit==".":
							matches+=1
							if x+missingcbit<5:
								cmatches+=1
							score+=1.0
						elif started and hit==" " and matchpos<=matchlen:
							score-=2
						elif started and hit=="N" and matchpos<=matchlen:
							score-=0.25
						elif started and matchpos<=matchlen:
							score-=1.0
					#print matches, cmatches, score, matchpos, started, matchlen
					if ((5-cmatches)<options.threeprime and (currseqlen-score)<options.ignore) or (currseqlen-score)<options.mismatch:
						primers[currseq]["CONTAMINANT_MIN_MISPRIMING_PASS"]=False
	
					if (currseqlen-score)<primers[currseq]["CONTAMINANT_MIN_MISPRIMING"]:	
						primers[currseq]["CONTAMINANT_MIN_MISPRIMING"]=(currseqlen-score)
						
					line=readfile.next().strip()
				
		
		for primername in primers.keys():
			if primers[primername]["CONTAMINANT_MIN_MISPRIMING"]==1000:
				primers[primername]["CONTAMINANT_MIN_MISPRIMING"]=0
			#print len(primers[currseq]["SEQUENCE"]), primers[primername]["CONTAMINANT_MAX_MISPRIMING"], primers[primername]["CONTAMINANT_2_2_PASS"], primers[primername]["CONTAMINANT_4_PASS"]
			if not primers[primername]["CONTAMINANT_MIN_MISPRIMING_PASS"]:#primers[primername]["CONTAMINANT_MAX_MISPRIMING"]>12:#(len(primers[primername]["SEQUENCE"])-3) and not primers[primername]["CONTAMINANT_2_2_PASS"] and not primers[primername]["CONTAMINANT_4_PASS"]:
				del primers[primername]
		primerlist=primers.keys()
		numseqs=len(primers.keys())
		
		rcount=0
		fcount=0
		for f in primers.keys():
			if primers[f]["strand"]=="f":
				fcount+=1
			elif primers[f]["strand"]=="r":
				rcount+=1
		
		print fcount, rcount
		#sys.exit()
		os.system('rm '+tmpname+'blast.out')
			
	elif options.contaminants:
		print "Cannot find contaminants file", options.contaminants, "skipping..."
	
	
#	handle=open(options.prefix+"testl_primers.tab", "w")
#	print_primers_to_tabfile(handle, primers)
#	handle.close()
	
	
	if options.human:
	
	
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
			print >> output, primers[primername]["SEQUENCE"]
			primers[primername]["HUMAN_MIN_MISPRIMING_PASS"]=True
			primers[primername]["HUMAN_MIN_MISPRIMING"]=0
	#		if x>99:
	#			break
		output.close()
	
	
		print "Running blast search of "+str(str(x+1))+" primers vs Human... may be slow!"
		sys.stdout.flush()
		
		
		blastallstring='echo \'blastall -e 100 -p blastn -d /lustre/scratch101/sanger/sh16/Chlamydia/primerdesigns/ncbi36_unmasked -i '+tmpname+'primerseqs.fasta${LSB_JOBINDEX} -o '+tmpname+'blast.tempout${LSB_JOBINDEX} -m 3 -W 7 -F F\''
		
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
		
		os.system('cat '+tmpname+'blast.tempout* >>'+tmpname+'blast.out; rm '+tmpname+'blast.tempout*')

		
		#os.system("blastall -e 100 -p blastn -a "+str(options.processors)+" -d /lustre/scratch103/sanger/sh16/Chlamydia/primerdesigns/softmasked_dusted.fa -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -m 3 -W 7 -F F")
		print "\nParsing blast results"
		sys.stdout.flush()
		stop=False
		currseq=""
		readfile=open(tmpname+"blast.out","rU")
		for line in readfile:
			stop=False
			line.strip()
			if len(line.split())==2 and line.split()[0]=="Query=":
				currseq=line.split()[1]
				currseqlen=len(primers[currseq]["SEQUENCE"])
				primers[currseq]["HUMAN_MIN_MISPRIMING"]=1000
				
				
				while not ( len(line.split())==4 and len(line.split()[0])>2 and line.split()[0][-2:]=="_0") and not (len(line.split())>0 and line.split()[0]=="BLASTN"):
					try:
						line=readfile.next().strip()
					except StopIteration:
						stop=True
						break
						
				if stop or (len(line.split())>0 and line.split()[0]=="BLASTN"):
					continue
				
				start=string.find(line, line.split()[2])
				end=start+len(line.split()[2])
				bestmatch=4
				missingcbit=currseqlen-int(line.strip().split()[3])
				line=readfile.next().strip()
				
				
				while len(line.split())==4 and line.split()[0]!="BLASTN":
					matches=0
					cmatches=0
					score=0.0
					matchlen=len(line[end:start:-1].strip())
					matchpos=0
					started=False
					#print matches, cmatches, score, matchpos, started, matchlen
					for x, hit in enumerate(line[end:start:-1]):
						if started:
							matchpos+=1
						elif hit!=" ":
							started=True
						if hit==".":
							matches+=1
							if x+missingcbit<5:
								cmatches+=1
							score+=1.0
						elif started and hit==" " and matchpos<=matchlen:
							score-=2
						elif started and hit=="N" and matchpos<=matchlen:
							score-=0.25
						elif started and matchpos<=matchlen:
							score-=1.0
					#print matches, cmatches, score, matchpos, started, matchlen
					if ((5-cmatches)<options.threeprime and (currseqlen-score)<options.ignore) or (currseqlen-score)<options.mismatch:
						primers[currseq]["HUMAN_MIN_MISPRIMING_PASS"]=False
	
					if (currseqlen-score)<primers[currseq]["HUMAN_MIN_MISPRIMING"]:	
						primers[currseq]["HUMAN_MIN_MISPRIMING"]=(currseqlen-score)
						
					line=readfile.next().strip()
				
		
		for primername in primers.keys():
			if primers[primername]["HUMAN_MIN_MISPRIMING"]==1000:
				primers[primername]["HUMAN_MIN_MISPRIMING"]=0
			
			#print len(primers[currseq]["SEQUENCE"]), primers[primername]["CONTAMINANT_MAX_MISPRIMING"], primers[primername]["CONTAMINANT_2_2_PASS"], primers[primername]["CONTAMINANT_4_PASS"]
			if not primers[primername]["HUMAN_MIN_MISPRIMING_PASS"]:#primers[primername]["CONTAMINANT_MAX_MISPRIMING"]>12:#(len(primers[primername]["SEQUENCE"])-3) and not primers[primername]["CONTAMINANT_2_2_PASS"] and not primers[primername]["CONTAMINANT_4_PASS"]:
				del primers[primername]
		primerlist=primers.keys()
		numseqs=len(primers.keys())
		rcount=0
		fcount=0
		for f in primers.keys():
			if primers[f]["strand"]=="f":
				fcount+=1
			elif primers[f]["strand"]=="r":
				rcount+=1
		
		print fcount, rcount
	#os.system("formatdb -i "+tmpname+"primerseqs.fasta -p F")
	
	os.system('mv '+tmpname+'blast.out test.blastout')
	handle=open(options.prefix+"candidate_primers.tab", "w")
	print_primers_to_tabfile(handle, primers)
	handle.close()	
	
	
	
#	print "Identifying primer-dimers..."
#	sys.stdout.flush()
#
	#count=0
	#total=0.0
	#hundredth=float(len(primerlist))/100
	
#	
##	checklist=[]
##	mindiffslist=[]
#	maxmispriminglist=[]
#	if options.human and options.contaminants:
##		checklist=checklist+[["HUMAN_2_4_PASS",1], ["CONTAMINANT_2_4_PASS",2], ["HUMAN_2_2_PASS",1], ["CONTAMINANT_2_2_PASS",1], ["HUMAN_4_PASS",1], ["CONTAMINANT_4_PASS",1]]
##		mindiffslist=mindiffslist+["CONTAMINANT_MIN_DIFFS", "HUMAN_MIN_DIFFS"]
#		maxmispriminglist=maxmispriminglist+["CONTAMINANT_MAX_MISPRIMING", "HUMAN_MAX_MISPRIMING"]
#	elif options.human:
##		checklist=checklist+[["HUMAN_2_4_PASS",1], ["HUMAN_2_2_PASS",1], ["HUMAN_4_PASS",1]]
##		mindiffslist=mindiffslist+["HUMAN_MIN_DIFFS"]
#		maxmispriminglist=maxmispriminglist+["HUMAN_MAX_MISPRIMING"]
#	elif options.contaminants:
##		checklist=checklist+[["CONTAMINANT_2_4_PASS",1], ["CONTAMINANT_2_2_PASS",1], ["CONTAMINANT_4_PASS",1]]
##		mindiffslist=mindiffslist+["CONTAMINANT_MIN_DIFFS"]
#		maxmispriminglist=maxmispriminglist+["CONTAMINANT_MAX_MISPRIMING"]
	
	genomelength=alignment.get_alignment_length()
	
	genomeblocks={"f":{}, "r":{}}
	blockstarts=[]
	
	for x in range(0,genomelength,5000):
		genomeblocks["f"][x]=0
		genomeblocks["r"][x]=0
		blockstarts.append(x)
	
	for primername in primers.keys():
		x=0
		#print x, blockstarts[x],genomelength,  seqFeature.location.nofuzzy_start
		while x<len(blockstarts) and int(primers[primername]["location"][0])>blockstarts[x]:
			#print x, blockstarts[x],genomelength,  seqFeature.location.nofuzzy_start
			x+=1
		genomeblocks[primers[primername]["strand"]][blockstarts[x-1]]+=1
		primers[primername]["block"]=blockstarts[x-1]
		
	
	
	print len(primers), "primers match specificity requirements\nChecking for primer-dimers..."

	
	discarded=0
	kept=0
	finalprimers={}
	
	for x, primer1 in enumerate(primerlist):
	
		#count=count+1
		#if count>=hundredth:
			#total=total+count
			#count=0
			#print "%.0f%% complete\r" % (100*(total/len(primerlist))),
		print kept, "kept,", discarded, "discarded\r",
		sys.stdout.flush()
	
#		primer1score=0
#		for test in checklist:
#			if primers[primer1][test[0]]:
#				primer1score+=test[1]
#		primer1mindiffs=len(primers[primer1]["SEQUENCE"])
#		for test in mindiffslist:
#			if int(primers[primer1][test])<primer1mindiffs:
#				primer1mindiffs=int(primers[primer1][test])
#		primer1maxmispriming=0
#		for test in maxmispriminglist:
#			if primers[primer1][test]>primer1maxmispriming:
#				primer1maxmispriming=primers[primer1][test]
		toremove=[]
		keep=True
		for primer2 in primerlist[x+1:]:
			handle=open(tmpname+"primer_3_input", "w")
				
			make_primer3_input_file(handle, primers[primer1]["SEQUENCE"], position+x, sequence2=primers[primer2]["SEQUENCE"], opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=25, self_end=25)
			handle.close()
			COMP_ANY, COMP_END=run_primer3_dimer_check(tmpname+"primer_3_input", primers[primer1]["SEQUENCE"], primers[primer2]["SEQUENCE"])
			#print COMP_ANY, COMP_END
			
			if len(primers[primer1]["SEQUENCE"])<=primers[primer2]["SEQUENCE"]:	
				minlen=len(primers[primer1]["SEQUENCE"])
			else:
				minlen=len(primers[primer2]["SEQUENCE"])
			
			#if (minlen-float(COMP_ANY))<options.mismatch or float(COMP_END)>5:
			if float(COMP_ANY)>options.self_any or float(COMP_END)>options.self_end:
#				print primer1, primer2, "form dimer"
				winner=None

#				primer2score=0
#				for test in checklist:
#					if primers[primer2][test[0]]:
#						primer2score+=test[1]
				
				
#				primer2mindiffs=len(primers[primer2]["SEQUENCE"])
#				for test in mindiffslist:
#					if int(primers[primer2][test])<primer2mindiffs:
#						primer2mindiffs=int(primers[primer2][test])
				
#				primer2maxmispriming=0
#				for test in maxmispriminglist:
#					if primers[primer2][test]>primer2maxmispriming:
#						primer2maxmispriming=primers[primer2][test]
				
#				if primer1score>primer2score:	
#					winner=primer1
#				elif primer2score>primer1score:	
#					winner=primer2
#				elif primer2mindiffs<primer1mindiffs:
#					winner=primer1
#				elif primer1mindiffs<primer2mindiffs:
#					winner=primer2
				if genomeblocks[primers[primer1]["strand"]][primers[primer1]["block"]]>genomeblocks[primers[primer2]["strand"]][primers[primer2]["block"]]:
					winner=primer2
				elif genomeblocks[primers[primer1]["strand"]][primers[primer1]["block"]]<genomeblocks[primers[primer2]["strand"]][primers[primer2]["block"]]:
					winner=primer1
				elif primers[primer1]["PENALTY"]<primers[primer2]["PENALTY"]:
					winner=primer1
				elif primers[primer2]["PENALTY"]<primers[primer1]["PENALTY"]:
					winner=primer2
				
#				print "Winner =", winner, primer1score, primer2score, primer1mindiffs, primer2mindiffs, primers[primer1]["PENALTY"], primers[primer2]["PENALTY"]
				
				if winner==primer2:
					keep=False
					genomeblocks[primers[primer1]["strand"]][primers[primer1]["block"]]-=1
					break
				elif winner==primer1:
					toremove.append(primer2)
					genomeblocks[primers[primer2]["strand"]][primers[primer2]["block"]]-=1
					
		
		
				
		if keep:
			finalprimers[primer1]=primers[primer1]
			kept+=1
			for loser in toremove:
				primerlist.remove(loser)
				discarded+=1
		else:
			for loser in toremove:
				genomeblocks[primers[loser]["strand"]][primers[loser]["block"]]+=1
			discarded+=1
	
#	os.system("blastall -p blastn -S 1 -e 100 -a 4 -d "+tmpname+"primerseqs.fasta -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -m 3 -W 7 -F F")
#	print "Parsing blast results"
#	sys.stdout.flush()
	
	
	
	
	
	print kept, "kept,", discarded, "discarded"
	print len(finalprimers), "primers found that do not form primer-dimers"
	
	
	rcount=0
	fcount=0
	for f in finalprimers.keys():
		if finalprimers[f]["strand"]=="f":
			fcount+=1
		elif finalprimers[f]["strand"]=="r":
			rcount+=1
	
	print fcount, rcount
	
	

	
	handle=open(options.prefix+"best_primers.tab", "w")
	print_primers_to_tabfile(handle, finalprimers)
	handle.close()
	
	print "Cleaning up..."
	os.system("rm "+tmpname+"*")
	print "Results can be found in "+options.prefix+"all_primers.tab and "+options.prefix+"best_primers.tab"
	
	print time.clock()-starttime
	sys.exit()
	
	
