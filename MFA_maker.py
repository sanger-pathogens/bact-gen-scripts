#!/usr/bin/env python
# !/usr/bin/env python
import string
import os, sys
from Bio import Seq
from Bio import AlignIO
from optparse import OptionParser
import tre

sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
from Si_SNPs_temp import *
from  multiprocessing import cpu_count



import time


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-a", "--alignment", action="store", dest="alignment", help="alignment file name", default="", metavar="FILE")
#	parser.add_option("-t", "--tree", action="store", dest="tree", help="tree file", default="", metavar="FILE")
#	parser.add_option("-o", "--outgroup", action="store", dest="outgroup", help="outgroup", default="")
	parser.add_option("-p", "--prefix", action="store", dest="prefix", help="prefix for output files", default="")
#	parser.add_option("-T", "--transformation", action="store", dest="transformation", help="transformation type (acctran, deltran or ML). [Default= %default]", default="acctran", type="choice", choices=["acctran","deltran", "ML"])
#	parser.add_option("-R", "--RAxML", action="store_true", dest="runtree", help="run phylogeny with RAxML [default=%default]", default=False)
#	parser.add_option("-m", "--model", action="store", dest="model", help="Model of evolution to use. [Default= %default]", default="GTRGAMMA", type="choice", choices=["GTRGAMMA","GTRGAMMAI", "GTRCAT", "GTRMIX", "GTRMIXI"])
#	parser.add_option("-b", "--bootstrap", action="store", dest="bootstrap", help="Number of bootstrap replicates (0 = do not run bootstrap). [Default= %default]", default=100, type="int", metavar="int")
#	parser.add_option("-e", "--embl", action="store", dest="embl", help="Embl/genbank annotation file for reference strain (for dN/dS etc.)", default="", metavar="FILE")
	parser.add_option("-c", "--contaminants", action="store", dest="contaminants", help="Name file containing contaminant accession numbers", default=False, metavar="FILE")
	parser.add_option("-A", "--processors", action="store", dest="processors", help="Number of processors to use for blast [default=%default]", default=1, type="int")
	parser.add_option("-H", "--human", action="store_true", dest="human", help="Blast primers against human genome", default=False)
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.embl!="" and options.reference=='':
#		DoError('No reference selected! If you give an embl file you must specify which sequence it is linked to (i.e. the reference)')
	if options.alignment=='':
		DoError('No alignment file selected')
	elif not os.path.isfile(options.alignment):
		DoError('Cannot find file '+options.alignment)
	if options.prefix!="":
		options.prefix=options.prefix+"_"
	if options.processors>cpu_count():
		DoError('You do not have '+str(options.processors)+' processors')
	elif options.processors<1:
		DoError('Number of processors must be > 0')
		
	return



def make_primer3_input_file(handle, sequence1, starting_pos, sequence2="", opt_size=17, min_size=15, max_size=24, opt_tm=42.5, min_tm=40, max_tm=45, min_product_size=1000, max_product_size=10000, min_gc=20, opt_gc=50, max_gc=80, gc_clamp=1, poly_x=5, self_any=8.00, self_end=3.00):
	
	
	
	print >> handle, "PRIMER_SEQUENCE_ID="+str(starting_pos)
	
	if sequence2=="":
		print >> handle, "SEQUENCE="+sequence1+"N"*1010+"agacgagtagcgttctctgagaatcaaaattctttctttgatggcttcccaacag"
	else:
		print >> handle, "SEQUENCE="+sequence1+"N"*1010+revcomp(sequence2)
		print >> handle, "PRIMER_LEFT_INPUT="+sequence1
		print >> handle, "PRIMER_RIGHT_INPUT="+sequence2

	print >> handle, "PRIMER_OPT_SIZE="+str(opt_size)
	print >> handle, "PRIMER_MIN_SIZE="+str(min_size)
	print >> handle, "PRIMER_MAX_SIZE="+str(max_size)
	print >> handle, "PRIMER_OPT_TM="+str(opt_tm)
	print >> handle, "PRIMER_MIN_TM="+str(min_tm)
	print >> handle, "PRIMER_MAX_TM="+str(max_tm)
	print >> handle, "PRIMER_NUM_NS_ACCEPTED=0"
	print >> handle, "PRIMER_PRODUCT_SIZE_RANGE="+str(min_product_size)+"-"+str(max_product_size)
	print >> handle, "PRIMER_MIN_GC="+str(min_gc)
	print >> handle, "PRIMER_GC_CLAMP="+str(gc_clamp)
	print >> handle, "PRIMER_OPT_GC_PERCENT="+str(opt_gc)
	print >> handle, "PRIMER_MAX_GC="+str(max_gc)
	print >> handle, "PRIMER_NUMBER_TO_RETURN=1"
	print >> handle, "PRIMER_MAX_POLY_X="+str(poly_x)
	print >> handle, "PRIMER_SELF_ANY="+str(self_any)
	print >> handle, "PRIMER_SELF_END="+str(self_end)
	#print >> handle, "PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION="+str(len(sequence))+",20 16,20"
	print >> handle, "="




def print_primers_to_tabfile(handle, primers):
	print >> handle, "ID   primers"
	for primer in primers.keys():
		if primers[primer]["strand"]=="f":
			print >> handle, "FT   primer          "+str(primers[primer]["location"][0])+".."+str(primers[primer]["location"][1])
		elif primers[primer]["strand"]=="r":
			print >> handle, "FT   primer          complement("+str(primers[primer]["location"][0])+".."+str(primers[primer]["location"][1])+")"
		print >> handle, "FT                   /length="+str(len(primers[primer]["SEQUENCE"]))
		for qualifier in primers[primer].keys():
			if qualifier!="location":
				print >> handle, "FT                   /"+str(qualifier)+"="+str(primers[primer][qualifier])
		
		print >> handle, "FT                   /colour="+str(primers[primer]["hit"])



def run_primer3(primer_3_input="primer_3_input"):

	primer3out=os.popen("primer3_core < "+primer_3_input)
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

			if not primerdata.has_key(hit):
				primerdata[hit]={}
			
			if "_".join(variablebits[prefixlen:])=="":
				primerdata[hit]["location"]=[int(value.split(',')[0]), int(value.split(',')[1])+int(value.split(',')[0])]
			else:
				primerdata[hit]["_".join(variablebits[prefixlen:])]=value
	return primerdata
	
def run_primer3_dimer_check(primer_3_input, sequence1, sequence2):

	primer3out=os.popen("primer3_core < "+primer_3_input)
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
			
		if phit==rhit and phit==lhit and sequence1==lseq and sequence2==rhit and COMP_ANY!="" and COMP_END!="":
			return COMP_ANY, COMP_END
				
	
	if COMP_ANY=="":
		COMP_ANY="0.0"
	if COMP_END=="":
		COMP_END="0.0"
	return COMP_ANY, COMP_END

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
	
	position=1
	primers={}
	
	count=0
	total=0.0
	hundredth=float(len(consensus_chunks))/100
	
	print "Identifying potential primers..."
	sys.stdout.flush()
	
	for chunk in consensus_chunks:
	
		count=count+1
		if count>=hundredth:
			total=total+count
			count=0
			print "%.0f%% complete\r" % (100*(total/len(consensus_chunks))),
			sys.stdout.flush()
	
		if len(chunk)>30:
			
			prev=0
			for x in range(0,len(chunk),250):
				
				#find potential forward primers
				
				handle=open(tmpname+"primer_3_input", "w")
				
				if x+500<len(chunk):
					make_primer3_input_file(handle, chunk[x:x+500], position+x)
					end=x+500
				elif len(chunk)>x+30:
					make_primer3_input_file(handle, chunk[x:], position+x)
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
					
					primername="f"+str(primerdata[key]["location"][0])
					
					if not primers.has_key(primername):
						primers[primername]={"hit":hit, "strand":"f"}
						hit+=1
						for value in primerdata[key].keys():
							primers[primername][value]=primerdata[key][value]
							
				#find potential reverse primers
				
				handle=open(tmpname+"primer_3_input", "w")
				
				if x+500<len(chunk):
					make_primer3_input_file(handle, revcomp(chunk[x:x+500]), position+x)
				elif len(chunk)>x+30:
					make_primer3_input_file(handle, revcomp(chunk[x:]), position+x)
				handle.close()
				
				primerdata=run_primer3(primer_3_input=tmpname+"primer_3_input")		

				keys=primerdata.keys()
				keys.sort()
				hit=1
				
				for key in keys:
					primerdata[key]["location"][0]=position+end-primerdata[key]["location"][0]-1
					primerdata[key]["location"][1]=position+end-primerdata[key]["location"][1]
					
					primername="r"+str(primerdata[key]["location"][0])
					
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

	
#	patterns={}
#	
#	for primername in primers.keys():
#		if patterns.has_key(primername):
#			print "Error: Two of your regions have the same name:", primername
#			print "Skipping..."
#			continue
#		patterns[primername]=[]
#		patterns[primername].append(tre.compile(primers[primername]["SEQUENCE"], tre.EXTENDED))
#		patterns[primername].append(tre.compile(revcomp(primers[primername]["SEQUENCE"]), tre.EXTENDED))
#	
#	
#	seqlines=open("contaminantdb.txt", "rU").readlines()
#	
#	
#	fz = tre.Fuzzyness(maxerr = 3)
#	for line in seqlines:
#		accession=line.strip()
#		
#		print "Searching primers against", accession
#		#seqobject=read_alignment(os.popen("pfetch "+accession).read())
#		seq=os.popen("pfetch "+accession).read()[1:].replace("\n","")
#		count=0
#		total=0.0
#		hundredth=float(len(patterns.keys()))/100
#		
#		
#		for primername in patterns.keys():
#			
#			count=count+1
#			if count>=hundredth:
#				total=total+count
#				count=0
#				print "%.0f%% complete\r" % (100*(total/len(patterns.keys()))),
#				sys.stdout.flush()
#			#search for first primer forward
#			matches = patterns[primername][0].match(seq, fz)
#			mincost1=4
#			primer1matches=[]
#			if matches:
#				mincost1=matches.cost
#				for match in matches.groups():
#					primer1matches.append([match,"f"])
#			
#			#search for first primer reerse
#			matches = patterns[primername][1].match(seq, fz)
#			if matches and matches.cost<=mincost1:
#				if matches.cost<mincost1:
#					primer1matches=[]
#					mincost1=matches.cost
#				for match in matches.groups():
#					primer1matches.append([match,"r"])
#			
#			
#			
#			if len(primer1matches)>0:
#				primers[primername][accession+"_matches"]=str(len(primer1matches))
#				primers[primername][accession+"_match_cost"]=str(mincost1)
#	
#	print "100% complete"
#	sys.stdout.flush()
	
	
	
	primerlist=primers.keys()
	numseqs=len(primers.keys())
	
	if numseqs==0:
		DoError("No primers found by primer3")
	
#		
#	
#	
	
	
	if options.contaminants and os.path.isfile(options.contaminants):
	
	
		output=open(tmpname+"primerseqs.fasta","w")
		for x, primername in enumerate(primers.keys()):
			print >> output, ">"+primername
			print >> output, primers[primername]["SEQUENCE"]
			primers[primername]["CONTAMINANT_2_2_PASS"]=True
			primers[primername]["CONTAMINANT_2_4_PASS"]=True
			primers[primername]["CONTAMINANT_4_PASS"]=True
			primers[primername]["CONTAMINANT_MIN_DIFFS"]=len(primers[primername]["SEQUENCE"])
		#		if x>99:
		#			break
		output.close()
	
		print "Creating database of specified contaminant sequences"
		sys.stdout.flush()
		
		if os.path.isfile(tmpname+"contaminantdb.fasta"):
			os.remove(tmpname+"contaminantdb.fasta")
	
		seqlines=open("contaminantdb.txt", "rU").readlines()
		
		for x, line in enumerate(seqlines):
			accession=line.strip()
			
			print "Downloading sequence for", accession
			sys.stdout.flush()
			os.system("pfetch "+accession+" >> "+tmpname+"contaminantdb.fasta")
#			if x>0:
#				break
		os.system("formatdb -i "+tmpname+"contaminantdb.fasta -p F")
		
		
		
		
		print "Running blast search vs specified contaminant sequences..."
		sys.stdout.flush()
		os.system("blastall -p blastn -a "+str(options.processors)+" -d contaminantdb.fasta -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -e 100 -m 3 -W 7 -F F")
		print "Parsing blast results"
		sys.stdout.flush()
		
		currseq=""
		readfile=open(tmpname+"blast.out","rU")
		for line in readfile:
			
			if len(line.split())==4 and len(line.split()[0])>2 and line.split()[0][-2:]=="_0":
				
				currseq=primerlist[int(line.split()[0].split("_")[0])-1]
				currseqlen=len(primers[currseq]["SEQUENCE"])
				start=string.find(line, line.split()[2])
				end=start+currseqlen
				bestmatch=4
				line=readfile.next()
				
				while len(line.split())==4 and line.split()[0]!="BLASTN":
					matches=0
					cmatches=0
					for x, hit in enumerate(line[end:start:-1]):
						if hit==".":
							matches+=1
							if x<5:
								cmatches+=1
					
					if (5-cmatches)<2:
						primers[currseq]["CONTAMINANT_2_2_PASS"]=False
						primers[currseq]["CONTAMINANT_2_4_PASS"]=False
					if ((end-start)-matches)<4:
						primers[currseq]["CONTAMINANT_4_PASS"]=False
						primers[currseq]["CONTAMINANT_2_4_PASS"]=False
	
					if ((end-start)-matches)<primers[currseq]["CONTAMINANT_MIN_DIFFS"]:
						primers[currseq]["CONTAMINANT_MIN_DIFFS"]=((end-start)-matches)
						
					line=readfile.next()
		
		for primername in primers.keys():
			if primers[primername]["CONTAMINANT_MIN_DIFFS"]==0:
				del primers[primername]
		primerlist=primers.keys()
		numseqs=len(primers.keys())
			
	elif options.contaminants:
		print "Cannot find contaminants file", options.contaminants, "skipping..."
	
	
	
	
	
	if options.human:
	
	
		output=open(tmpname+"primerseqs.fasta","w")
		for x, primername in enumerate(primers.keys()):
			print >> output, ">"+primername
			print >> output, primers[primername]["SEQUENCE"]
			primers[primername]["HUMAN_2_2_PASS"]=True
			primers[primername]["HUMAN_2_4_PASS"]=True
			primers[primername]["HUMAN_4_PASS"]=True
			primers[primername]["HUMAN_MIN_DIFFS"]=len(primers[primername]["SEQUENCE"])
	#		if x>99:
	#			break
		output.close()
	
	
		print "Running blast search vs Human... may be slow!"
		sys.stdout.flush()
		os.system("blastall -e 100 -p blastn -a "+str(options.processors)+" -d /lustre/scratch103/sanger/sh16/Chlamydia/primerdesigns/softmasked_dusted.fa -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -m 3 -W 7 -F F")
		print "Parsing blast results"
		sys.stdout.flush()
		
		currseq=""
		readfile=open(tmpname+"blast.out","rU")
		for line in readfile:
			
			if len(line.split())==4 and len(line.split()[0])>2 and line.split()[0][-2:]=="_0":
				#print line.strip()
				currseq=primerlist[int(line.split()[0].split("_")[0])-1]
				currseqlen=len(primers[currseq]["SEQUENCE"])
				start=string.find(line, line.split()[2])
				end=start+currseqlen
				bestmatch=4
				line=readfile.next()
				
				while len(line.split())==4 and line.split()[0]!="BLASTN":
					matches=0
					cmatches=0
					for x, hit in enumerate(line[end:start:-1]):
						if hit==".":
							matches+=1
							if x<5:
								cmatches+=1
					
					if (5-cmatches)<2:
						#print "cmatch fail", (5-cmatches)
						primers[currseq]["HUMAN_2_2_PASS"]=False
						primers[currseq]["HUMAN_2_4_PASS"]=False
					if ((end-start)-matches)<4:
						primers[currseq]["HUMAN_4_PASS"]=False
						primers[currseq]["HUMAN_2_4_PASS"]=False
	
					if ((end-start)-matches)<primers[currseq]["HUMAN_MIN_DIFFS"]:
						primers[currseq]["HUMAN_MIN_DIFFS"]=((end-start)-matches)
						
					line=readfile.next()
	
		for primername in primers.keys():
			if primers[primername]["CONTAMINANT_MIN_DIFFS"]==0:
				del primers[primername]
		primerlist=primers.keys()
		numseqs=len(primers.keys())
	#os.system("formatdb -i "+tmpname+"primerseqs.fasta -p F")
	
	
	
	
	
	print "Identifying primer-dimers..."
	sys.stdout.flush()

	count=0
	total=0.0
	hundredth=float(len(primerlist))/100
	
	
	checklist=[]
	mindiffslist=[]
	if options.human and options.contaminants:
		checklist=checklist+[["HUMAN_2_4_PASS",1], ["CONTAMINANT_2_4_PASS",2], ["HUMAN_2_2_PASS",1], ["CONTAMINANT_2_2_PASS",1], ["HUMAN_4_PASS",1], ["CONTAMINANT_4_PASS",1]]
		mindiffslist=mindiffslist+["CONTAMINANT_MIN_DIFFS", "HUMAN_MIN_DIFFS"]
	elif options.human:
		checklist=checklist+[["HUMAN_2_4_PASS",1], ["HUMAN_2_2_PASS",1], ["HUMAN_4_PASS",1]]
		mindiffslist=mindiffslist+["HUMAN_MIN_DIFFS"]
	elif options.contaminants:
		checklist=checklist+[["CONTAMINANT_2_4_PASS",1], ["CONTAMINANT_2_2_PASS",1], ["CONTAMINANT_4_PASS",1]]
		mindiffslist=mindiffslist+["CONTAMINANT_MIN_DIFFS"]
	
	finalprimers={}
	
	for x, primer1 in enumerate(primerlist):
	
		count=count+1
		if count>=hundredth:
			total=total+count
			count=0
			print "%.0f%% complete\r" % (100*(total/len(primerlist))),
			sys.stdout.flush()
	
		primer1score=0
		for test in checklist:
			if primers[primer1][test[0]]:
				primer1score+=test[1]
		primer1mindiffs=len(primers[primer1]["SEQUENCE"])
		for test in mindiffslist:
			if int(primers[primer1][test])<primer1mindiffs:
				primer1mindiffs=int(primers[primer1][test])
		keep=True
		for primer2 in primerlist[x+1:]:
			handle=open(tmpname+"primer_3_input", "w")
				
			make_primer3_input_file(handle, primers[primer1]["SEQUENCE"], position+x, sequence2=primers[primer2]["SEQUENCE"], self_any=25.00, self_end=25.00)
			handle.close()
			COMP_ANY, COMP_END=run_primer3_dimer_check(tmpname+"primer_3_input", primers[primer1]["SEQUENCE"], primers[primer2]["SEQUENCE"])
			
			if float(COMP_ANY)>8 or float(COMP_END)>3:
#				print primer1, primer2, "form dimer"
				winner=None

				primer2score=0
				for test in checklist:
					if primers[primer2][test[0]]:
						primer2score+=test[1]
				
				
				primer2mindiffs=len(primers[primer2]["SEQUENCE"])
				for test in mindiffslist:
					if int(primers[primer2][test])<primer2mindiffs:
						primer2mindiffs=int(primers[primer2][test])
				
				if primer1score>primer2score:	
					winner=primer1
				elif primer2score>primer1score:	
					winner=primer2
				elif primer2mindiffs<primer1mindiffs:
					winner=primer1
				elif primer1mindiffs<primer2mindiffs:
					winner=primer2
				elif primers[primer1]["PENALTY"]<primers[primer2]["PENALTY"]:
					winner=primer1
				elif primers[primer2]["PENALTY"]<primers[primer1]["PENALTY"]:
					winner=primer2
				
#				print "Winner =", winner, primer1score, primer2score, primer1mindiffs, primer2mindiffs, primers[primer1]["PENALTY"], primers[primer2]["PENALTY"]
				
				if winner==primer2:
					keep=False
					break
				elif winner==primer1:
					primerlist.remove(primer2)
					
					
		if keep:
			finalprimers[primer1]=primers[primer1]
	
#	os.system("blastall -p blastn -S 1 -e 100 -a 4 -d "+tmpname+"primerseqs.fasta -i "+tmpname+"primerseqs.fasta -o "+tmpname+"blast.out -m 3 -W 7 -F F")
#	print "Parsing blast results"
#	sys.stdout.flush()
	
	
	
	
	
	
	
	
	
	
	
	
	handle=open(options.prefix+"all_primers.tab", "w")
	print_primers_to_tabfile(handle, primers)
	handle.close()
	
	handle=open(options.prefix+"best_primers.tab", "w")
	print_primers_to_tabfile(handle, finalprimers)
	handle.close()
	
#	print "Cleaning up..."
#	os.system("rm "+tmpname+"*")
	print "Results can be found in "+options.prefix+"all_primers.tab and "+options.prefix+"best_primers.tab"
	
	print time.clock()-starttime
	sys.exit()
	
	
	
#	
#	print "Running fasta search vs specified contaminant sequences..."
#	os.system("fasta34_t -q -T 4 -E 1 -H -m 9i primerseqs.fasta contaminantdb.fasta > blast.out")
#	print "Parsing fasta results"
#	
#	currseq=""
#	readfile=open("blast.out","rU")
#	for line in readfile:
#		
#		if len(line)>3 and line[:3]==">>>":
#			#print line.strip()
#			currseq=line.split()[0][3:-1]
#			currseqlen=len(primers[currseq]["SEQUENCE"])
#			bestmatch=4
#		elif  len(line)>2 and line[:2]==">>":
#			currmatch=line.split()[0][2:]
#		elif len(line.split())>0 and line.split()[0]=="Smith-Waterman":
#			percent=float(line.split()[3][:-1])
#			overlap=int(line.split()[8])
#			nonmatches=currseqlen-int((percent/100)*overlap)
#			if nonmatches==bestmatch and nonmatches<4:
#				primers[primername][accession+"_matches"]+=1
#			if nonmatches<bestmatch:
#				primers[primername][accession+"_matches"]=1
#				primers[primername][accession+"_match_cost"]=nonmatches
#				bestmatch=nonmatches
#		elif len(line.split())>0 and line.split()[0]==currseq:
#				end=len(line.strip())
#				bits=line.strip()[len(currseq):].split(" ")
#				start=len(bits)-1+len(currseq)
#				#print line[start:end]
#				#print line
#				line=readfile.next()
#				#print line[start:end],
#				matches=0
#				cmatches=0
#				for x, hit in enumerate(line[end:start:-1]):
#					if hit==":":
#						matches+=1
#						if x<5:
#							cmatches+=1
#				
#				if (5-cmatches)<2:
#					#print "cmatch fail", (5-cmatches)
#					primers[currseq]["CONTAMINANT_2_2_PASS"]=False
#					primers[currseq]["CONTAMINANT_2_4_PASS"]=False
#				if ((end-start)-matches)<4:
#					primers[currseq]["CONTAMINANT_4_PASS"]=False
#					primers[currseq]["CONTAMINANT_2_4_PASS"]=False
#
#				if ((end-start)-matches)<primers[currseq]["CONTAMINANT_MIN_DIFFS"]:
#					primers[currseq]["CONTAMINANT_MIN_DIFFS"]=((end-start)-matches)
#					#print "match fail", ((end-start)-matches)
#				
#				#print matches, cmatches, , 5-cmatches
##				for x in print line[start:end]:
#	
#	
#	output=open("primerseqs.fasta","w")
#	numseqs=0
#	for x, primername in enumerate(primers.keys()):
#		if primers[primername]["CONTAMINANT_2_4_PASS"] or primers[primername]["CONTAMINANT_2_2_PASS"]:
#			numseqs+=1
#			print >> output, ">"+primername
#			print >> output, primers[primername]["SEQUENCE"]
#			primers[primername]["HUMAN_2_2_PASS"]=True
#			primers[primername]["HUMAN_2_4_PASS"]=True
#			primers[primername]["HUMAN_4_PASS"]=True
#			primers[primername]["HUMAN_MIN_DIFFS"]=len(primers[primername]["SEQUENCE"])
#		if x>99:
#			break
#	output.close()
#	
#	
#	print "Running fasta search vs Human... may be slow!"
#		
#	for x in range(1,32):
#		if numseqs==0:
#			break
#		blastdb="/data/blastdb/embl_con_hum-"+str(x)
#		os.system("rm blast.out")
#		print "Running fasta search vs "+blastdb
#		os.system("fasta34_t -T 4 -q -d 100 -b 1 -E 0.5 -H -m 9i primerseqs.fasta "+blastdb+" >> blast.out")
#	
#		print "Parsing fasta results"
#		currseq=""
#		readfile=open("blast.out","rU")
#		for line in readfile:
#		
#			if len(line)>3 and line[:3]==">>>":
#				#print line.strip()
#				currseq=line.split()[0][3:-1]
#				currseqlen=len(primers[currseq]["SEQUENCE"])
#				bestmatch=4
#			elif  len(line)>2 and line[:2]==">>":
#				currmatch=line.split()[0][2:]
#			elif len(line.split())>0 and line.split()[0]=="Smith-Waterman":
#				percent=float(line.split()[3][:-1])
#				overlap=int(line.split()[8])
#				nonmatches=currseqlen-int((percent/100)*overlap)
#				if nonmatches==bestmatch and nonmatches<4:
#					primers[primername][accession+"_matches"]+=1
#				if nonmatches<bestmatch:
#					primers[primername][accession+"_matches"]=1
#					primers[primername][accession+"_match_cost"]=nonmatches
#					bestmatch=nonmatches
#			elif len(line.split())>0 and line.split()[0]==currseq:
#					end=len(line.strip())
#					bits=line.strip()[len(currseq):].split(" ")
#					start=len(bits)-1+len(currseq)
#					#print line[start:end]
#					#print line
#					line=readfile.next()
#					#print line[start:end],
#					matches=0
#					cmatches=0
#					for x, hit in enumerate(line[end:start:-1]):
#						if hit==":":
#							matches+=1
#							if x<5:
#								cmatches+=1
#					
#					if (5-cmatches)<2:
#						#print "cmatch fail", (5-cmatches)
#						primers[currseq]["HUMAN_2_2_PASS"]=False
#						primers[currseq]["HUMAN_2_4_PASS"]=False
#					if ((end-start)-matches)<4:
#						primers[currseq]["HUMAN_4_PASS"]=False
#						primers[currseq]["HUMAN_2_4_PASS"]=False
#					if ((end-start)-matches)<primers[currseq]["HUMAN_MIN_DIFFS"]:
#						primers[currseq]["HUMAN_MIN_DIFFS"]=((end-start)-matches)
#						#print "match fail", ((end-start)-matches)
#				
#				#print matches, cmatches, , 5-cmatches
##				for x in print line[start:end]:
#		output=open("primerseqs.fasta","w")
#		numseqs=0
#		for x, primername in enumerate(primers.keys()):
#			if (primers[primername]["HUMAN_2_4_PASS"] or primers[primername]["HUMAN_2_2_PASS"]) and (primers[primername]["CONTAMINANT_2_4_PASS"] or primers[primername]["CONTAMINANT_2_2_PASS"]):
#				numseqs+=1
#				print >> output, ">"+primername
#				print >> output, primers[primername]["SEQUENCE"]
#			if x>99:
#				break
#		output.close()
#
#	handle=open("primers.tab", "w")
#	print_primers_to_tabfile(handle, primers)
#	handle.close()
#	print time.clock()-starttime
#	sys.exit()
#	
#	
#	
#	
#	
#	
#	
#	
#	
