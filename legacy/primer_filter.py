#!/usr/bin/env python

import os, sys
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_SeqIO import *
from Si_general import *
from Bio import SeqIO
from Bio.SeqUtils import GC
import gzip
from optparse import OptionParser, OptionGroup


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
	group.add_option("-i", "--input", action="store", dest="tabfile", help="Input primer tab file name", default="", metavar="FILE")
	group.add_option("-p", "--prefix", action="store", dest="prefix", help="Prefix for output files", default="")
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
	
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.tabfile=='':
		DoError('No primer tab file selected')
	elif not os.path.isfile(options.tabfile):
		DoError('Cannot find file '+options.tabfile)
	
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

def make_primer3_input_file(handle, sequence1, starting_pos, sequence2="", opt_size=17, min_size=15, max_size=24, opt_tm=42.5, min_tm=40, max_tm=45, min_gc=20, opt_gc=50, max_gc=80, gc_clamp=1, poly_x=5, self_any=8.00, self_end=3.00):
	
	
	
	print >> handle, "SEQUENCE_ID="+str(starting_pos)
	
	if sequence2=="":
		print >> handle, "PRIMER_TASK=pick_detection_primers"
		print >> handle, "SEQUENCE_TEMPLATE="+sequence1
		print >> handle, "PRIMER_PICK_LEFT_PRIMER=1"
		print >> handle, "PRIMER_PICK_INTERNAL_OLIGO=0"
		print >> handle, "PRIMER_PICK_RIGHT_PRIMER=0"
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
	print >> handle, "PRIMER_NUM_RETURN=5"
	print >> handle, "PRIMER_MAX_POLY_X="+str(poly_x)
	print >> handle, "PRIMER_SELF_ANY="+str(self_any)
	print >> handle, "PRIMER_SELF_END="+str(self_end)
	print >> handle, "PRIMER_TM_FORMULA=1"
	print >> handle, "PRIMER_SALT_CORRECTIONS"
	#print >> handle, "PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION="+str(len(sequence))+",20 16,20"
	print >> handle, "="
	
	
	
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
			
		if phit==rhit and phit==lhit and sequence1==lseq and sequence2==rhit and COMP_ANY!="" and COMP_END!="":
			return COMP_ANY, COMP_END
				
	
#	if COMP_ANY=="" or COMP_END=="":
#		print sequence1, sequence2
#		for line in primer3out:
#			print line.strip()
#		sys.exit()
	
	if COMP_ANY=="":
		COMP_ANY=float(len(sequence1))
	if COMP_END=="":
		COMP_END=5.0
	return COMP_ANY, COMP_END

################
# Main program #
################		

if __name__ == "__main__":


	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	
	#make random name for temporary files
	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))	
	

	try:
		emblrecord=open_annotation(options.tabfile)
	except StandardError:
		DoError("Cannot open annotation file "+options.tabfile+" please check the format")
	
	
	
	genomelength=len(emblrecord.seq)
	
	genomeblocks={}
	blockstarts=[]
	
	for x in range(0,genomelength,10000):
		genomeblocks[x]=0
		blockstarts.append(x)
	
	#print blockstarts
	#sys.exit()
	newfeatures=[]
	toremove=[]
	
	for seqFeature in emblrecord.features:
		if seqFeature.type=="primer":
			
			if float(seqFeature.qualifiers["HUMAN_MAX_MISPRIMING"][0])<14 and float(seqFeature.qualifiers["CONTAMINANT_MAX_MISPRIMING"][0])<14:
				x=0
				#print x, blockstarts[x],genomelength,  seqFeature.location.nofuzzy_start
				while x<len(blockstarts) and int(seqFeature.location.nofuzzy_start)>blockstarts[x]:
					#print x, blockstarts[x],genomelength,  seqFeature.location.nofuzzy_start
					x+=1
				genomeblocks[blockstarts[x-1]]+=1
				seqFeature.qualifiers["block"]=blockstarts[x-1]
			else:
				toremove.append(seqFeature)
			
	
		else:
			toremove.append(seqFeature)
	for seqFeature in toremove:
		emblrecord.features.remove(seqFeature)
	
	
	print len(emblrecord.features), "primers match specificity requirements\nChecking for primer dimers..."
	finalprimers=[]
	
	count=0.0
	hundredth=float(len(emblrecord.features))/100
	
	
	
	print genomeblocks
	#sys.exit()
	discarded=0
	kept=0
	
	for x, seqFeature in enumerate(emblrecord.features):
		count=count+1
		print kept, "kept,", discarded, "discarded\r",
		#print count 
		if count==50:
			print genomeblocks
			count=0
		sys.stdout.flush()
		toremove=[]
		keep=True
		for seqFeature2 in emblrecord.features[x+1:]:
			handle=open(tmpname+"primer_3_input", "w")
						
			make_primer3_input_file(handle, seqFeature.qualifiers["SEQUENCE"][0], seqFeature.location.nofuzzy_start, seqFeature2.qualifiers["SEQUENCE"][0], opt_size=options.opt_size, min_size=options.min_size, max_size=options.max_size, opt_tm=options.opt_tm, min_tm=options.min_tm, max_tm=options.max_tm, min_gc=options.min_gc, opt_gc=options.opt_gc, max_gc=options.max_gc, gc_clamp=options.gc_clamp, poly_x=options.poly_x, self_any=25, self_end=25)
			handle.close()
			COMP_ANY, COMP_END=run_primer3_dimer_check(tmpname+"primer_3_input", seqFeature.qualifiers["SEQUENCE"], seqFeature2.qualifiers["SEQUENCE"])
			
#			if seqFeature.qualifiers["SEQUENCE"][0] in ['TAGTCGTACTGCTCG', 'GCAGTACGACTACTTG'] or seqFeature2.qualifiers["SEQUENCE"][0] in ['TAGTCGTACTGCTCG', 'GCAGTACGACTACTTG']:
#				print count, COMP_ANY, COMP_END, seqFeature.qualifiers["SEQUENCE"][0], seqFeature2.qualifiers["SEQUENCE"][0]
			if float(COMP_ANY)>8 or float(COMP_END)>3:
#				print primer1, primer2, "form dimer"
				winner=None

				if genomeblocks[seqFeature.qualifiers["block"]]>genomeblocks[seqFeature2.qualifiers["block"]]:
					winner='primer2'
				elif genomeblocks[seqFeature.qualifiers["block"]]<genomeblocks[seqFeature2.qualifiers["block"]]:
					winner='primer1'
				elif float(seqFeature.qualifiers["PENALTY"][0])<float(seqFeature2.qualifiers["PENALTY"][0]):
					winner='primer1'
				elif float(seqFeature.qualifiers["PENALTY"][0])>float(seqFeature2.qualifiers["PENALTY"][0]):
					winner='primer2'
				
#				if seqFeature.qualifiers["SEQUENCE"][0] in ['TAGTCGTACTGCTCG', 'GCAGTACGACTACTTG'] or seqFeature2.qualifiers["SEQUENCE"][0] in ['TAGTCGTACTGCTCG', 'GCAGTACGACTACTTG']:
#					print "Winner =", winner, genomeblocks[seqFeature.qualifiers["block"]], genomeblocks[seqFeature2.qualifiers["block"]], float(seqFeature.qualifiers["PENALTY"][0]), float(seqFeature2.qualifiers["PENALTY"][0])
				
				
				
				if winner=='primer2':
					genomeblocks[seqFeature.qualifiers["block"]]-=1
					keep=False
					break
				elif winner=='primer1':
					genomeblocks[seqFeature2.qualifiers["block"]]-=1
					toremove.append(seqFeature2)
					
		if keep:
			#print seqFeature
			kept+=1
			finalprimers.append(seqFeature)
			for loser in toremove:
				discarded+=1
				emblrecord.features.remove(loser)
		else:
			for loser in toremove:
				genomeblocks[loser.qualifiers["block"]]+=1
			discarded+=1
				
		
		
	#print len(goodprimers), primernum

	print len(finalprimers), "primers found that do not form primer-dimers"
	
	handle =open(options.prefix+"_filtered_primers.gb","w")
	SeqIO.write(finalprimers,handle, "gb")
	handle.close()
	os.system("rm "+tmpname+"*")