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
	group.add_option("-a", "--assembler", action="store", type="choice", dest="assembler", choices=["velvetOptimiser","my_script","my_script-sc", "SOAPdenovo"], help="Assembly program to use (choose from my_script, velvetOptimiser or SOAPdenovo) [default= %default]", default="my_script")
	group.add_option("-n", "--number", action="store", type="int", dest="readspersplit", help="number of reads to put in each partition (0=use one partition) [default=%default]", default=0)
	group.add_option("-m", "--maxreads", action="store", type="int", dest="maxreads", help="number of reads from the fastq file to use for assembly (these will be randomly selected). This will override the -n option", default=0)
	group.add_option("-S", "--scaffold", action="store_true", dest="scaffold", help="Create scaffolds rather than contigs [default=%default]", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Read Filtering Options")
	group.add_option("-d", "--database", action="store", dest="database", help="database of contaminants to be removed (optional)", default="")
	group.add_option("-F", "--filter", action="store_true", dest="filter", help="Filter reads in fastq files for quality [default=%default]", default=False)
	group.add_option("-q", "--quality", action="store", type="int", dest="quality_cutoff", help="Quality cutoff to trim reads to [default=%default]", default=15)
	group.add_option("-l", "--length", action="store", type="int", dest="length_cutoff", help="Length cutoff: remove reads if at least one of the pair is less than this length [default=%default]", default=36)
	group.add_option("-v", "--save", action="store_true", dest="savefiltered", help="Save filtered fastq files [default=%default]", default=False)
	group.add_option("-g", "--gc", action="store", dest="gcfilter", help="Filter reads by GC (format must be a letter a (above) or b (below) followed by a percentage) [default=%default]", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Contig Filtering Options")
	group.add_option("-L", "--contiglength", action="store", type="int", dest="contig_length_cutoff", help="Contig length cutoff: remove contigs if they are less than this length (0 for no filtering) [default=%default]", default=100)
	group.add_option("-G", "--contiggc", action="store", dest="contiggcfilter", help="Filter contigs by GC (format must be a letter a (above) or b (below) followed by a percentage) [default=%default]", default=False)
	parser.add_option_group(group)
	
	group = OptionGroup(parser, "Reference options")
	#group.add_option("-R", "--reference", action="store", dest="reference", help="Reference file (for calculating reference length) (optional)", default="") #Need to add this functionality
	group.add_option("-R", "--reflength", action="store", type="int", dest="reflen", help="Length of reference [default=%default]", default=False)
	
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

	if options.database and not os.path.isfile(options.database):
		print "Cannot find file", options.database
		sys.exit()

	if options.gcfilter and not options.gcfilter[0] in ["a", "b"]:
		DoError("Invalid read GC cutoff. Must start with a (above) or b (below) followed by a float")
	if options.gcfilter:
		try:
			value=float(options.gcfilter[1:])
		except StandarError:
			DoError("Invalid read GC cutoff. Must start with a (above) or b (below) followed by a float")
	
	if options.maxreads<0:
		options.maxreads=0
	
	if options.contiggcfilter and not options.contiggcfilter[0] in ["a", "b"]:
		DoError("Invalid contig GC cutoff. Must start with a (above) or b (below) followed by a float")
	if options.contiggcfilter:
		try:
			value=float(options.contiggcfilter[1:])
		except StandarError:
			DoError("Invalid contig GC cutoff. Must start with a (above) or b (below) followed by a float")


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
    

###############################################
# Function to filter contigs for size and GC% #
###############################################

def filter_contigs(contigfile, outputname, l=100, gc=False):
	
	if gc and gc[0] not in ["a","b"]:
		DoError("Invalid gc a/b value")
	
	elif gc and gc[0]=='a':
		greater=True
		gcstring="above"
	else:
		greater=False
		gcstring="below"
	
	if gc:
		try:
			gcvalue=float(gc[1:])
		except StandardError:
			DoError("Invalid gc float value")
		
	
	if l>0 and not gc:
		print "Filtering contigs less than", l, "bases long."
	elif l>0 and gc:
		print "Filtering contigs less than", l, "bases long and with GC%", gcstring, gcvalue
	elif gc:
		print "Filtering contigs with GC%", gcstring, gcvalue
	else:
		return
		
		

	seq_records=[]
	
	for seq_record in SeqIO.parse(open(contigfile), "fasta"):
		seq=str(seq_record.seq)
		seqlen=len(seq.upper().replace("N",""))
		seqgc=(float(len(seq.upper().replace("A","").replace("T","")))/seqlen)*100
		
		
		keep=True
		
		if gc:
			if greater and seqgc>gcvalue:
				keep=False
			elif not greater and seqgc<gcvalue:
				keep=False
		
		if keep and l>0 and seqlen<l:
			keep=False
		
		if keep:
			seq_records.append(seq_record)
		
	SeqIO.write(seq_records, open(outputname,"w"), "fasta")


##############################
# Function to filter repeats #
##############################
def filter_repeats(contigfile, tmpname):
	ide = 97;
	print "Filtering contigs contained in larger contigs with at least "+str(ide)+"% match"
	
	seq_records=[]
	lengths={}
	
	for seq_record in SeqIO.parse(open(contigfile), "fasta"):
		seq_records.append(seq_record)
		lengths[seq_record.id]=len(seq_record.seq)

	if len(seq_records)==0:
		DoError("Amos has made an empty fasta file...bah")
		return
	
	os.system("nucmer --maxmatch --minmatch 200 --mincluster 200 --nosimplify --prefix="+tmpname+" "+contigfile+" "+contigfile+" >  /dev/null 2>&1")
	os.system("show-coords -r -T -o "+tmpname+".delta > "+tmpname+".coords")
	coords_lines=open(tmpname+".coords", "rU").readlines()[4:]

	matchlengths={}
	toremove=[]
	for line in coords_lines:
		words=line.strip().split()
	   
		if words[8] not in toremove and words[7]!=words[8] and lengths[words[8]]<lengths[words[7]]:

			if not matchlengths.has_key(words[7]):
				matchlengths[words[7]]={}
			if not matchlengths[words[7]].has_key(words[8]):
				matchlengths[words[7]][words[8]]=[0.0,0.0]

			matchlengths[words[7]][words[8]][0]+=int(words[5])
			matchlengths[words[7]][words[8]][1]+=(float(words[6])/100)*int(words[5])
		  
			
			if matchlengths[words[7]][words[8]][1]>(matchlengths[words[7]][words[8]][0]/100)*ide and matchlengths[words[7]][words[8]][0]>(float(lengths[words[8]])/100)*ide:
				toremove.append(words[8])
	   
			 

	
	
	count=1
	my_new_records=[]
	for record in seq_records:
		if record.id not in toremove:
			record.id="contig_"+str(count)
			my_new_records.append(record)
			count+=1
	SeqIO.write(my_new_records, open(contigfile,"w"), "fasta")
	

	return

	

def get_read_length(fastqfile):
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
	db=options.database
	
	
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
	
	
	if not shuffled and db!="" and os.path.isfile(db):
		print "Mapping reads vs", db
		os.system("cp "+db+" "+tmpname+"_db.fasta")
		os.system("bwa index "+tmpname+"_db.fasta")
		os.system(BWA_DIR+"bwa aln -n 3 -q 15 "+tmpname+"_db.fasta "+forward+" > "+tmpname+".F.sai")
		os.system(BWA_DIR+"bwa aln -n 3 -q 15 "+tmpname+"_db.fasta "+reverse+" > "+tmpname+".R.sai")
		#Join both aligments
		os.system(BWA_DIR+"bwa sampe "+tmpname+"_db.fasta "+tmpname+".F.sai "+tmpname+".R.sai "+forward+" "+reverse+" > "+tmpname+".sam")
		#os.system(SAMTOOLS_DIR+"samtools view -b -S "+tmpname+".sam > "+tmpname+".bam")
		print "Finding reads that do not map against", db
		remove_unmapping_reads(tmpname+".sam")
		#os.system("~sh16/scripts/get_unmapped_fastq_from_bam.py "+tmpname+".sam")
#		os.system("cp "+tmpname+"_unmapped.fastq "+forward.replace("_1.fastq","")+".fastq")
#		os.system("mv "+tmpname+"_unmapped.fastq "+tmpname+".fastq")
		os.system("rm "+tmpname+".sam")
		if not shuffled and options.assembler in ["velvetOptimiser","my_script","my_script-sc"]:
			print "Shuffling sequences"
			os.system("~sh16/scripts/shufflefastqSequences.pl "+tmpname+"_1_unmapped.fastq "+tmpname+"_2_unmapped.fastq "+tmpname+".fastq")
			os.system("rm "+tmpname+"_1_unmapped.fastq "+tmpname+"_2_unmapped.fastq")
		else:
			os.system("mv "+tmpname+"_1_unmapped.fastq "+tmpname+"_1.fastq")
			os.system("mv "+tmpname+"_2_unmapped.fastq "+tmpname+"_2.fastq")
		
		#remove the next 3 commands and put in the last 2
		#os.system("mv "+tmpname+"_unmapped.fastq "+forward.replace("_1.fastq","")+".fastq")
		#os.system("rm "+tmpname+"*")
		#sys.exit()
	
	elif not shuffled and options.assembler in ["velvetOptimiser","my_script","my_script-sc"]:
		print "Shuffling sequences"
		
		os.system("~sh16/scripts/shufflefastqSequences.pl "+forward+" "+reverse+" "+tmpname+".fastq")
	elif options.assembler in ["velvetOptimiser","my_script","my_script-sc"]:
		if shuffled.split(".")[-1]=="gz":
			os.system("zcat "+shuffled+" > "+tmpname+".fastq")
		else:
			os.system("cp "+shuffled+" "+tmpname+".fastq")
	elif not shuffled:
		if forward.split(".")[-1]=="gz":
			os.system("zcat "+forward+" > "+tmpname+"_1.fastq")
		else:
			os.system("cp "+forward+" "+tmpname+"_1.fastq")
		if reverse.split(".")[-1]=="gz":
			os.system("zcat "+reverse+" > "+tmpname+"_2.fastq")
		else:
			os.system("cp "+reverse+" "+tmpname+"_2.fastq")
		
	else:
		os.system("rm "+tmpname+"*")
		DoError("Cannot use "+options.assembler+" with shuffled data")
		
		
	
#	if db!="" and os.path.isfile(db):
#		print "Blasting vs contaminant database (this may take some time!)"
#		os.system("formatdb -p F -i "+db)
#		os.system("~sh16/scripts/fastq2fasta.pl "+tmpname+".fastq "+tmpname+".fasta")
#		os.system("blastall -p blastn -m 8 -b 1 -v 1 -e 1e-5 -i "+tmpname+".fasta -d "+db+" -o "+tmpname+".blast >  /dev/null 2>&1")
#		print "removing contaminant blast hits"
#		os.system("~sh16/scripts/remove_blast_hits_from_fastq.py "+tmpname+".fastq "+tmpname+".blast")
#		os.system("rm "+tmpname+".fasta formatdb.log "+tmpname+".blast")
#		os.system("mv "+tmpname+"_filtered.fastq "+tmpname+".fastq")
	
	
	#Filter reads for length and GC
	
	if options.filter and options.assembler in ["velvetOptimiser","my_script","my_script-sc"]:
		filter_fastq(tmpname+".fastq", reversefile=False, shuffled=True, quality_cutoff=options.quality_cutoff, length_cutoff=options.length_cutoff, GCcutoff=options.gcfilter)
		os.system("mv "+tmpname+"_filtered.fastq "+tmpname+".fastq")
		
	elif options.filter:
		filter_fastq(tmpname+"_1.fastq", reversefile=tmpname+"_2.fastq", shuffled=False, quality_cutoff=options.quality_cutoff, length_cutoff=options.length_cutoff, GCcutoff=options.gcfilter)
		os.system("cp "+tmpname+"_1_filtered.fastq Si_1.fastq")
		os.system("cp "+tmpname+"_2_filtered.fastq Si_2.fastq")
		os.system("mv "+tmpname+"_1_filtered.fastq "+tmpname+"_1.fastq")
		os.system("mv "+tmpname+"_2_filtered.fastq "+tmpname+"_2.fastq")
		
	
	
	
	
	#Count the number of reads in the files
	
	if options.assembler in ["velvetOptimiser","my_script","my_script-sc"]:
		linecount=bufcount(tmpname+".fastq")
		
	else:
		linecount1=bufcount(tmpname+"_1.fastq")
		linecount2=bufcount(tmpname+"_2.fastq")
		if linecount1!=linecount2:
			os.system("rm "+tmpname+"*")
			DoError("Forward and reverse line counts are not the same")
		linecount=linecount1+linecount2
			
	
	
	#print linecount, (linecount/4)/readspersplit
	 
	if readspersplit==0 or maxreads>0:
		numfiles=1
		if maxreads>0 and maxreads<=float(linecount)/8:
			readsperfile=maxreads/2
		else:
			readsperfile=int(float(linecount)/8)
	else:
		numfiles=(linecount/4)/readspersplit
		readsperfile=int((float(linecount)/8)/(numfiles))
	if numfiles<1:
		numfiles=1
	
	print "There are", linecount/8, "paired reads in your fastq file. Splitting fastq into", numfiles, "files containing approximately", readsperfile, "paired reads each"
	sys.stdout.flush()
	
	outfiles=[]
	print maxreads, readsperfile
	if maxreads>0 and readsperfile<float(linecount)/8:
		if options.assembler in ["velvetOptimiser","my_script","my_script-sc"]:
			outfile = open(tmpname+"_1.fastq","w")
			readstoadd=sample(xrange(linecount/8), readsperfile)
			
			readstoadd.sort()
			readsset=set(readstoadd)
			fastqfile=open(tmpname+".fastq","rU")
			linenum=0
			count=0
			for x in fastqfile:
				lines=[x.strip()]
				linenum+=1
				for y in range(0,7):
					lines.append(fastqfile.next().strip())
					linenum+=1
				
				if linenum/8 in readsset:
					for line in lines:
						count+=1
						print >> outfile, line
					readsset.remove(linenum/8)
					
			outfile.close()
			print count
		else:
			fastqfilef=open(tmpname+"_1.fastq","w")
			fastqfiler=open(tmpname+"_2.fastq","w")
			readstoadd=sample(xrange(linecount/4), readsperfile)
			readstoadd.sort()
			readsset=set(readstoadd)
			fastqfile=open(tmpname+".fastq","rU")
			linenum=0
			for x in fastqfile:
				linesf=[x.strip()]
				linesr=[]
				linenum+=1
				for y in range(0,3):
					linesf.append(fastqfilef.next().strip())
					linenum+=1
				for y in range(0,3):
					linesr.append(fastqfiler.next().strip())
				
				if linenum/4 in readsset:
					for line in linesf:
						print >> fastqfilef, line
					for line in linesr:
						print >> fastqfiler, line
					readsset.remove(linenum/8)
			fastqfilef.close()
			fastqfiler.close()
	
	elif options.assembler in ["velvetOptimiser","my_script","my_script-sc"]:
		fastqfile=open(tmpname+".fastq","rU")
		for x in range(0,numfiles):
			outfiles.append(open(tmpname+"_"+str(x+1)+".fastq","w"))
			
		for x in fastqfile:
		
			lines=[x.strip()]
				
			for y in range(0,7):
				lines.append(fastqfile.next().strip())
			
			rannum=randint(0, len(outfiles)-1)
			
			for line in lines:
				print >> outfiles[rannum], line
			
		
		for output in outfiles:
			output.close()
			
			
	else:
		fastqfilef=open(tmpname+"_1.fastq","rU")
		fastqfiler=open(tmpname+"_2.fastq","rU")
		for x in range(0,numfiles):
			outfiles.append([open(tmpname+"_"+str(x+1)+"_1.fastq","w"),open(tmpname+"_"+str(x+1)+"_2.fastq","w")])

		for x in fastqfilef:
		
			linesf=[x.strip()]
			linesr=[]
			
			for y in range(0,3):
				linesf.append(fastqfilef.next().strip())
			for y in range(0,3):
				linesr.append(fastqfiler.next().strip())
			
			rannum=randint(0, len(outfiles)-1)
			
			for line in linesf:
				print >> outfiles[rannum][0], line
			for line in linesr:
				print >> outfiles[rannum][1], line
			
		
		for output in outfiles:
			output[0].close()
			output[1].close()
	

	if options.assembler in ["velvetOptimiser","my_script","my_script-sc"]:
		if options.savefiltered and (options.filter or os.path.isfile(db)):
			os.system("mv "+tmpname+".fastq "+nameprefix+"_filtered.fastq")
		else:
			os.system("rm "+tmpname+".fastq")
	else:
		if options.savefiltered and (options.filter or os.path.isfile(db)):
			os.system("mv "+tmpname+"_1.fastq "+nameprefix+"_1_filtered.fastq")
			os.system("mv "+tmpname+"_2.fastq "+nameprefix+"_2_filtered.fastq")
		else:
			os.system("rm "+tmpname+"_1.fastq")
			os.system("rm "+tmpname+"_2.fastq")
		
		
		
		
	
	if options.assembler=="velvetOptimiser":
		outputlocation="velvet_data_*/contigs.fa"
		print "Using VelvetOptimiser.pl to assemble each chunk"
		
		maxkmer=int(float(readlength*3/4))
		if maxkmer>61:
			maxkmer=61
		
		if numfiles>1:
			if options.scaffold:
				velvetstring="/software/pathogen/external/applications/VelvetOptimiser-2.1.7/VelvetOptimiser.pl --p "+tmpname+"_${LSB_JOBINDEX}_velvet -e "+str(int(float(readlength*3/4)))+" -f \"-fastq -shortPaired "+tmpname+"_${LSB_JOBINDEX}.fastq\""
			else:
				velvetstring="/software/pathogen/external/applications/VelvetOptimiser-2.1.7/VelvetOptimiser.pl --p "+tmpname+"_${LSB_JOBINDEX}_velvet -e "+str(int(float(readlength*3/4)))+" -f \"-fastq -shortPaired "+tmpname+"_${LSB_JOBINDEX}.fastq\" --o \"-scaffolding no\""
		else:
			if options.scaffold:
				velvetstring="/software/pathogen/external/applications/VelvetOptimiser-2.1.7/VelvetOptimiser.pl --p "+tmpname+"_1_velvet --e "+str(int(float(readlength*3/4)))+" -f \"-fastq -shortPaired "+tmpname+"_1.fastq\""
			else:
				velvetstring="/software/pathogen/external/applications/VelvetOptimiser-2.1.7/VelvetOptimiser.pl --p "+tmpname+"_1_velvet --e "+str(int(float(readlength*3/4)))+" -f \"-fastq -shortPaired "+tmpname+"_1.fastq\" --o \"-scaffolding no\""
		#print velvetstring
		
	elif options.assembler in ["my_script","my_script-sc"]:
		outputlocation="velvet/contigs.fa"
		print "Using velvet_assembly.sh to assemble each chunk"
		kmC=15
		if options.reflen:
			maxcoverage=(readlength*int((float(linecount)/4)/(numfiles)))/options.reflen
			#kmC = maxcoverage * (readlength - (20)) / readlength
			kmC = maxcoverage * ((readlength+1)/2) / readlength
			#print maxcoverage, kmC,
		if kmC<12:
			kmC=12
		elif kmC>20:
			kmC=20
		#print kmC
		
		kmC=10
		
		if options.assembler=="my_script-sc":
			wga="-w"
		else:
			wga=""
		
		if numfiles>1:
			if options.scaffold:
				velvetstring="~sh16/scripts/velvet_assembly.sh -f "+tmpname+"_${LSB_JOBINDEX}.fastq "+wga+"  -p -e "+str(kmC)+" -i 250"
			else:
				velvetstring="~sh16/scripts/velvet_assembly.sh -f "+tmpname+"_${LSB_JOBINDEX}.fastq "+wga+"  -n -p -e "+str(kmC)+" -i 250"
				
		else:
			if options.scaffold:
				velvetstring="~sh16/scripts/velvet_assembly.sh -f "+tmpname+"_1.fastq "+wga+"  -p -e "+str(kmC)+" -i 250"
			else:
				velvetstring="~sh16/scripts/velvet_assembly.sh -f "+tmpname+"_1.fastq "+wga+"  -n -p -e "+str(kmC)+" -i 250"
			
	elif options.assembler=="SOAPdenovo":
		if options.scaffold:
			outputlocation="SOAP.scafSeq"
		else:
			outputlocation="SOAP.contig"
		print "Using SOAPdenovo to assemble each chunk"
#		kmC=15
#		if options.reflen:
#			maxcoverage=(readlength*int((float(linecount)/4)/(numfiles)))/options.reflen
#			#kmC = maxcoverage * (readlength - (20)) / readlength
#			kmC = maxcoverage * ((readlength+1)/2) / readlength
#			#print maxcoverage, kmC,
#		if kmC<12:
#			kmC=12
#		elif kmC>20:
#			kmC=20
		#print kmC
		
		for x in range(1,numfiles+1):
			pwd=os.getcwd()
			contigout=open(tmpname+"_configfile_"+str(x),"w")
			print >> contigout, "#maximal read length"
			print >> contigout, "max_rd_len="+str(readlength)
			print >> contigout, "[LIB]"
			print >> contigout, "avg_ins=300"
			print >> contigout, "reverse_seq=0"
			print >> contigout, "asm_flags=3"
			print >> contigout, "rank=1"
			print >> contigout, "map_len=32"
			print >> contigout, "q1="+pwd+"/"+tmpname+"_"+str(x)+"_1.fastq"
			print >> contigout, "q2="+pwd+"/"+tmpname+"_"+str(x)+"_2.fastq"
			contigout.close()
		
		
		#could add kmer stuff here
		if numfiles>1:
			velvetstring=SOAPdir+"SOAPdenovo all -s "+tmpname+"_configfile_${LSB_JOBINDEX} -o "+tmpname+"_${LSB_JOBINDEX}_SOAP -K 31 -d 2 -D 2 -p 1"
		else:
			velvetstring=SOAPdir+"SOAPdenovo all -s "+tmpname+"_configfile_1 -o "+tmpname+"_1_SOAP -K 31 -d 2 -D 2 -p 1"

	
	if numfiles>1:
		os.system('echo \''+velvetstring+'\' | bsub -J "'+tmpname+'[1-'+str(numfiles)+']" > '+tmpname+'jobstring')
	
		
		jobnum=open(tmpname+'jobstring', "rU").read().split(">")[0].split("<")[1]
	
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
	else:
		os.system(velvetstring)
	
	for x in range(1,numfiles+1):
		if options.assembler in ["velvetOptimiser","my_script"]:
			os.system("rm "+tmpname+"_"+str(x)+".fastq")
		else:
			os.system("rm "+tmpname+"_"+str(x)+"_1.fastq "+tmpname+"_"+str(x)+"_2.fastq")
	
	print

	
	if numfiles>1:
		print "\nJoining contigs with Amos"
		
		os.system("mkdir "+tmpname+"_joined")
		os.system("cp "+tmpname+"_1_"+outputlocation+" "+tmpname+"_joined/current_assembly.fa")
		
		pwd=os.getcwd()
		
		for x in range(2,numfiles+1):
			os.system("cp "+tmpname+"_"+str(x)+"_"+outputlocation+" "+tmpname+"_joined/contigs.fa")
			os.chdir(pwd+"/"+tmpname+"_joined/")
			refcount=os.popen("grep -c '^>' current_assembly.fa").readlines()[0].strip()
			os.system("cat current_assembly.fa contigs.fa > both.fa")
			os.system("toAmos -s both.fa -o both.afg")
			os.system("minimus2 both -D REFCOUNT="+refcount)
			#os.system("minimus2 both -D REFCOUNT=0")
			os.system("cat both.fasta both.singletons.seq > current_assembly.fa")
			
			
			
			filter_repeats("current_assembly.fa", tmpname)
			
			
			os.system('~as9/bin/new.stats current_assembly.fa | grep "total_length\|number\|mean_length\|longest\|n50"')
			os.system("rm -rf both* contigs.fa")
			
	#		os.system("/nfs/users/nfs_a/as9/bin/maximus_test.pl -1 joined.fasta -2 contigs.fa -n 1 -o temp -l 200 -m 200")
	#		os.system("cat results*/temp_all_merged.fasta results*/temp.singletons.seq > joined.fasta")
	#		os.system("rm -rf joined.fasta.* contigs.fa* ite1 ite2 results*")
			os.chdir(pwd)
		
		if options.contig_length_cutoff!=0 or options.contiggcfilter:
			filter_contigs(tmpname+"_joined/current_assembly.fa", outputname, l=options.contig_length_cutoff, gc=options.contiggcfilter)
		else:	
			os.system("mv "+tmpname+"_joined/current_assembly.fa "+outputname)
	else:
		if options.contig_length_cutoff!=0 or options.contiggcfilter:
			os.system("mv "+tmpname+"_1_"+outputlocation+" "+tmpname+"_1")
			filter_contigs(tmpname+"_1", outputname, l=options.contig_length_cutoff, gc=options.contiggcfilter)
		else:		
			os.system("mv "+tmpname+"_1_"+outputlocation+" "+outputname)
	
	os.system("sed 's/"+tmpname+"/"+outputname.split("/")[-1].split(".")[0]+"/g' "+outputname+" > "+tmpname)
	os.system("mv "+tmpname+" "+outputname)
	os.system("rm -rf "+tmpname+"*")
	print "Final stats:"
	os.system('~as9/bin/new.stats '+outputname+' | grep "total_length\|number\|mean_length\|longest\|n50"')
	print "Done"
