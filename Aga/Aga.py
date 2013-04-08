#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence

##################
# Import modules #
##################

import string, re
import os, sys, random, math, time
from optparse import OptionParser, OptionGroup
from random import *
from Bio import SeqIO
from Bio.Seq import Seq


##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()

##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-r", "--reference", action="store", dest="ref", help="Reference fasta", default="", metavar="FILE")
	group.add_option("-e", "--embl", action="store", dest="embl", help="Reference embl fasta", default="", metavar="FILE")
	group.add_option("-t", "--tab", action="store", dest="tab", help="Reference tab file of MGEs", default="", metavar="FILE")
	group.add_option("-o", "--output_prefix", action="store", dest="prefix", help="Output prefix", default="")
	group.add_option("-m", "--no_map", action="store_false", dest="mapping", help="Do not map data against reference using smalt first [default = do mapping]", default=True)
	group.add_option("-s", "--single", action="store_true", dest="single", help="Fastq files are not paired", default=False)
	parser.add_option_group(group)
#	group.add_option("-d", "--contaminant_database", action="store", dest="contaminants", help="Name file containing contaminant accession numbers", default=False, metavar="FILE")
#	group.add_option("-H", "--human", action="store_true", dest="human", help="Blast primers against human genome", default=False)
#	parser.add_option_group(group)
#if len(sys.argv)!=3:
#	print "Usage: create_pan_genome.py <ssaha_folders>"
#	sys.exit()


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.ref=='':
		DoError('No reference file selected')
	elif not os.path.isfile(options.ref):
		DoError('Cannot find file '+options.ref)
	if options.embl!='' and not os.path.isfile(options.embl):
		DoError('Cannot find file '+options.embl)
	if options.prefix=='':
		DoError('No output prefix specified')
	
	
	return
	
##############################
# Function to filter repeats #
##############################
def filter_repeats(contigfile, tmpname):
	ide = 90;
	print "Filtering contigs contained in larger contigs with at least "+str(ide)+"% match"
	
	seq_records=[]
	lengths={}
	
	for seq_record in SeqIO.parse(open(contigfile), "fasta"):
		seq_records.append(seq_record)
		lengths[seq_record.id]=len(seq_record.seq)

	if len(seq_records)==0:
		DoError("Amos has made an empty fasta file...bah")
		return
	
	os.system("nucmer --maxmatch --minmatch 200 --mincluster 200 --nosimplify --prefix="+tmpname+" "+contigfile+" "+contigfile)#+" >  /dev/null 2>&1")
	os.system("show-coords -r -T -o "+tmpname+".delta > "+tmpname+".coords")
	coords_lines=open(tmpname+".coords", "rU").readlines()[4:]

	matchlengths={}
	toremove=[]
	for line in coords_lines:
		words=line.strip().split()
	   
		if words[8] not in toremove and words[7]!=words[8] and lengths[words[8]]<lengths[words[7]]:

			if not words[7] in matchlengths:
				matchlengths[words[7]]={}
			if not words[8] in matchlengths[words[7]]:
				matchlengths[words[7]][words[8]]=[0.0,0.0]

			matchlengths[words[7]][words[8]][0]+=int(words[5])
			matchlengths[words[7]][words[8]][1]+=(float(words[6])/100)*int(words[5])
		  
			
			if matchlengths[words[7]][words[8]][1]>(matchlengths[words[7]][words[8]][0]/100)*ide and matchlengths[words[7]][words[8]][0]>(float(lengths[words[8]])/100)*ide:
				toremove.append(words[8])
	   
			 

	
	
	count=1
	my_new_records=[]
	for record in seq_records:
		if record.id not in toremove:
			#record.id="contig_"+str(count)
			my_new_records.append(record)
			count+=1
	SeqIO.write(my_new_records, open(contigfile,"w"), "fasta")
	

	return

################
# Main program #
################		

if __name__ == "__main__":
	
	starttime=time.clock()

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	chars = string.ascii_letters + string.digits
	
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	
	
	acc_file=tmpname+"_acc.mfa"
	core_file=tmpname+"_core.mfa"
	
	
	if options.tab!="":
		try:
			tablines=open(options.tab,"rU").readlines()
		except StandardError:
			print "Failed to read tab file"
			sys.exit()
		ftloc=-1
		note=""
		regions=[]
		for line in tablines:
#			         111111111122222222223333333333
#			123456789012345678901234567890123456789
#			FT   misc_feature    2825393..2826315
			line=line.strip()
			if len(line)>22 and line[:2]=="FT" and len(line[2:21].strip())>0:
				try: ftloc=map(int,(line[21:].split("..")))
				except StandardError:
					print line
					sys.exit()
				note=""
			elif len(line)>22 and line[:2]=="FT" and len(line[2:21].strip())==0 and ftloc!=-1:
				if line[21:].split("=")[0]=="/note":
					note='='.join(line[22:].split("=")[1:])
					if len(ftloc)==1:
						regions.append([ftloc[0], ftloc[0], note])
					else:
						regions.append([ftloc[0], ftloc[1], note])
						
					ftloc=-1

		
		refseq=''.join(open(options.ref, "rU").read().split('\n')[1:])
		
		cur=0
		count=1
		core_count=1
		acc_count=1
		core_out=open(core_file, "w")
		acc_out=open(acc_file, "w")
		
		for region in regions:
			if region[0]-cur>1:
				print "ref", cur, "..",region[0]-1
				print >> core_out, ">Ref_"+str(count), cur, "..",region[0]-1
				print >> core_out, refseq[cur+1:region[0]]
				core_count+=1
				count+=1
				
			print "acc", region[0], "..", region[1], region[2]
			print >> acc_out, ">Ref_"+str(count), region[0], "..",region[1]-1, region[2]
			print >> acc_out, refseq[region[0]:region[1]]
			cur=region[1]
			acc_count+=1
			count+=1
		if cur<len(refseq):
			print >> core_out, ">Ref_"+str(count), cur, "..",len(refseq)
			print >> core_out, refseq[cur:]
		core_out.close()
		acc_out.close()
		
	else:
		os.system("cp "+options.ref+" "+core_file)
			
	sys.exit()

	if options.mapping:
		os.system("~sh16/scripts/multiple_mappings_to_bam.py -M 2 -z 0.9 -p smalt -v latest -X -r "+core_file+" "+' '.join(args)+' > '+tmpname+'jobstring')
		
		jobnum=open(tmpname+'jobstring', "rU").readlines()[-4].split(">")[0].split("<")[1]
		
		
		print "Mapping reads against reference"
		
		todo=1
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
			print exited, "jobs exited without finishing"
			sys.exit()
	
	
	

	
	depths=[]
	donefolders=[]
	for arg in args:
		folder="_".join('.'.join(arg.split("/")[-1].split(".")[:-1]).split("_")[:-1])
		if folder in donefolders:
			continue
		donefolders.append(folder)
		
		if options.mapping or (not os.path.isfile(folder+"_SMALT/"+folder+"_unmapped_1.fastq") and not os.path.isfile(folder+"_SMALT/"+folder+"_unmapped_1.fastq.gz")) or (not os.path.isfile(folder+"_SMALT/"+folder+"_unmapped_2.fastq") and not os.path.isfile(folder+"_SMALT/"+folder+"_unmapped_2.fastq.gz")):
			print "Identifying reads not in proper pairs for", folder
			sys.stdout.flush()
			os.system("~sh16/scripts/Aga/bam_filter.py -t atleastoneunmapped -b "+folder+"_SMALT/"+folder+".bam -o "+folder+"_SMALT/"+folder+"_unmapped")

		contigdepths=os.popen("~sh16/scripts/Aga/contig_stats.py -b "+folder+"_SMALT/"+folder+".bam -H" )#.readlines()
		totlen=0.0
		totdepth=0.0
		for contigdepth in contigdepths:
			words=contigdepth.strip().split()
			totlen+=float(words[1])
			totdepth+=float(words[1])*float(words[2])

		exp=totdepth/totlen
		print folder, "mean depth =", exp
		sys.stdout.flush()
		depths.append([exp, folder])

		contigdepths.close()
	
	depths.sort()
	#depths.reverse()


	for strains in depths:
	
		exp=strains[0]
		folder=strains[1]
		
		print folder, exp
		
		if exp >40:
			exp=40
		elif exp<30:
			exp=30
		
		
		
		#outfile1=open(tmpname+"_1.fastq","w")
		#outfile2=open(tmpname+"_2.fastq","w")
		
		if os.path.isfile(folder+"/"+folder+"_SMALT/"+folder+"_unmapped_1.fastq.gz"):
			os.system("gunzip -f "+folder+"_SMALT/"+folder+"_unmapped_1.fastq.gz")
		if os.path.isfile(folder+"/"+folder+"_SMALT/"+folder+"_unmapped_2.fastq.gz"):
			os.system("gunzip -f "+folder+"_SMALT/"+folder+"_unmapped_2.fastq.gz")
		
		if os.path.isfile(folder+"_SMALT/"+folder+"_unmapped_1.fastq") and os.path.isfile(folder+"_SMALT/"+folder+"_unmapped_2.fastq"):
			os.system("cp "+folder+"_SMALT/"+folder+"_unmapped_1.fastq "+tmpname+"_1.fastq")
			os.system("cp "+folder+"_SMALT/"+folder+"_unmapped_2.fastq "+tmpname+"_2.fastq")
		else:
			print "Cannot find "+folder+"_SMALT/"+folder+"_unmapped.fastq"
			sys.exit()

		
		#os.system("~sh16/scripts/multiple_mappings_to_bam.py -r "+acc_file+" -f -X -y -p smalt -L "+tmpname+"_[12].fastq")# >  /dev/null 2>&1")

		os.system("samtools faidx "+acc_file)
		os.system("smalt index -k 13 -s 1 "+acc_file+".index "+acc_file)
		os.system("smalt map -y 0.9 -r "+str(randrange(1,99999))+" -f samsoft -o "+tmpname+".sam "+acc_file+".index "+tmpname+"_1.fastq "+tmpname+"_2.fastq")
		os.system("samtools view -b -S "+tmpname+".sam -t "+acc_file+".fai > "+tmpname+".1.bam")
		os.system("samtools sort "+tmpname+".1.bam "+tmpname)
		os.system("samtools index "+tmpname+".bam")
		os.system("rm -f "+tmpname+".1.bam "+tmpname+".sam "+acc_file+".fai "+acc_file+".index.*")
		
		
		contigstats=os.popen("~sh16/scripts/Aga/contig_stats.py -b "+tmpname+".bam -H" )#.readlines()
		contiglist=[]
		
		for contigstat in contigstats:
			words=contigstat.strip().split()
			print words
			if float(words[7])<90 or float(words[2])<(0.2*strains[0]):
				contiglist.append(words[0])
				
		contigstats.close()
		
		print "Keeping all reads that map to", contiglist
		
		if len(contiglist)>0:
			contigstring=','.join(contiglist)
			print contigstring
			os.system("~sh16/scripts/Aga/bam_filter.py -t aga -c "+contigstring+" -b "+tmpname+".bam -o "+tmpname+"_unmapped")
		else:
			os.system("~sh16/scripts/Aga/bam_filter.py -t atleastoneunmapped -b "+tmpname+".bam -o "+tmpname+"_unmapped")
			
		
		
		#os.system('~sh16/scripts/iterative_assembler.py -L 500 -n 0 -f '+tmpname+"_unmapped_1.fastq -r "+tmpname+"_unmapped_2.fastq -o "+folder+"_contigs.mfa")
		os.system('~sh16/scripts/velvet_assembly.sh -p -n -i 300 -e '+str(exp)+' -f '+tmpname+"_unmapped_1.fastq -r "+tmpname+"_unmapped_2.fastq -s "+folder+".fastq")
		os.system("rm -rf "+tmpname+".fasta "+folder+"_contigs.mfa")
		os.system("cp "+folder+"_velvet/contigs.fa "+folder+"_contigs.mfa")
		os.system("rm -rf "+folder+"_velvet")
		os.system("cat "+core_file+" "+acc_file+" > "+tmpname+".fasta")
		os.system("~sh16/scripts/fasta2fastq_shredder.py "+tmpname+".fasta "+tmpname+" 76 3 l 250")
		os.system("rm -rf "+tmpname+".fasta")
		os.system("smalt index -k 13 -s 1 "+folder+"_contigs.mfa.index "+folder+"_contigs.mfa")
		os.system("smalt map -r "+str(randrange(1,99999))+" -f samsoft -o "+tmpname+".sam "+folder+"_contigs.mfa.index "+tmpname+"_1.fastq "+tmpname+"_2.fastq")
		os.system("samtools view -b -o "+tmpname+".bam -t "+folder+"_contigs.mfa -S "+tmpname+".sam")
		os.system("samtools sort "+tmpname+".bam "+tmpname+"_sort")
		os.system("samtools index "+tmpname+"_sort.bam")
		
		os.system("~sh16/scripts/Contig_summary.py "+folder+"_contigs.mfa")
		os.system("~sh16/scripts/Contig_summary.py "+acc_file)
		
		contigstats=os.popen("~sh16/scripts/Aga/contig_stats.py -b "+tmpname+"_sort.bam -H" )#.readlines()
		os.system("rm -rf "+tmpname+".[sb]am")
		contiglist=[]
		
		for contigstat in contigstats:
			words=contigstat.strip().split()
			print words
			if float(words[7])<90:# or float(words[2])<(0.2*strains[0]):
				contiglist.append(words[0])
				
		contigstats.close()
		output=open(tmpname+".fasta", "w")
		lines=open(folder+"_contigs.mfa", "rU").read().split(">")[1:]
		for line in lines:
			seqname=line.split("\n")[0].split()[0]
			seq=''.join(line.split("\n")[1:])
			if seqname in contiglist and len(seq)>500:
				print >> output, ">"+seqname
				print >> output, ''.join(line.split("\n")[1:])
		output.close()
		os.system("cat "+tmpname+".fasta >> "+acc_file)
#		os.system("~sh16/scripts/Contig_summary.py "+folder+"_contigs.mfa "+acc_file)
#		filter_repeats(acc_file, tmpname)
		os.system("~sh16/scripts/Contig_summary.py "+acc_file)
	
	
	os.system("mv "+acc_file+" "+options.prefix+"_accessory_genome.fasta")	
	os.system("cat "+core_file+" "+options.prefix+"_accessory_genome.fasta > "+options.prefix+"_pan_genome.fasta")
	#filter_repeats(options.prefix+"_pan_genome.fasta", tmpname)
	os.system("rm -rf "+tmpname+"*")
	
	#Make pseudofastq for reference
	if not os.path.isfile('.'.join(options.ref.split('.')[:-1])+"_1.fastq") and not os.path.isfile('.'.join(options.ref.split('.')[:-1])+"_2.fastq"):
		refname='.'.join(options.ref.split('.')[:-1])
	else:
		refname=tmpname
	
	os.system("~sh16/scripts/fasta2fastq_shredder.py "+options.ref+" "+refname+" 76 3 c 250")
		
	#final mapping of all isolates against ref+accessory
	if options.embl!="":
		os.system("~sh16/scripts/multiple_mappings_to_bam.py -p smalt -v latest -G -O "+options.prefix+" -E -f -x -r "+options.prefix+"_pan_genome.fasta -e "+options.embl+" "+' '.join(args)+' '+refname+'_[12].fastq')
	else:
		os.system("~sh16/scripts/multiple_mappings_to_bam.py -p smalt -v latest -G -O "+options.prefix+" -E -f -x -r "+options.prefix+"_pan_genome.fasta "+' '.join(args)+' '+refname+'_[12].fastq')
	
	
	
	
	
	
