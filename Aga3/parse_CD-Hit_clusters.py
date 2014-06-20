#!/usr/bin/env python


import string
import os, sys
from optparse import OptionParser


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] <list of mfa files>"
	parser = OptionParser(usage=usage)
	parser.add_option("-f", "--fasta file", action="store", dest="fasta", help="Input fasta file", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file prefix. [Default = %default]", default="Resistome")
	parser.add_option("-m", "--minmatch", action="store", dest="minmatch", help="minimum match length to retain", default=25.0, type="float")
	#Could add more options in here so people can specify similarities etc.
	
	
	return parser.parse_args()

def DoError(errorstring):
	print "Error:", errorstring
	sys.exit()


def parse_CDHitEST_clstr_file(filehandle):
	clusters={}
	clusternumber=0
	
	print "Reading clusters from CD-Hit Output..."
	
	for line in filehandle:
		line=line.strip()
		if len(line.split())==2 and line.split()[0]==">Cluster":
			try:
				clusternumber=int(line.split()[1])+1
				clusters[clusternumber]=[]
			except ValueError:
				print "Error: Expecting integer in >Cluster line of CD-Hit cluster file."
				sys.exit()
		
		else:
			if clusternumber==0 or not clusternumber in clusters:
				print "Expecting line starting with >Cluster"
				print line
				sys.exit()
			
			words=line.split("...")[0].split()
			if line.split()[-1]=="*":
				strand="+"
				ref_seq=True
				percentid=100
			else:
				ref_seq=False
				if len(line.split()[-1].split('/'))!=3:
				       print "Cannot read strand/percentid"
				       print line
				       sys.exit()
				strand=line.split()[-1].split('/')[-2]

				if not strand in ["-", "+"]:
					print "Cannot read strand"
					print line
					sys.exit()
				try:
				       percentid=float(line.split()[-1].split('/')[-1].replace("%",""))
				except ValueError:
				       print "Cannot read percentid"
				       print line
				       sys.exit()
			
			if len(words)>2 and words[2][0]==">":
				clusters[clusternumber].append({"name":words[2][1:], "strand":strand, "key":ref_seq, "percent_id":percentid})
				
			else:
				print "Error: Expecting to find gene name in third column of line"
		
		
	print "Found", len(clusters), "clusters"

	return clusters


################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	

	if options.fasta=="":
		DoError("Fasta file (-f) required")
	elif not os.path.isfile(options.fasta):
		DoError("Cannot find file "+options.fasta)
	
		
	fastaname=options.fasta
	

#	print "Clustering reads to genes with CD-Hit-est..."
#	
#	returnval=os.system(' '.join(["cd-hit-est -i", options.db, "-G 0 -M 0",  "-o", options.output+"_clusters -d 10000"]))
#	
#	if returnval!=0 or not os.path.isfile(options.output+"_clusters.clstr"):
#		print "Error: Clustering failed. Please check the format of your sequence database file."
#		sys.exit()
	
	clusters=parse_CDHitEST_clstr_file(open(options.output+"_clusters.clstr", "rU"))
	
	
	print clusters
	sys.exit()
	
	reads=set([])
	genes=set([])
	
	for cluster in clusters:
		if len(clusters[cluster])>1:
			for x in clusters[cluster]:
				if x["key"]:
					genes.add(x["name"])
					print x["name"]+":", len(clusters[cluster])-1, "read matches"
				else:
					reads.add(x["name"].rstrip("1").rstrip("2").rstrip("/"))
		else:
			print "No matches found"
			sys.exit()

	print len(reads), "reads in", len(genes), "genes"
	
	output=open("tmp.lst", "w")
	for read in reads:
		print >> output, read+"/1"
		print >> output, read+"/2"
		
	output.close()
	print "seqtk subseq -l 1000 ", options.ffastq, "tmp.lst", ">", ffastqname.split(".")[0]+"_subset.fastq"
	returnval=os.system(' '.join(["seqtk subseq -l 1000 ", options.ffastq, "tmp.lst", ">", ffastqname.split(".")[0]+"_subset.fastq"]))
	
	print "seqtk return value:", returnval
	if returnval!=0:
		print "Error: seqtk on forward reads failed."
		print "COMMAND: seqtk subseq -l 1000 ", options.ffastq, "tmp.lst", ">", ffastqname.split(".")[0]+"_subset.fastq"
		sys.exit()
	
	returnval=os.system(' '.join(["seqtk subseq -l 1000 ", options.rfastq, "tmp.lst", ">", rfastqname.split(".")[0]+"_subset.fastq"]))
	
	print "seqtk return value:", returnval
	if returnval!=0 :
		print "Error: seqtk on reverse reads failed."
		print "COMMAND: seqtk subseq -l 1000 ", options.rfastq, "tmp.lst", ">", rfastqname.split(".")[0]+"_subset.fastq"
		sys.exit()
	#options.db="test.fasta"
	returnval=os.system("smalt index -k 13 -s 1 index "+options.db)
	
	print "smalt index return value:", returnval
	if returnval!=0:
		print "Error: smalt index of reference failed."
		sys.exit()
	
	returnval=os.system(' '.join(["smalt map", "-d 0", "-y", idcutoff, "-o", "tmp.bam", "-f bam", "index", ffastqname.split(".")[0]+"_subset.fastq", rfastqname.split(".")[0]+"_subset.fastq"]))
	
	print "smalt map return value:", returnval
	
	if returnval!=0:
		print "Error: smalt mapping failed."
		sys.exit()
	
	returnval=os.system("samtools sort tmp.bam "+options.output)
	
	print "samtools sort return value:", returnval
	
	if returnval!=0:
		print "Error: samtools sorting failed."
		sys.exit()
	
#	returnval=os.system("~sh16/scripts/resistome/filter_doubleclipped_reads.py -b tmp_sort.bam -o "+options.output+"_filtered.bam")
#
#	print "filter_doubleclipped_reads.py return value:", returnval
#	
#	if returnval!=0:
#		print "Error: filtering of double-clipped reads failed."
#		sys.exit()
#	
	returnval=os.system("samtools index "+options.output+".bam")

	print "samtools index return value:", returnval
	
	if returnval!=0:
		print "Error: samtools indexing of bam failed."
		sys.exit()
	
	
	for proposed in genes:
		print "Analysing "+proposed
		
		returnval=os.system('grep -A 1 ">'+proposed+'" '+options.db+' > '+proposed+'_ref.fasta')
		
		returnval=os.system('~sh16/scripts/resistome/bam_filter.py -f pairedfastq -t contigsonemapped -c '+proposed+' -b '+options.output+'.bam -o tmp_proposed')
		returnval=os.system('smalt index index '+proposed+'_ref.fasta')
		if returnval!=0:
			continue
		returnval=os.system('smalt map -x -r 0 -o '+proposed+'.sam index tmp_proposed_1.fastq tmp_proposed_2.fastq')
		if returnval!=0:
			continue
		returnval=os.system('samtools view -S -b -o tmp.bam -h '+proposed+'.sam')
		if returnval!=0:
			continue
		returnval=os.system('samtools sort tmp.bam '+proposed)
		if returnval!=0:
			continue
		returnval=os.system('samtools index '+proposed+'.bam')
		if returnval!=0:
			continue
		returnval=os.system('~sh16/scripts/resistome/extract_clipping_info.py -b '+proposed+'.bam')
		if returnval!=0:
			continue
		returnval=os.system('~sh16/scripts/iCANDY.py '+proposed+'.bam starts.plot '+proposed+'_ref.fasta -d area -l 1 -g 0 -w -Y 0 -o '+proposed+'_coverage.pdf')
		if returnval!=0:
			continue
			
		
#		sys.exit()
		
		
		
		
		
#		returnval=os.system('sga preprocess --pe-mode 1 -o SGA.fastq tmp_proposed_1.fastq tmp_proposed_1.fastq')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga index -a ropebwt --no-reverse SGA.fastq')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga correct -k 41 --discard --learn -o SGA.ec.k41.fastq SGA.fastq')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga index -a ropebwt SGA.ec.k41.fastq')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga filter --homopolymer-check --low-complexity-check SGA.ec.k41.fastq')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga fm-merge -m 41 -o '+proposed+'.merged.fasta SGA.ec.k41.filter.pass.fa')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga index -d 1000000 '+proposed+'.merged.fasta')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga rmdup '+proposed+'.merged.fasta')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga overlap -m 41 '+proposed+'.merged.rmdup.fa')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga assemble -m 55 -g 0 -r 10 -o '+proposed+'.assemble.55 '+proposed+'.merged.rmdup.asqg.gz')
#		if returnval!=0:
#			continue
#		numcontigs=0
#		for line in open(proposed+'.merged.rmdup.fa',"rU"):
#			if len(line)>0 and line[0]==">":
#				numcontigs+=1
#		if numcontigs==0:
#			print "No contigs found for", proposed
#			continue
#		else:
#			print numcontigs, "found for", proposed
#			
#		returnval=os.system('smalt index SGA.merged.rmdup.index '+proposed+'.merged.rmdup.fa')
#		if returnval!=0:
#			continue
#		returnval=os.system('smalt map -x -r 0 -o '+proposed+'.sam SGA.merged.rmdup.index tmp_proposed_1.fastq tmp_proposed_1.fastq')
#		if returnval!=0:
#			continue
#		returnval=os.system('samtools view -S -b -o tmp.bam -h '+proposed+'.sam')
#		if returnval!=0:
#			continue
#		returnval=os.system('samtools sort tmp.bam '+proposed)
#		if returnval!=0:
#			continue
#		returnval=os.system('samtools index '+proposed+'.bam')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga-bam2de.pl -n 5 --prefix libPE '+proposed+'.bam')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga-astat.py -b '+str(numcontigs)+' -m 20 '+proposed+'.bam > libPE.astat')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga scaffold -m 50 --pe libPE.de -a libPE.astat -o scaffolds.n5.scaf '+proposed+'.merged.rmdup.fa')
#		if returnval!=0:
#			continue
#		returnval=os.system('sga scaffold2fasta -m 55 -a SGA.merged.rmdup.asqg.gz -o '+proposed+'scaffolds.n55.fa -d 1 --use-overlap --write-unplaced scaffolds.n5.scaf')
		

