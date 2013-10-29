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
	parser.add_option("-f", "--forwardfastq", action="store", dest="ffastq", help="Input forward fastq files", default="")
	parser.add_option("-r", "--reversefastq", action="store", dest="rfastq", help="Input reverse fastq files", default="")
	parser.add_option("-d", "--database", action="store", dest="db", help="Database of files to search for", default="")
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
	
	if options.ffastq=="":
		DoError("Forward fastq file required")
	elif not os.path.isfile(options.ffastq):
		DoError("Cannot find forward fastq file")
	
	if options.rfastq=="":
		print "No reverse fastq file selected. Assuming data is single ended"
		DoError("Single ended data not supported yet")
	elif not os.path.isfile(options.rfastq):
		DoError("Cannot find reverse fastq file")

	if options.db=="":
		DoError("Database file (-d) required")
	elif not os.path.isfile(options.db):
		DoError("Cannot find file "+options.db)
	

	fastaname=options.ffastq.split(".")[0]+"_resistome.fasta"
	
	print "Extracting fastq files to fasta"

	#returnval=os.system(' '.join(["/nfs/users/nfs_s/sh16/scripts/resistome/fastool --to-fasta", options.ffastq, options.rfastq, ">", fastaname]))

	returnval=os.system(' '.join(["seqtk seq", "-l 0",  "-A", options.ffastq, ">", fastaname]))

	if returnval!=0 or not os.path.isfile(fastaname):
		print "Error: Conversion failed. Please check the format of your forward fasttq files."
		sys.exit()

	returnval=os.system(' '.join(["seqtk seq", "-l 0", "-A", options.rfastq, ">>", fastaname]))

	if returnval!=0:
		print "Error: Conversion failed. Please check the format of your reverse fasttq files."
		sys.exit()
	

	print "Calculating read length (assuming all reads are equal)"
	readlength=0
	for x, line in enumerate(open(fastaname, "rU")):
		if x>1:
			break
		elif x==1:
			readlength=len(line.strip())

	if readlength<=20:
		DoError("Read length must be greater than 20")

	print "Read length is", readlength
	idcutoff=str(options.minmatch/float(readlength))
	print "CD-Hit matches will be retained when greater than", idcutoff

	print "Clustering reads to genes with CD-Hit-est-2d..."
	
	returnval=os.system(' '.join(["cd-hit-est-2d -i", options.db, "-i2", fastaname, "-aS", idcutoff, "-G 0 -M 0",  "-o", options.db+"_clusters -d 10000"]))
	
	if returnval!=0 or not os.path.isfile(options.db+"_clusters.clstr"):
		print "Error: Clustering failed. Please check the format of your sequence database file."
		sys.exit()
	
	clusters=parse_CDHitEST_clstr_file(open(options.db+"_clusters.clstr", "rU"))
	
	reads=set([])
	genes=set([])
	
	for cluster in clusters:
		if len(clusters[cluster])>1:
			for x in clusters[cluster]:
				if x["key"]:
					genes.add(x["name"])
					print x["name"]+":", len(clusters[cluster])-1, "read matches"
				else:
					reads.add(x["name"].strip("/1").strip("/2"))

	print len(reads), "reads in", len(genes), "genes"

	output=open("tmp.lst", "w")
	for read in reads:
		print >> output, read+"/1"
		print >> output, read+"/2"
		
	output.close()
	
	returnval=os.system(' '.join(["seqtk subseq -l 1000 ", options.ffastq, "tmp.lst", ">", options.ffastq.split(".")[0]+"_subset.fastq"]))
	
	print "seqtk return value:", returnval
	
	returnval=os.system(' '.join(["seqtk subseq -l 1000 ", options.rfastq, "tmp.lst", ">", options.rfastq.split(".")[0]+"_subset.fastq"]))
	
	print "seqtk return value:", returnval
	options.db="test.txt"
	returnval=os.system("smalt index -k 13 -s 1 index "+options.db)
	
	print "smalt index return value:", returnval
	
	returnval=os.system(' '.join(["smalt map", "-d 0", "-o", "tmp.bam", "-f bam", "index", options.ffastq.split(".")[0]+"_subset.fastq", options.rfastq.split(".")[0]+"_subset.fastq"]))
	
	print "smalt map return value:", returnval
	
	returnval=os.system("samtools sort tmp.bam tmp_sort")
	
	print "samtools sort return value:", returnval
	
	returnval=os.system("samtools index tmp_sort.bam")

	print "samtools index return value:", returnval
	
	
	
