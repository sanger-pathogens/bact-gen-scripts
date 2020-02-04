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

	parser.add_option("-d", "--database", action="store", dest="db", help="Database of files to search for", default="")
	#Could add more options in here so people can specify similarities etc.
	
	
	return parser.parse_args()


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
				if len(line.split()[-1].split('/'))!=2:
				       print "Cannot read strand/percentid"
				       print line
				       sys.exit()
				strand=line.split()[-1].split('/')[0]
				try:
				       percentid=float(line.split()[-1].split('/')[1].replace("%",""))
				except ValueError:
				       print "Cannot read strand/percentid"
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
	
	if options.db=="":
		DoError("Database file (-d) required")
	
	if not os.path.isfile(options.db):
		DoError("Cannot find file "+options.db)
	
	
	print "Clustering genes with CD-Hit est..."
	
	returnval=os.system("cd-hit-est -aS 0.9 -i "+options.db+" -o "+options.db+"_clusters -d 10000")
	
	if returnval!=0 or not os.path.isfile(options.db+"_clusters.clstr"):
		print "Error: Clustering failed. Please check the format of your sequence database file."
		sys.exit()
	
#	clusters=parse_CDHitEST_clstr_file(open(options.db+"_clusters.clstr", "rU"))
#
#	print clusters
#	
#	print "Creating cluster fasta directory..."
#	
#	os.system("mkdir resistome_cluster_directory")
	
	
	
	
	
	
