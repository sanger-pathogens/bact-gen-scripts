#!/usr/bin/env python
#Converts tfa fiels to fasta, runs Gblocks, converts to phylip format and runs phymlboot.py for a list of KOGs
#./runtrees.py [inputfile]
#Where inputfile is a text file containing a list of KOGs in the directory /research/gene_hunter/mitochondrial1/

import string, re
import os, sys, commands, getopt

def Usage():
	print '\nchangetype.py, Written by Simon Harris, Newcastle University, 2007\n'
	print 'Changes format of alignments using readseq\n'
	print 'Usage:\n'
	print 'changetype.py [-i] {input} [-t] {output file type} [-h]'
	print 'or'
	print 'changetype.py [--input] {input} [--type] {output file type} [--help]'
	print '\n\t-i\tinput file name or directory\n\t-t\tOutput file type:\n\t\t1 = IG/Stanford\n\t\t2 = Genbank/GB'
	print '\t\t3 = NBFR\n\t\t4 = EMBL\n\t\t5 = GCG\n\t\t6 = DNAStrider\n\t\t7 = Fitch\n\t\t8 = Pearson/Fasta'
	print '\t\t9 = Zucker\n\t\t10 = Olsen (in only)\n\t\t11 = Phylip3.2 (sequencial)\n\t\t12 = Phylip (interleaved)\n\t\t13 = Plain/Raw\n\t\t14 = PIR/CODATA'
	print '\t\t15 = MSF\n\t\t16 = ASN.1\n\t\t17 = NEXUS (PAUP)\n\t\t18 = Pretty (out only)\n\t-h\thelp\n'
	print 'Copyright Simon R Harris, Newcastle University, Newcastle Upon Tyne, UK. 2007\n'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def main(arg):


	try:
		opts, args = getopt.getopt(argv, "i:t:h", ["input=", "type=", "help"])
	except getopt.GetoptError:
		Usage()
		sys.exit(2)
	
	typeno=''
	inputfile=''

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-i", "--input"):
			inputfile=arg
		elif opt in ("-t", "--type"):
			typeno=arg
	
	if inputfile=='':
		print 'no input file/directory selected'
		Usage()
		sys.exit()
		
	if typeno=='':
		print 'no output file type selected'
		Usage()
		sys.exit()

	return inputfile, typeno

#------------------------------------------------------------------------------------os.system(line)
# Run main() with command line arguments
#------------------------------------------------------------------------------------
if __name__ == "__main__":
	argv=sys.argv[1:]
	directory, typeno=main(argv)
	
#------------------------------------------------------------------------------------
# Check if input is a file or directory
#------------------------------------------------------------------------------------

if os.path.isdir(directory)==1:
	filenames=os.listdir(directory)
elif os.path.isfile(directory)==1:
	filenames=[]
	filenames.append(directory)
else:
	print 'invalid input file/directory'
	sys.exit()
	
#------------------------------------------------------------------------------------
# Call readseq to convert file to desired format
#------------------------------------------------------------------------------------
	
for x in range(0,len(filenames)):
	if filenames[x][0]!='.':
		print 'Converting '+filenames[x]
		textcall = ("").join(['java -cp /nfs/users/nfs_s/sh16/scripts/readseq.jar run -a -f '+typeno+' '+filenames[x] ])
		os.system(textcall)
		
		#If converting to type 11, edit output file to remove I on first line
		
		if typeno=='11':
			changedfile=open(filenames[x]+'.phylip2','rU')
			lines=changedfile.read()
			changedfile.close()
			if lines.split()[2]=='I':
				templines=lines.split(' ')[:3]+lines.split(' ')[4:]
				lines=' '.join(templines)
				changedfile=open(filenames[x]+'.phylip2','w')
				print >> changedfile, lines



