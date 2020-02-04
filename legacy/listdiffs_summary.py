#!/usr/bin/env python

#Extracts orthologs from tab files and 



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, commands, getopt, gzip


def Usage():
	print 'listdiffs.py Usage:'
	print 'listdiffs.py -i [input alignment] -o [output file name] -f [frame 123, default 0 = unknown] -r [name of sequence to use as reference] -a [do alignment] {-h}'
	print 'or'
	print 'listdiffs.py --in [input alignment] --out [output file name] --frame [frame 123, default 0 = unknown] --reference [name of sequence to use as reference] --align [do alignment] --help'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'



#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------

def getOptions(arg):


	try:
		opts, args = getopt.getopt(argv, "hi:o:f:ar:", ["help", "in=", "out=", "orfs=", "align", "reference="])
	except getopt.GetoptError:
		Usage()
		sys.exit(2)

	inputfile=''
	outfile=''
	orfs=''
	align='n'
	ref=""

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-i", "--in"):
			inputfile=arg
                elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-f", "--orfs"):
			orfs=arg
		elif opt in ("-a", "--align"):
			align='y'
		elif opt in ("-r", "--reference"):
			ref=arg

	if inputfile=='':
		print 'No input fasta file selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=inputfile+"_diversity.txt"
	if orfs!='' and ref=='':
		print 'Reference sequence required if list of orfs given'
		Usage()
		sys.exit()

	return inputfile, outfile, orfs, align, ref




if __name__ == "__main__":
	argv=sys.argv[1:]
	inputfile, outfile, orfs, align, ref=getOptions(argv)

print "Summarising diversity in "+inputfile

#Open the input file and align it if necessary

if align=='y':
	os.system("muscle -in "+inputfile+" -out "+inputfile+".aln")
	lines=open(inputfile+".aln", "rU").readlines()
else:
	lines=open(inputfile, "rU").readlines()


sequs={}
currseq=''
seqorder=[]

#Read in each sequence

for line in lines:
	if len(line)>0 and line[0]==">":
		sequs[line.strip().split()[0][1:]]=''
		currseq=line.strip().split()[0][1:]
		seqorder.append(line.strip().split()[0][1:])
	elif len(line)>0:
		sequs[currseq]=sequs[currseq]+line.strip().upper()

output= open(outfile,'w')
print >> output, "Sequence1\tSequence2\tNumber of SNPs\tNumber of indels\tTotal Indel length\tNumber of indels\tTotal Indel length"



for seq in sequs.keys():
	print >> output, seq,
	for seq2 in sequs.keys():
		identity=0
		SNPs=0
		inindel=['n','n']
		indelcount=[0,0]
		indellen=[0,0]
		if seq!=seq2:
			print >> output, "\t"+seq2,
			for base in range(len(sequs[seq])):
				if sequs[seq][base].upper() in ['A', 'C', 'G', 'T']:
					inindel[0]='n'
					if sequs[seq2][base].upper() in ['A', 'C', 'G', 'T']:
						inindel[1]='n'
						if sequs[seq][base].upper()==sequs[seq2][base].upper():
							identity=identity+1
						else:
							SNPs=SNPs+1
					else:
						indellen[1]=indellen[1]+1
						if inindel[1]=='n':
							indelcount[1]=indelcount[1]+1
							inindel[1]='y'
				elif sequs[seq2][base].upper() in ['A', 'C', 'G', 'T']:
					indellen[0]=indellen[0]+1
					if inindel[0]=='n':
						indelcount[0]=indelcount[0]+1
						inindel[0]='y'
			print >> output, "\t"+str(SNPs)+"\t"+str(indelcount[0])+"\t"+str(indellen[0])+"\t"+str(indelcount[1])+"\t"+str(indellen[1])





output.close()



print "Done."










