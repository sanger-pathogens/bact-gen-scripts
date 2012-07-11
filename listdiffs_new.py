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

#if a reference sequence has been given, check it is in the alignment

if ref not in sequs.keys():
	print "Cannot find reference sequence "+ref+" in alignment file"

#Check all sequences are the same length

seqlen=len(sequs[currseq])


for key in sequs.keys():
	if len(sequs[key])!=seqlen:
		print "Sequences of different lengths! Exciting..."
		sys.exit()


#open the orfs file if present

frame={}
if orfs!='':
	orflines=open(orfs,"rU").readlines()
	for line in orflines:
		words=line.strip().split(',')
		count=2
		for x in range(int(words[1]), int(words[2])+1):
			if x+count<int(words[2]):
				frame[x-1]=[words[0],words[3]]
			count=count-1
			if count==-1:
				count=2
	

#Find sites that contain a gaps

nonconstsites=[]
indel=0

geneticcode={'TTT':'Phe', 'TTC':'Phe', 'TTA':'Leu', 'TTG':'Leu', 'TCT': 'Ser', 'TCC': 'Ser','TCA': 'Ser','TCG': 'Ser', 'TAT': 'Tyr','TAC': 'Tyr', 'TAA': 'STOP', 'TAG': 'STOP', 'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp', 'CTT': 'Leu','CTC': 'Leu','CTA': 'Leu','CTG': 'Leu', 'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'CAT': 'His', 'CAC': 'His', 'CAA': 'Gin', 'CAG': 'Gin', 'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met', 'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val', 'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu', 'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'}
aas={}
for key in sequs.keys():
	aas[key]=[]

refpos=0
count=-1
comp={'A':'T', 'T':'A', 'C':'G', 'G':'C'}
lastorf=''

for i in range(0,seqlen):
	states={}
	gapsite='n'
	newindel='n'
	numgaps=0
	
	if orfs!='':
		if frame.has_key(refpos):
			if frame[refpos][0]!=lastorf:
				count=0
				lastorf=frame[refpos][0]
			else:
				count=count+1
		else:
			count=-1

	if ref!='' and sequs[ref][i] not in ['-','X','?']:
		refpos=refpos+1


	for key in sequs.keys():

		if sequs[key][i]=='-':
			numgaps=numgaps+1
			gapsite='y'
			if i==0 or sequs[key][i-1]!='-':
				newindel='y'
		elif sequs[key][i] not in ['-','X','?'] and states.has_key(sequs[key][i]):
                        states[sequs[key][i]]=states[sequs[key][i]]
                elif sequs[key][i] not in ['-','X','?']:
                        states[sequs[key][i]]=1
		
		if orfs!='':
			if seqlen>i+2 and count==0:
				if frame[refpos-1][1]=='+':
					codons=sequs[key][i]+sequs[key][i+1]+sequs[key][i+2].upper()
				elif frame[refpos-1][1]=='-':
					codons=comp[sequs[key][i+2]]+comp[sequs[key][i+1]]+comp[sequs[key][i]].upper()
				else:
					print "Error, 4th column of orf file must indicate direction (- or +). Found "+frame[refpos-1][1]
				for j in range(0,3):
					if geneticcode.has_key(codons):
						aas[key].append(geneticcode[codons])
					else:
						aas[key].append('-')
			elif count==-1:
				aas[key].append('-')

	nostates=len(states.keys())
			
	if newindel=='y':
		indel=indel+1
	if gapsite=='y':
		nonconstsites.append([i,indel,nostates, numgaps, refpos])
	elif gapsite=='n' and nostates>1:
		nonconstsites.append([i,-1,nostates, numgaps, refpos])

	if count==2:
		count=-1


#Print sites with gaps or SNPs to output file

output=open(outfile, 'w')

print >> output, 'Position in alignment',
if ref!='':
	print >> output, '\tPosition in '+ref,
if orfs!='':
	print >> output, '\tORF',
for i in seqorder:
	print >> output, '\t'+i,
if orfs!='':
	print >> output, '\tAmino acids',

print >> output, '\n',

lastindel=-2
length=0

for i in nonconstsites:

	donepos='n'
	currentaas=[]
	for key in seqorder:
		
		if orfs!='':
			if aas[key][i[0]] not in currentaas:
				currentaas.append(aas[key][i[0]])

		if sequs[key][i[0]]=='-' and (i[0]==0 or sequs[key][i[0]-1]!='-'):
			
			if donepos=='n':
				print >> output, i[0]+1,
				if refpos!='':
					if sequs[ref][i[0]]=='-':
						print >> output, '\t-',
					else:
						print >> output, '\t'+str(i[4]),
				if orfs!='':
					if frame.has_key(i[4]):
						print >> output, '\t'+frame[i[4]][0],
					else:
						print >> output, '\t',
				donepos='y'

			j=i[0]
			length=0
			while j<len(sequs[key]) and sequs[key][j]=='-':
				length=length+1
				j=j+1
			if i[3]>(len(seqorder)/2):
				print >> output, '\t-',
			else:
				print >> output, '\t'+str(length)+' bp deletion',
		
		elif i[3]>(len(seqorder)/2) and (i[1]!=lastindel) and sequs[key][i[0]]!='-':
			if donepos=='n':
				print >> output, i[0]+1,
				if ref!='':
					if sequs[ref][i[0]]=='-':
						print >> output, '\t-',
					else:
						print >> output, '\t'+str(i[4]),
				if orfs!='':
					if frame.has_key(i[4]):
						print >> output, '\t'+frame[i[4]][0],
					else:
						print >> output, '\t',
				donepos='y'

			print >> output, '\t'+str(length)+' bp insertion '+sequs[key][i[0]], 


		elif i[2]>1 or (i[1]!=lastindel and i[1]!=-1):
			if donepos=='n':
                                print >> output, i[0]+1,
				if ref!='':
					if refpos!='':
						if sequs[ref][i[0]]=='-':
							print >> output, '\t-',
						else:
							print >> output, '\t'+str(i[4]),
				if orfs!='':
					if frame.has_key(i[4]):
						print >> output, '\t'+frame[i[4]][0],
					else:
						print >> output, '\t',
                                donepos='y'

			print >> output, '\t'+sequs[key][i[0]],
	
	if '-' in currentaas:
		currentaas.remove('-')
	if donepos=='y':
		
		if orfs!='' and frame!=0:
			print >> output,'\t', ', '.join(currentaas)
		else:
			print >> output, '\n',
	lastindel=i[1]
output.close()

print "Done."
