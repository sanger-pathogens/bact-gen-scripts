#!/usr/bin/env python
import string, re, gzip
import os, sys, getopt, random, math
from Bio import SeqIO
from Bio.Seq import Seq
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_SeqIO import *
from Si_general import *

def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g'}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp


def Usage():
	print 'remove_blocks_from_aln.py Usage:'
	print 'Removes regions from a fasta alignment based on a tab file input. You can choose to remove or keep the regions in the tab file'
	print 'remove_blocks_from_aln.py [options]'
	print 'Options:'
	print '-a <file name>\talignment file name'
	print '-k\t\tkeep regions in tab file (default is to remove them)'
	print '-o <file name>\toutput file name'
	print '-t <file name>\ttab file name (containing regions to keep/remove)'
	print '-r <name>\treference name (optional, but required if there are gaps in the reference sequence relative to the tab file)'
	print '-R\t\tDo not remove blocks from reference sequence (default is to remove from all sequences)'
	print '-s <char>\tSymbol to use for removed regions (default = N)'
	print '-h\t\tshow this help'
	print 'Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2010'



##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "ho:a:t:kr:Rs:c", ["align=", "out=", "=tab", "keep", "reference=", "refrem", "symbol=", "cut"])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	outfile=''
	alnfile=''
	tabfile=''
	keepremove='r'
	reference=''
	refrem=True
	symbol="N"

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-a", "--align"):
			alnfile=arg
		elif opt in ("-t", "--tab"):
			tabfile=arg
		elif opt in ("-k", "--keep"):
			keepremove='k'
		elif opt in ("-c", "--cut"):
			keepremove='c'
		elif opt in ("-r", "--reference"):
			reference=arg
		elif opt in ("-R", "--refrem"):
			refrem=False
		elif opt in ("-s", "--symbol"):
			symbol=arg

	
	if alnfile=='' or not os.path.isfile(alnfile):
		print 'Error: Alignment file not found!'
		Usage()
		sys.exit()
	elif tabfile=='' or not os.path.isfile(tabfile):
		print 'Error: Tab file not found!'
		Usage()
		sys.exit()
	elif outfile=='':
		print 'Error: No output file specified!'
		Usage()
		sys.exit()
	elif keepremove not in ['k','r', 'c']:
		print "What the???"
		sys.exit()
	symbol=symbol.upper()
	if symbol not in ["N", "X", "?", "-"]:
		print 'Error: Symbol must be N, X, ? or -!'
		Usage()
		sys.exit()
		
	
	
	
	return alnfile, outfile, tabfile, keepremove, reference, refrem, symbol




if __name__ == "__main__":
	argv=sys.argv[1:]
	
	alnfile, outfile, tabfile, keepremove, reference, refrem, symbol=getOptions(argv)
	
	if keepremove=='c':
		refrem=True
	
	regions=[]
	
	for line in open(tabfile, 'rU'):
		try:
			if len(line.split())>2 and line.split()[0]=="FT" and (line.split()[1].lower() in ["misc_feature", "cds", "mobile_element", "fasta_record"]):
				if len(line.split()[2].split('..'))==1:
					start=int(line.split()[2])
					end=start
					regions.append([start-1,end-1,"f"])
				elif line.split()[2][:10]=="complement":
					line=line.replace(")","").replace("complement(", "")
					
					if int(line.split()[2].split('..')[0])<int(line.split()[2].split('..')[1]):
						start=int(line.split()[2].split('..')[0])
						end=int(line.split()[2].split('..')[1])
					else:
						start=int(line.split()[2].split('..')[1])
						end=int(line.split()[2].split('..')[2])
					
					regions.append([start-1,end-1,"r"])
				else:
					if int(line.split()[2].split('..')[0])<int(line.split()[2].split('..')[1]):
						start=int(line.split()[2].split('..')[0])
						end=int(line.split()[2].split('..')[1])
					else:
						start=int(line.split()[2].split('..')[1])
						end=int(line.split()[2].split('..')[0])
					regions.append([start-1,end-1,"f"])
		except StandardError:
			print line.split()
			sys.exit()
	
	print "Found", len(regions), "regions"
	
	
	regions.sort()
	
	
	sequences={}
	currseq=''
	foundref=False
	if reference!="":
		lines=[]
		count=-1
		curseqlist=[]
		append=curseqlist.append
		for linea in open(alnfile, "rU"):
			linea=linea.strip()
			if linea[0]==">":
				if count>-1 and foundref:
					sequence=''.join(curseqlist)
					break
				count=count+1
				#lines.append(linea.split()[0][1:]+'\n')
				name=linea.split()[0][1:]
				
				if name==reference:
					curseqlist=[]
					append=curseqlist.append
					foundref=True
				
			elif foundref:	
	#				lines[count]=lines[count]+linea
				append(linea)
		if count>-1 and foundref:
			sequence=''.join(curseqlist)
		curseqlist=[]
		if not foundref:
			DoError("Cannot find reference in alignment")
		else:
			reflen=len(sequence)
			print "Found reference:", reference, "("+str(reflen)+"bp)"
	
		print "Adjusting region locations to alignment positions"
		reftoaln={}
		refnum=0
		for alnnum, base in enumerate(sequence):
			if base!="-":# and base!="N":
				reftoaln[refnum]=alnnum
				refnum+=1
		for region in regions:
			if region[0] not in reftoaln or region[1] not in reftoaln:
				DoError("One of your regions has locations outside of the reference sequence length")
			region[0]=reftoaln[region[0]]
			region[1]=reftoaln[region[1]]
			if region[0]>region[1]:
				DoError("coordinated are the wrong way around")
	regions.sort()
	reftoaln=''
	#sys.exit()
	
	def remove_blocks_from_sequence(name, sequence, regions):
		if not refrem and name==reference:
			return name, sequence
		newsequence=''
		lastregionend=0
		distsincelastgoodblock=1
		for x, region in enumerate(regions):
			if region[1]<=lastregionend:
				distsincelastgoodblock+=1
				continue
			if region[0]<=lastregionend:
				region[0]=lastregionend
			
			if keepremove=='k':
				if region[2]=="f":
					newsequence=newsequence+sequence[region[0]:region[1]+1]
				else:
					newsequence=newsequence+revcomp(sequence[region[0]:region[1]+1])
			elif keepremove=='c':
				if x==0:
					newsequence=newsequence+sequence[:region[0]]
				else:
					newsequence=newsequence+sequence[regions[x-distsincelastgoodblock][1]+1:region[0]]
				if x==len(regions)-1:
					newsequence=newsequence+sequence[region[1]+1:]
			else:
				if x==0:
					newsequence=newsequence+sequence[:region[0]]+symbol*(region[1]+1-region[0])
				else:
					newsequence=newsequence+sequence[regions[x-distsincelastgoodblock][1]+1:region[0]]+symbol*(region[1]+1-region[0])
				if x==len(regions)-1:
					newsequence=newsequence+sequence[region[1]+1:]
#			if len(newsequence)!=region[1]+1:
#				print region, len(newsequence), lastregionend
#				sys.exit()
#			else:
#				print region, len(newsequence), "OK"
			lastregionend=region[1]
			distsincelastgoodblock=1
		return name, newsequence
	
	lines=[]
	count=-1
	curseqlist=[]
	append=curseqlist.append
	alnout=open(outfile,"w")
	names=set([])
	for linea in open(alnfile, "rU"):
		linea=linea.strip()
		if linea[0]==">":
			if count>-1:
				sequence=''.join(curseqlist)
				if count==0:
					reflen=len(sequence)
				else:
					if len(sequence)!=reflen:
						DoError("Input sequences of different lengths. Is your file an alignment?")
				newname, newseq=remove_blocks_from_sequence(name, sequence, regions)
				if len(newseq)!=len(sequence) and keepremove!='c' and keepremove!="k":
					DoError("Output and input sequences of different lengths. Do you have overlapping features in your inputfile?")
				print >> alnout, ">"+newname
				print >> alnout, newseq
			count=count+1
			#lines.append(linea.split()[0][1:]+'\n')
			curseqlist=[]
			append=curseqlist.append
			name=linea.split()[0][1:]
			if name in names:
				DoError("Found duplicate sequence name: "+name)
			else:
				names.add(name)
		else:	
#				lines[count]=lines[count]+linea
			append(linea)
	if count>-1:
		sequence=''.join(curseqlist)
		if len(sequence)!=reflen:
			DoError("Input sequences of different lengths. Your file is not an alignment")
		newname, newseq=remove_blocks_from_sequence(name, sequence, regions)
		if len(newseq)!=len(sequence) and keepremove!='c':
			DoError("Output and input sequences of different lengths. Do you have overlapping features in your inputfile?")
		print >> alnout, ">"+newname
		print >> alnout, newseq
	curseqlist=[]
	alnout.close()
	print "Adjusted", len(names), "sequences"
			
	
	
#	if len(newsequences[sequence])!=len(sequences[sequence]) and keepremove!='c':
#		print "ERROR: Output and input sequences of different lengths. Do you have overlapping features in your inputfile?"
#	
	print "Original alignment length:", len(sequence), "New alignment length:", len(newseq)
	print "Done."
