#!/usr/bin/env python
#/usr/bin/python


#Compares mummer snp output files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re, gzip
import os, sys, getopt, random, math
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.5/site-packages/']))
from scipy.stats import chi2


#########
# Usage #
#########

def Usage():
	os.system('clear')
	print '\nsummarise_snps.py Usage:\n'
	print 'nsummarise_snps.py [OPTIONS] <list of fastq files>'
	print "\nINPUT OPTIONS:"
	print "-r Reference dna sequence:\t<fasta file> (required)"
	print "-e Reference embl file:\t\t<embl or tab file>"
	print "-s reads are multiplexed with stupid naming convention e.g. 2610_4_1_10 and 2610_4_2_10 for pairs of strain 10 in pool 4"
	print "-p reads are paired-end"
	print "\nSNP calling OPTIONS:"
	print "-q Minimum mapping quality:\t<integer between 1 and 60> [default=30]"	
	print "-d Minimum mapping depth:\t<integer between 1 and 100,000> [default=1]"
	print "\nOUTPUT OPTIONS:"
	print "-o Prefix for output files:\t<file name prefix>"
	print "-a Create SNP alignment file"
	print "-g Create coverage plots"
	print "-t Create SNP tab file"
	print "-P Run phylogeny with RAxML"
	print "-M Model of evolution\t\t<GTR, GTRGAMMA, etc.> [default=GTRGAMMA]"
	print "-f Use fast bootstrap [recommended for datasets with large numbers of strains or snps]"
	print "-b No. bootstrap replicates\t<value between 0 and 10,000> [default=100]"
	print "\nUSAGE OPTIONS:"
	print "-h Show this help"
	print "-I Run interactively using a phylip-style menu"
	print '\nCopyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009\n'

	

##############################
# Get command line arguments #
##############################


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "hd:e:r:o:tai:Rcpgb:M:q:h:svIm:Pn:xfT:", ["depth=", "model=", "bootstrap=", "help", "phylogeny", "embl=", "ref=", "out=", "align", "chisquared=", "tabfile", "graphs", "quality=", "interactive", "pairedend", "rtype="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	ref=''
	outfile=''
	inputdirs=[]
	i=0
	tabfile='n'
	align='n'
	embl=''
	graphs='n'
	raxml='n'
	bootstrap=100
	model='GTRGAMMA'
	quality=30
	depth=1
	interactive='n'
	pairedend='n'
	chisquared='n'
	fastboot='n'
	rtype='solexa'
	multiplexnamed='n'
	heteropercent=75

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-t", "--tabfile"):
			tabfile='y'
		elif opt in ("-f", "--fastboot"):
			fastboot='y'
		elif opt in ("-P", "--phylogeny"):
			raxml='y'
		elif opt in ("-a", "--align"):
			align='y'
		elif opt in ("-g", "--graphs"):
			graphs='y'
		elif opt in ("-e", "--embl"):
			embl=arg
		elif opt in ("-b", "--bootstrap"):
			bootstrap=int(arg)
		elif opt in ("-M", "--model"):
			model=arg
		elif opt in ("-s", "--stupid"):
			multiplexnamed='y'
		elif opt in ("-q", "--quality"):
			quality=int(arg)
		elif opt in ("-d", "--depth"):
			depth=int(arg)
		elif opt in ("-I", "--interactive"):
			interactive='y'
		elif opt in ("-p", "--pairedend"):
			pairedend='y'
		elif opt in ("-x", "--chisquared"):
			chisquared='y'
		elif opt in ("-T", "--rtype"):
			rtype=int(arg)

	inputdirs=args
	
	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=ref.split('/')[-1].split('.')[0]+"_q"+str(quality)+"_d"+str(depth)
	if align=='n' and raxml=='y':
		align='y'
	
	
	if interactive=='y':
		ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed=menusystem(ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed)

	if ref=='':
		print 'Error: No reference dna file selected!'
		Usage()
		sys.exit()
	
	
	return ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed



####################
# Interactive Menu #
####################

def menusystem(ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed):
	
	os.system('clear')
	
	if outfile==ref.split('/')[-1].split('.')[0]+"_q"+str(quality)+"_d"+str(depth):
		outorig='y'
	else:
		outorig='n'
	
	print "\nsummarise_snps.py: Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009"
	
	print "\nINPUT OPTIONS:"
	
	if ref=='':
		print "r: Reference dna sequence:\t\tNone selected (required)"
	else:
		print "r: Reference dna sequence:\t\t"+ref
		if embl=='':
			print "e: Reference embl file:\t\t\tNone selected"
		else:
			print "e: Reference embl file:\t\t\t"+embl
		print "\nOPTIONS:"
		if pairedend=='n':
			print "p: Use paired-end reads:\t\tno"
		else:
			print "p: Use paired-end reads:\t\tyes"
		print "q: Minimum snp quality:\t\t\t"+str(quality)
		print "d: Minimum mapping depth:\t\t"+str(depth)
		print "D: Reset to default mapping parameters"
		print "\nOUTPUT OPTIONS:"
		print "o: Prefix for output files:\t\t"+outfile
		if align=='n':
			print "a: Create SNP alignment file:\t\tno"
		else:
			print "a: Create SNP alignment file:\t\tyes"
		if graphs=='n':
			print "g: Create coverage plots:\t\tno"
		else:
			print "g: Create coverage plots:\t\tyes"
		if tabfile=='n':
			print "t: Create SNP tab file:\t\t\tno"
		else:
			print "t: Create SNP tab file:\t\t\tyes"
		if chisquared=='n':
			print "x: Calculate chi-squared:\t\tno"
		else:
			print "x: Calculate chi-squared:\t\tyes"
		if raxml=='n':
			print "P: Run phylogeny with RAxML:\t\tno"
		else:
			print "P: Run phylogeny with RAxML:\t\tyes"
		if raxml=='y':
			print "M:  Model of evolution:\t\t\t"+model
			if fastboot=='n':
				print "f:  Use fast bootstrap:\t\t\tno"
			else:
				print "f:  Use fast bootstrap:\t\t\tyes"
			if bootstrap==0:
				print "b:  Run bootstrap:\t\t\tno"
			else:
				print "b:  No. bootstrap replicates:\t\t"+str(bootstrap)
		
	print "\nQ: QUIT"
	
	if ref=="":
		message="\nPlease select an option:"
		inputlist=['r', 'Q']
	else:
		message="\nPlease select an option or type y to run:"
		inputlist=['r', 'e', 'R', 'q','d','x','v','o','a','g','t','P','y','Q','D','p']
		if pairedend=='y':
			inputlist=inputlist+['i']
		if raxml=='y':
			inputlist=inputlist+['M','b','f']
			
	ui=''
	while ui not in inputlist:
		ui=raw_input(message+' ')

	if ui=='y':
		os.system('clear')
		return ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed
		
	elif ui=='r':
		oldref=ref
		ref=''
		while not os.path.isfile(ref):
			ref=raw_input('Enter reference file name or Q to go back to the menu: ')
			if ref=='Q':
				ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed=menusystem(oldref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed)
			elif not os.path.isfile(ref):
				print "File not found"
				
	elif ui=='e':
		oldembl=embl
		embl=''
		while not os.path.isfile(embl):
			embl=raw_input('Enter embl file name or Q to go back to the menu: ')
			if embl=='Q':
				ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed=menusystem(ref, inputdirs, outfile, tabfile, align, oldembl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed)
			elif not os.path.isfile(embl):
				print "File not found"

	
	elif ui=='q':
		quality=0
		while quality > 60 or quality < 1:
			quality=int(raw_input('Enter minimum snp quality (1-60): '))
	
	elif ui=='d':
		depth=0
		while depth > 100000 or depth < 1:
			depth=int(raw_input('Enter minimum read depth for mapping and SNP calling (1-100,000): '))
	
	elif ui=='o':
		outfile=''
		while outfile=='':
			outfile=raw_input('Enter prefix for output file names, or D to use the default: ')
			if outfile!='D' and outfile!='':
				outorig='n'
			elif outfile=='D':
				outorig='y'
				
	elif ui=='a':
		if align=='n':
			align='y'
		else:
			align='n'
			
	elif ui=='g':
		if graphs=='n':
			graphs='y'
		else:
			graphs='n'
			
	elif ui=='t':
		if tabfile=='n':
			tabfile='y'
		else:
			tabfile='n'
			
	elif ui=='x':
		if chisquared=='n':
			chisquared='y'
		else:
			chisquared='n'
			
	elif ui=='f':
		if fastboot=='n':
			fastboot='y'
		else:
			fastboot='n'
			
	elif ui=='P':
		if raxml=='n':
			raxml='y'
			align='y'
		else:
			raxml='n'
	
	elif ui=='M':
		if model=='GTRGAMMA':
			model='GTR'
		elif model=='GTR':
			model='GTRCAT'
		else:
			model='GTRGAMMA'
	
	elif ui=='b':
		bootstrap=-1
		while bootstrap > 10000 or bootstrap < 0:
			bootstrap=int(raw_input('Enter number of bootstrap replicates (0-10,000): '))
	
	elif ui=='D':
		quality=30
		depth=5
	
	elif ui=='p':
		if pairedend=='n':
			pairedend='y'
		else:
			pairedend='n'
			
	elif ui=='Q':
		sys.exit()
	
	if outorig=='y':
		outfile=ref.split('/')[-1].split('.')[0]+"_q"+str(quality)+"_d"+str(depth)
		
	ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed=menusystem(ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed)

	return ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, depth, pairedend, chisquared, fastboot, rtype, multiplexnamed


##########################################################
# Algorithm to count steps needed to move between codons #
##########################################################

def countcodonchanges(codon, SNPcodon, geneticcode, sd, nd, loopsd=0, loopnd=0, pathcount=0):
	
	for x in range(3):
		if codon[x]!=SNPcodon[x]:
			newSNPcodon=SNPcodon[:x]+codon[x]+SNPcodon[x+1:]
			
			#print  SNPcodon, newSNPcodon, geneticcode[SNPcodon], geneticcode[newSNPcodon]
			
			if geneticcode[newSNPcodon]=='*':
				continue
			elif geneticcode[SNPcodon]==geneticcode[newSNPcodon]:
				newloopnd=loopnd
				newloopsd=loopsd+1
			else:
				newloopnd=loopnd+1
				newloopsd=loopsd
			
			
			#print SNPcodon, newSNPcodon, codon, sd, nd, newloopsd, newloopnd, pathcount
			if newSNPcodon!=codon:
				sd, nd, pathcount=countcodonchanges(codon, newSNPcodon, geneticcode, sd, nd, newloopsd, newloopnd, pathcount)
			
			else:
				sd=sd+newloopsd
				nd=nd+newloopnd
				pathcount=pathcount+1
						
	return sd, nd, pathcount


############################################################################
# Algorithm to calculate dN/dS and related stats for two aligned sequences #
############################################################################

def dnbyds(CDS, SNPseq, CDSbasenumbers):
	
	#geneticcode={'TTT':'Phe', 'TTC':'Phe', 'TTA':'Leu', 'TTG':'Leu', 'TCT': 'Ser', 'TCC': 'Ser','TCA': 'Ser','TCG': 'Ser', 'TAT': 'Tyr','TAC': 'Tyr', 'TAA': 'STOP', 'TAG': 'STOP', 'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp', 'CTT': 'Leu','CTC': 'Leu','CTA': 'Leu','CTG': 'Leu', 'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'CAT': 'His', 'CAC': 'His', 'CAA': 'Gin', 'CAG': 'Gin', 'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met', 'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val', 'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu', 'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'}
	geneticcode={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
	
	codonsynonyms={}
	
	for codon in geneticcode.keys():
		thiscodon={}
		codonsynonyms[codon]=0.0
		
		for x in range(3):
			numsyn=0.0
			numnotstop=3
			for y in ['A', 'C', 'G', 'T']:
				if codon[x]!=y:
					newcodon=codon[:x]+y+codon[x+1:]
					if geneticcode[newcodon]==geneticcode[codon]:
						numsyn=numsyn+1
					elif geneticcode[newcodon]=='*':
						numnotstop=numnotstop-1
			codonsynonyms[codon]=codonsynonyms[codon]+(numsyn/numnotstop)
	
	S1=0.0
		
	for x in range(0,len(CDS),3):
		codon=CDS[x:x+3]
		if not '-' in codon:
			S1=S1+(float(codonsynonyms[codon]))
	
	N=(len(CDS))-S1
	
	#if N > 1000000:
	#	print "Number of synonymous sites in reference =", S1, ", Number of nonsynonymous sites in reference =", N,

	
	S=0.0
	N=0.0
	S1=0.0
	N1=0.0
	S2=0.0
	N2=0.0
	Sd=0.0
	Nd=0.0
	pS=0.0
	pN=0.0
	gapcount=0
	numcodons=0
	varianceS=0.0
	varianceN=0.0
	z=0.0
	dN=0.0
	dS=0.0
	SNPtype={}
	AAtype={}
	
	if len(CDS)!=len(SNPseq):
		print "Error: sequences must be the same length to calculate dN/dS!"
		return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}, SNPtype, AAtype
	
	for x in range(0,len(CDS),3):
		if x+3>len(CDS):
			break
		numcodons=numcodons+1
		codon=CDS[x:x+3]
		SNPcodon=SNPseq[x:x+3]
		codonposn=[CDSbasenumbers[x], CDSbasenumbers[x+1], CDSbasenumbers[x+2]]
		newSNPcodon='---'
		
		if '-' in codon:
			gapcount=gapcount+3
			continue
		
		if codon!=SNPcodon:
			
			for y,z in enumerate(codon):
			
				if SNPcodon[y]!=z and SNPcodon[y]!='-':
					newSNPcodon=codon[:y]+SNPcodon[y]+codon[y+1:]
			
					if geneticcode[newSNPcodon]=='*':
						SNPtype[codonposn[y]]='2'
					elif geneticcode[codon]=='*':
						SNPtype[codonposn[y]]='3'
					elif geneticcode[newSNPcodon] == geneticcode[codon]:
						SNPtype[codonposn[y]]='S'
					else:
						SNPtype[codonposn[y]]='N'
						if not '-' in SNPcodon:
							AAtype[codonposn[y]]=geneticcode[SNPcodon]
					
					#print codon, SNPcodon, newSNPcodon, codonposn[y], SNPtype[codonposn[y]]
					
				
				elif SNPcodon[y]=='-':
					SNPtype[codonposn[y]]='1'
					
			
			
		
		if '-' in codon or '-' in SNPcodon or geneticcode[codon]=='*' or geneticcode[SNPcodon]=='*':
			gapcount=gapcount+3
			continue

		
		#s=float(codonsynonyms[codon])/3
		#n=float(3-s)
		S1=S1+(float(codonsynonyms[codon]))
		S2=S2+(float(codonsynonyms[SNPcodon]))
		
		
		
		sd=0.0
		nd=0.0
		
		pathcount=0
		if codon!=SNPcodon:
			sd, nd, pathcount=countcodonchanges(codon, SNPcodon, geneticcode, sd, nd)
		
		if pathcount>0:
			sd=float(sd)/pathcount
			nd=float(nd)/pathcount
		
		
		Sd=Sd+sd
		Nd=Nd+nd
	
	S=(S1+S2)/2
	N=(len(CDS)-gapcount)-S
	
	#pNb=pN/numcodons
	#pSb=pS/numcodons
	
	if N!=0:
		pN=Nd/N
	if S!=0:
		pS=Sd/S
	
	
	if pS==0:
		#print "No sites are synonymous."
		return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}, SNPtype, AAtype
	if pS<0.75 and pN<0.75:
		dS=(-3*(math.log(1 - ((pS*4)/3))))/4
		dN=(-3*(math.log(1 - ((pN*4)/3))))/4
		if dN==-0.0:
			dN=0.0
		varianceS=(9 * pS * (1 - pS))/(((3 - (4 *pS)) **2) * (len(CDS)-gapcount));
		varianceN=(9 * pN * (1 - pN))/(((3 - (4 *pN)) **2) * (len(CDS)-gapcount));
		z=(dN - dS) / math.sqrt(varianceS + varianceN)
		
	else:
		#print "Too divergent for JC! Using pN/pS instead."
		dS=pS
		dN=pN

	return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}, SNPtype, AAtype


###############################################################
# Run a bootstrapped phyogeny of an alignment file with RAxML #
###############################################################

def RAxMLphylogeny(alignfile, model, bootstrap, fastboot, force='y'):
	outlen=0
	if force=='y':
		userinput='y'
	else:
		userinput='x'
	
	alignfileprefix=alignfile.split('.')[0]
	
	if '/' in outfile:
		outlen=-1*len(alignfileprefix)
	if os.path.isfile('RAxML_info.'+alignfileprefix):
		print '\nRAxML files with extension '+alignfileprefix+' already exist!'
		while userinput not in ['y','n']:
			userinput=raw_input('Overwrite? (y/n): ')
			userinput.lower()
		if userinput=='y':
			os.system('rm RAxML_*'+alignfileprefix)
	print "Running RAxML phylogeny of "+alignfile+" with the "+model+" model of evolution and "+str(bootstrap)+" bootstrap replicates..."
	if fastboot=='n':
		os.system("RAxML -f i -b "+str(random.randrange(1,99999))+" -# "+str(bootstrap)+" -m "+model+" -s "+alignfile+" -n "+alignfileprefix+"boot")
		os.system("RAxML -f d -m "+model+" -s "+alignfile+" -n "+alignfileprefix+"ml")
		os.system("RAxML -f b -t RAxML_result."+alignfileprefix+"ml -z RAxML_bootstrap."+alignfileprefix+"boot -m "+model+" -s "+alignfile+" -n "+alignfileprefix)
	else:
		os.system("RAxML -f a -x "+str(random.randrange(1,99999))+" -# "+str(bootstrap)+" -m "+model+" -s "+alignfile+" -n "+alignfileprefix)


###############################################################
# Recursive algorithm to reconstruct CDSs from EMBL CDS lines #
###############################################################

def getCDSseq(embldata, CDSstring, direction, CDScount, feature, CDSposns, CorP):
	
	if '(' in CDSstring:
		if CDSstring.split('(')[0]=='complement':
			if direction=='f':
				direction='r'
			else:
				direction='f'
		parts='('.join(CDSstring.split('(')[1:])
		parts=')'.join(parts.split(')')[:-1])
		getCDSseq(embldata, parts, direction, CDScount, feature, CDSposns, CorP)
		
	else:
		joinlist=CDSstring.split(',')
		for join in joinlist:
			a=int(join.split('..')[0])
			b=int(join.split('..')[1])
			if direction=='f':
				start=a
				end=b+1
				y=1
			else:
				start=b
				end=a-1
				y=-1
			for x in range(start,end,y):
				
				embldata[donelen+x]=CorP
				if words[1]=='CDS':
					CDSposns[CDScount].append((x, direction))


	return embldata, CDSposns






############################################################################
# Get amino acid sequence from concatnated CDS sequences and snp locations #
############################################################################

def getAAsequence(CDS, geneticcode, snplocations):
	
	AAseq=""
	nextsnplocation=0
	
	for x in range(0,len(CDS),3):
		if x+3>len(CDS):
			break
		codon= CDS[x-1]+CDS[x]+CDS[x+1]
		
		if '-' in codon:
			AAseq=AAseq+'-'
		else:
			AAseq=AAseq+geneticcode[codon]
				
			
	return AAseq





#####################
# SNPanalysis class #
#####################

class SNPanalysis:
	def __init__(self, fastq='', name='', mapped='n', CDSseq=''):
		self.fastq=fastq
		self.name=name
		self.runname=''
		self.mapped=mapped
		self.CDSseq=CDSseq
		self.dNdSstats={}
		self.SNPtypes={}
		self.AAtypes={}
		self.goodlen=0
		self.CDSseq=''
		self.snpsummary={}
		self.nummapped=0
		self.percentmapped=0.0
		self.intragenic=0
		self.intergenic=0
		self.SNPs=0
		self.sequence=""
		self.AAsequence=""

def replace_all(text, outlist, incharacter):
	for f in outlist:
		text=text.replace(f, incharacter)
	
	return text


########
# Main #
########


if __name__ == "__main__":
        argv=sys.argv[1:]
        ref, inputdirs, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, quality, mindepth, pairedend, chisquared, fastboot, rtype, multiplexnamed=getOptions(argv)		

	geneticcode={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

	snps={}
	
	refbases={}
	bases={}#0=A 1=C 2=G 3=T
	nostates=0
	converter={}
	convertback={}
	snpbases=[]
	
	pools=[]
	count=0
	poolsort=[]
	runpileups='n'
	
	print '\nChecking input files...'
	sys.stdout.flush()
	
	for pool in inputdirs:
		
		if not os.path.isfile(pool):
			print "File "+pool+" not found! Exiting...\n"
			continue
		
		if pool.split('.')[-1]=='gz':
			pool='.'.join(pool.split('.')[:-1])
		if pool[-1]=='/':
			pool=pool[:-1]
			
			
		pool='.'.join(pool.split('/')[-1].split('.')[:-1])
		
		
		if multiplexnamed=='y':
			if rtype=='solexa' and pairedend=='y':
				if pool.split('_')[-2]=='1':
					pool='_'.join(pool.split('_')[:2])+'_'+pool.split('_')[-1]
				else:
					continue
		
		else:
			if pairedend=='y':
				if pool[-2:]=='_1':
					pool=pool[:-2]
				else:
					continue		

		name=pool
		pool=pool+'_ssaha'

		
		print pool+'...',
		sys.stdout.flush()
		if not os.path.isdir(pool):
		 	print "Directory "+pool+" not found!"
		 	continue
	
		userinput='x'

		if not os.path.isfile(pool+'/all.snp') and not os.path.isfile(pool+'/all.snp.gz'):
			print "\nError: No Ssaha snp file ("+pool+"/all.snp or "+pool+"/all.snp.gz) found!"
			continue

		elif not os.path.isfile(pool+'/all.pileup') and not os.path.isfile(pool+'/all.pileup.gz'):
			print "\nError: No Ssaha snp file ("+pool+"/all.pileup or "+pool+"/all.pileup.gz) found!"
			continue
	
				
		pools.append(SNPanalysis())
		pools[count].runname=pool
		if pools[count].runname=='':
			print "Error: Illegal run name!"
			sys.exit()
		pools[count].name=name
		poolsort.append(pool)
		print 'ok'
		sys.stdout.flush()
		count=count+1
		
		
	if len(pools)==0:
		print "\nError: No valid input folders!"
		sys.exit()
	
#	if snptype=='ssaha':
#		mindepth=float(mindepth)/5
	
	

	
	#Loading reference sequence
	if not os.path.isfile(ref):
		print "\nError: Cannot open reference fasta file "+ref
		sys.exit()
	reflines=open(ref, "rU").readlines()
	refseq=''
	for line in reflines:
		if len(line.strip())>0 and line.strip()[0]!='>':
			refseq=refseq+line.strip().upper()
	
	refseq=replace_all(refseq,['N','W','R','S','B','Y','D','K','N','H','M','V'], '-')
	
	ref=ref.split('/')[-1].split('.')[0]
	reflen=len(refseq)
	reflennogaps=len(refseq.replace('-',''))
	
	# Reading embl file
	
	if embl!='':
		print "\nReading EMBL file...",
		sys.stdout.flush()
		CDSposns={}
		CDSnames={}
		embldata=[]
		donelen=0
		CDScount=0
		pseudoposns={}
		pseudonames={}
		pseudocount=0
		for i in range(reflennogaps+1):
			embldata.append('I')
	
		lines=open(embl,'rU').readlines()
		
		for x, line in enumerate(lines):
			words=line.strip().split()
			if len(words)>2 and words[1] in ('CDS', 'rRNA', 'tRNA') and '..' in words[2]:
				#print words
				location=words[2].replace('<','').replace('>','')
				extralocationlines=1
				foundpool='n'
				y=x+1
				locustag=''
				pseudo='n'
				while foundpool=='n' and y<len(lines) and (len(lines[y].split())<=1 or lines[y].split()[1] not in ('CDS', 'rRNA', 'tRNA', 'gene', 'misc_feature', 'sig_peptide', 'repeat_unit', 'repeat_region', 'misc_RNA')):
					if len(lines[y].split())<=1:
						continue
					if lines[y].split()[1][0]!='/' and y==x+extralocationlines:
						location=location+lines[y].split()[1].strip()
						extralocationlines=extralocationlines+1
					if '/colour=11' in lines[y].split()[1]:
						pseudo='y'
					if '/systematic_id' in lines[y].split()[1] or '/locus_tag' in lines[y].split()[1]:
						locustag=lines[y].split()[1].split('"')[1]
						#print locustag
						foundpool='y'
					y=y+1
				
				direction='f'
				feature=words[1]
				if feature=='CDS' and pseudo=='n':
					CDScount=CDScount+1
					CDSposns[CDScount]=[]
					CDSnames[CDScount]=locustag
				elif feature=='CDS':
					pseudocount=pseudocount+1
					pseudoposns[pseudocount]=[]
					pseudonames[pseudocount]=locustag
					embldata, pseudoposns=getCDSseq(embldata, location, direction, pseudocount, feature, pseudoposns, 'P')
				if pseudo=='n'and feature in ['CDS', 'rRNA', 'tRNA']:
					embldata, CDSposns=getCDSseq(embldata, location, direction, CDScount, feature, CDSposns, feature[0])
				
	
		
		snptoCDSname={}
		for i in CDSposns.keys():
			for j in CDSposns[i]:
				snptoCDSname[j[0]]=CDSnames[i]
		
		
#		print pseudoposns.keys()
#		print pseudonames
		
		snptopseudoname={}
		for i in pseudoposns.keys():
			for j in pseudoposns[i]:
				snptopseudoname[j[0]]=pseudonames[i]
			
		print "Done"
		sys.stdout.flush()
	
	
	poolsort.sort()
	
	avcoverage=[0.0]*(reflen+1)
	covcount=[0]*(reflen+1)
	
	#hassnps={}
	allsnps={}
	
	alssnpsummary={}
	
	print "\nIdentifying SNPs and % of reference mapped..."
	sys.stdout.flush()

	for pool in pools:
		print pool.name+'...',
		sys.stdout.flush()
		tmpout=open(pool.runname+'/snp_sites.txt', "w")
		if os.path.isfile(pool.runname+'/all.snp'):
			lines=open(pool.runname+'/all.snp', 'rU').readlines()
		elif os.path.isfile(pool.runname+'/all.snp.gz'):
			lines=gzip.open(pool.runname+'/all.snp.gz', 'r').readlines()
		if graphs=='y':
			poolcovout=gzip.open(pool.runname+'/'+pool.name+'_q'+str(quality)+'_cov.plot.gz', 'w')
			poolcovynout=gzip.open(pool.runname+'/'+pool.name+'_q'+str(quality)+'_covyn.plot.gz', 'w')
		snpcount=0
		tempsnps={}
		tempseq=['-']*reflen
		tempsnpsummary={'A':{'C':0, 'G':0, 'T':0}, 'C':{'A':0, 'G':0, 'T':0}, 'G':{'C':0, 'A':0, 'T':0}, 'T':{'C':0, 'G':0, 'A':0}}
		hetsnpcount=0
		for line in lines:
			words=line.split()
			het_hez=words[0].split('_')[1].replace(':','')
			snplocation=int(words[3])
			refbase=words[5]
			cnsbase=words[6]
			snpscore=int(words[2])
			snpdepth=int(words[4])
			tempseq[snplocation-1]='?'
			#tempseq[snplocation-1]='-'#add this back in to get a really fast list of snps, but no other output
			if refbase not in ["A","C","G","T","a","c","t","g"]:
				tempseq[snplocation-1]='-'
				continue
			
			if het_hez=='hez':
				hetsnpcount=hetsnpcount+1
			
			Big={'A':int(words[7])-int(words[13]), 'C':int(words[8])-int(words[14]),'G':int(words[9])-int(words[15]),'T':int(words[10])-int(words[16])}
			Small={'A':int(words[13]), 'C':int(words[14]), 'G':int(words[15]), 'T':int(words[16])}
			
			if snpscore<quality or snpdepth<mindepth:
				continue
			
			if het_hez=='hez':
				bestscore=0
				bestsnp=''
				for x in Big.keys():
					if Big[x]>=bestscore:
						bestsnp=x
						bestscore=Big[x]
					
				if float(bestscore)/snpdepth<0.75:
					continue
				elif bestsnp!=refbase:
					cnsbase=bestsnp
				elif bestsnp==refbase:
					cnsbase=refbase
					continue

			
			if allsnps.has_key(snplocation):
				allsnps[snplocation][0]=allsnps[snplocation][0]+1
				allsnps[snplocation][2][cnsbase]=allsnps[snplocation][2][cnsbase]+1
			else:
				allsnps[snplocation]=[1,refbase,{'A':0, 'C':0, 'G':0, 'T':0}]
				allsnps[snplocation][2][cnsbase]=1
			tempsnpsummary[refbase][cnsbase]=tempsnpsummary[refbase][cnsbase]+1
			tempseq[snplocation-1]=cnsbase
			print >> tmpout, snplocation
			snpcount=snpcount+1

		
		#snpcount=len(tempsnps.keys())
		
		print '%d SNPs found,' % (snpcount),
		sys.stdout.flush()
		tmpout.close()
		#pool.sequence=tempseq;continue #add this back in to get a really fast list of snps, but no other output
		
		if os.path.isfile(pool.runname+'/all.pileup'):
			lines=open(pool.runname+'/all.pileup', 'rU')
		elif os.path.isfile(pool.runname+'/all.pileup.gz'):
			lines=gzip.open(pool.runname+'/all.pileup.gz', 'r')
		mapped=0.0
		hetcount=0
		#for x in hassnps.keys():
			#line=lines[x-1]
			
		if quality/5>mindepth:
			mindepth=quality/5
		
		for line in lines:
			#print line
			words=line.split()
			snplocation=int(words[2])
			refbase=words[4]
			if refbase not in ["A","C","G","T","a","c","t","g"]:
				tempseq[snplocation-1]=='-'
				continue
			readdepth=int(words[3])
			
			baselocs={"A":5, "C":6, "G":7, "T":8}
			
			#remove from here down to next commented line (which needs to be added back in) to get old version
			locqual=0.0
			if readdepth>=mindepth and tempseq[snplocation-1]!='?':
				#print line
				locqual=((int(words[baselocs[refbase]])-int(words[baselocs[refbase]+6]))*5)+((float(words[baselocs[refbase]+6]))*2.5)
		
			avcoverage[snplocation]=avcoverage[snplocation]+readdepth
	
			if locqual>=quality:
	#		if readdepth>=mindepth and tempseq[snplocation-1]!='?':
				
				if tempseq[snplocation-1]=='-':
					tempseq[snplocation-1]=refbase
				covcount[snplocation]=covcount[snplocation]+1
				mapped=mapped+1
				
				
				basefoundcount=0
				for base in words[5:9]:
					if int(base)!=0:
						basefoundcount=basefoundcount+1
				if basefoundcount>1:
					hetcount=hetcount+1
				if graphs=='y':
					print >> poolcovynout, 1
			else:
				if tempseq[snplocation-1]=='?':
					tempseq[snplocation-1]='-'
				if graphs=='y':
					print >> poolcovynout, 0
			if graphs=='y':
				print >> poolcovout, readdepth
			if tempseq[snplocation-1]=='?':
				print snplocation
			
		if graphs=='y':
			poolcovout.close()
			poolcovynout.close()
		pool.SNPs=snpcount
		pool.sequence=tempseq
		pool.mapped=tempsnps
		pool.snpsummary=tempsnpsummary
		pool.nummapped=mapped
		pool.hetcount=hetcount
		pool.hetsnpcount=hetsnpcount
		pool.percentmapped=((mapped/reflen)*100)
		print '%.2f%% of reference mapped' % (pool.percentmapped)
		sys.stdout.flush()
	
	
	
	#Chi-squared looking for snp clustering. Window size=1000
	if chisquared=='y':	
		print "\nCalculating chi-squared for snps..."
		sys.stdout.flush()
		
		snpsort=allsnps.keys()
		snpsort.sort()
		window=100
		windowcounts=[0]*(reflennogaps)
		X2=[]
		pvalues=[]
		currwindowcount=0
		currposn=0
		X2converter={}
		
		expectedsnps=float(len(snpsort))/reflennogaps
		expectednonsnps=window*(1-expectedsnps)
		expectedsnps=expectedsnps*window
		
		while (expectedsnps<10 or expectednonsnps<10) and window<(reflen/1000):
			print "Expected frequencies too low for chi-squared test using a window size of "+str(window)+" (",expectedsnps,"and", expectednonsnps,"). Increasing window size to "+str(window*10)
			expectedsnps=expectedsnps*10
			expectednonsnps=expectednonsnps*10
			window=window*10
			
			
		print "Calculating moving window snp counts..."
		sys.stdout.flush()
		
		for x in snpsort:
			if refseq[x]!='-':
				y=x-(window/2)
				curposn=y
				while (y+window)>curposn:
					if curposn<0:
						actualposn=reflennogaps+curposn
					elif curposn>(reflennogaps-1):
						actualposn=curposn-reflennogaps
					else:
						actualposn=curposn
						
					windowcounts[actualposn]=windowcounts[actualposn]+1
	
					curposn=curposn+1
		
		
		print "Calculating chi-squared values..."
		sys.stdout.flush()	
				
		for x in windowcounts:
			if not X2converter.has_key(x):
				X2converter[x]=[ (((x-expectedsnps)**2)/expectedsnps) + ((((window-x)-expectednonsnps)**2)/expectednonsnps),0]
				if (window-x)>expectednonsnps:
					X2converter[x][1]=(1-chi2.cdf(X2converter[x][0],1))
				else:
					X2converter[x][1]=1
			
			X2.append(X2converter[x][0])
			pvalues.append(X2converter[x][1])	
		
		
		print "\nPrinting chi-squared files..."
		sys.stdout.flush()
		
		output=open(outfile+'_chi_squared_pvalue.plot','w')
		outputb=open(outfile+'_chi_squared.plot','w')
		outputc=open(outfile+'_moving_window_snps'+str(window)+'.plot','w')
		
		for x, pvalue in enumerate(pvalues):
			print >> output, pvalue
			print >> outputb, X2[x]
			print >> outputc, windowcounts[x]
			
		output.close()
		outputb.close()
		outputc.close()
		
		print "Done."

	#DNDS stuff
	comp={'A':'T','T':'A','G':'C','C':'G', '-':'-'}
	
	if embl!='' and ref!='':
		
		print "\nLooking for SNPs in pseudogene stop codons..."
		sys.stdout.flush()
		for i in pseudoposns.keys():
			dirn=pseudoposns[i][0][1]
			startposn=pseudoposns[i][0][0]
			x=startposn
			refpseudoseq=''
			while embldata[x] in ['P', 'I']:
				if dirn=='f':
					refpseudoseq=refpseudoseq+refseq[x-1].upper()
					x=x+1
				else:
					refpseudoseq=refpseudoseq+comp[refseq[x-1].upper()]
					x=x-1
					
			for pool in pools:
				if pool.name==ref:
					continue
				
				x=startposn
				querypseudoseq=''
				while embldata[x] in ['P', 'I']:
				
					if pool.sequence[x-1]=='-':
						poolbase=refseq[x-1]
					else:
						poolbase=pool.sequence[x-1]
				
					if dirn=='f':
						querypseudoseq=querypseudoseq+poolbase.upper()
						x=x+1
					else:
						querypseudoseq=querypseudoseq+comp[poolbase.upper()]
						x=x-1
				
				if refpseudoseq!=querypseudoseq:
					for x in range(0, len(refpseudoseq), 3):
						if len(refpseudoseq)<(x+3):
							break
						refcodon=refpseudoseq[x]+refpseudoseq[x+1]+refpseudoseq[x+2]
						querycodon=querypseudoseq[x]+querypseudoseq[x+1]+querypseudoseq[x+2]
						if refcodon!=querycodon:
							if geneticcode[refcodon]=='*':
								print 'SNP in pseudogene stop codon in', pool.name, pseudonames[i], refcodon, '->', querycodon, '('+geneticcode[querycodon]+')'
								sys.stdout.flush()
						

	
	
	if embl!='' and ref!='':
		print "\nCalculating dN/dS values...\nConcatenating reference CDSs...",
		sys.stdout.flush()
		
		dnbydsstats={}
		
		dndsout=open(outfile+'_dnds.out','w')
		dnout=open(outfile+'_dn.out','w')
		dsout=open(outfile+'_ds.out','w')
		zout=open(outfile+'_dnds_z.out','w')
		
		
		
		print >> dndsout, 'Total',
		print >> dnout, 'Total',
		print >> dsout, 'Total',
		print >> zout, 'Total',
			
		refCDSseq=''
		CDSbasenumbers=[]
		for key in CDSposns.keys():
			
			if float(len(CDSposns[key]))/3!=len(CDSposns[key])/3:
				print "\nError! CDS "+CDSnames[key]+" has a length that is not divisible by 3. Is it a pseudogene? Skipping..."
				continue
		
			for x in CDSposns[key]:
				CDSbasenumbers.append(x[0]-1)
				if x[1]=='f':
					refCDSseq=refCDSseq+refseq[x[0]-1].upper()
				else:
					refCDSseq=refCDSseq+comp[refseq[x[0]-1].upper()]
		
		
		print 'Done'
		sys.stdout.flush()
		tmpseqs={}
		snptypes={}
		
		
		
		for poolnum, pool in enumerate(pools):
			if pool.name==ref:
				#if align=='y':
					#pool.AAsequence=getAAsequence(refCDSseq, geneticcode, allsnps.keys())
					#refAAseqlen=len(pool.AAsequence)
				continue
			
			print pool.name+'...',
			sys.stdout.flush()
			
			tmpCDSseq=''
			
			for key in CDSposns.keys():
				if float(len(CDSposns[key]))/3!=len(CDSposns[key])/3:
					continue
				for x in CDSposns[key]:
					if x[1]=='f':
						tmpCDSseq=tmpCDSseq+pool.sequence[x[0]-1].upper()
					else:
						tmpCDSseq=tmpCDSseq+comp[pool.sequence[x[0]-1].upper()]
			
			#if align=='y':
			#	pool.AAsequence=getAAsequence(tmpCDSseq, geneticcode, allsnps.keys())
			
			pool.dNdSstats, pool.SNPtypes, pool.AAtypes=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers)
			
			
			
			
#Add this back in to get pairwise dNdS...very slow
			#dndstempoutput=open("pairwise_dNdS.txt","w")	
			
#			for poolb in pools[poolnum:]:
#				if pool.name==poolb.name:
#				#if align=='y':
#					#pool.AAsequence=getAAsequence(refCDSseq, geneticcode, allsnps.keys())
#					#refAAseqlen=len(pool.AAsequence)
#					continue
#				
#				sys.stdout.flush()
#				
#				tmpCDSseqb=''
#				
#				for keyb in CDSposns.keys():
#					if float(len(CDSposns[keyb]))/3!=len(CDSposns[keyb])/3:
#						continue
#					for y in CDSposns[keyb]:
#						if y[1]=='f':
#							tmpCDSseqb=tmpCDSseqb+poolb.sequence[y[0]-1].upper()
#						else:
#							tmpCDSseqb=tmpCDSseqb+comp[pool.sequence[y[0]-1].upper()]
#				
#				
#				dNdSstats, SNPtypes, AAtypes=dnbyds(tmpCDSseq, tmpCDSseqb, CDSbasenumbers)
#				if dNdSstats['dS']!=0:
#					print >> dndstempoutput, pool.name, poolb.name, "%.2f" % (dNdSstats['dN']/dNdSstats['dS'])
#				else:
#					print >> dndstempoutput, pool.name, poolb.name, "-"





			
			#N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd
			if pool.dNdSstats['dS']!=0:
				print "%.2f" % (pool.dNdSstats['dN']/pool.dNdSstats['dS'])
				print >> dndsout, '\t'+str(pool.dNdSstats['dN']/pool.dNdSstats['dS']),
				print >> dnout, '\t'+str(pool.dNdSstats['dN']),
				print >> dsout, '\t'+str(pool.dNdSstats['dS']),
				print >> zout, '\t'+str(pool.dNdSstats['z']),
				sys.stdout.flush()
			else:
				print '-'
			
		#dndstempoutput.close()	
			
		#Testing doing each gene seperately

		print "\nCalculating CDS dN/dS values...",
		
		
		

		for key in CDSposns.keys():
		
			refCDSseq=''
			
			if float(len(CDSposns[key]))/3!=len(CDSposns[key])/3:
				continue
			
			for x in CDSposns[key]:
				if x[1]=='f':
					refCDSseq=refCDSseq+refseq[x[0]-1].upper()
				else:
					refCDSseq=refCDSseq+comp[refseq[x[0]-1].upper()]
			
			print >> dndsout, '\n'+CDSnames[key],
			print >> dnout, '\n'+CDSnames[key],
			print >> dsout, '\n'+CDSnames[key],
			print >> zout, '\n'+CDSnames[key],
			#print CDSnames[key]
			#sys.stdout.flush()
			
			for pool in pools:
				if pool.name==ref:
					continue
					
				tmpCDSseq=''
				
				for y in CDSposns[key]:
					if y[1]=='f' and pool.sequence[y[0]-1]!='-':
						tmpCDSseq=tmpCDSseq+pool.sequence[y[0]-1].upper()
					elif pool.sequence[y[0]-1]!='-':
						tmpCDSseq=tmpCDSseq+comp[pool.sequence[y[0]-1].upper()]
					else:
						tmpCDSseq=tmpCDSseq+'-'
				#print refCDSseq+'\n'+tmpCDSseq+'\n'
				dnbydsstats, gene_snptypes, gene_AAs=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers)
								
				if dnbydsstats['dS']!=0:
					print >> dndsout, '\t'+str(dnbydsstats['dN']/dnbydsstats['dS']),
					print >> dnout, '\t'+str(dnbydsstats['dN']),
					print >> dsout, '\t'+str(dnbydsstats['dS']),
					print >> zout, '\t'+str(dnbydsstats['z']),
				else:
					print >> dndsout, '\t-',
					print >> dnout, '\t'+str(dnbydsstats['dN']),
					print >> dsout, '\t'+str(dnbydsstats['dS']),
					print >> zout, '\t'+str(dnbydsstats['z']),
				dndsout.flush()
				
		
		
		dndsout.close()
		dnout.close()
		dsout.close()
		zout.close()
			
			
		print "Done"
		sys.stdout.flush()
	
	#print summary file and alignment file
	
	print "\nWriting output file(s)...",
	sys.stdout.flush()
	
	output=open(outfile+'.out','w')
	
	#if len(allsnps.keys())>1:
	#	print >> output, 'Ref_sequence\tPosition_in_ref',
	#else:
	print >> output, 'Position_in_'+ref,
	if embl!='':
	        print >> output, '\tCDS/rRNA/tRNA/Intergenic\tCDS_name\tSynonymous/Non-synonymous',
	print >> output, '\tRef_base\tSNP_base\tTotal',
	
	for pool in pools:
		print >> output, '\t'+pool.name,
	print >> output, '\tPattern\n',
	
	
	if tabfile=='y':
	        tabout=open(outfile+'_snps.tab','w')
	if graphs=='y':
		snpgraph=gzip.open(outfile+'_snps.plot.gz','w')
		covgraph=gzip.open(outfile+'_cov.plot.gz','w')
		covcountgraph=gzip.open(outfile+'_covcount.plot.gz','w')
	
	#minquality=10
	#minqualityperdepth=4
	
	snpsequence={}
	straintabfiles={}
	for pool in pools:
		snpsequence[ref]=''
		snpsequence[pool.name]=''
		if tabfile=='y':
			straintabfiles[pool.name]=open(pool.runname+'/snps.tab','w')
	
	if ref=='':
		ref=pools[0].name
	
	if tabfile=='y':
		print >> tabout, 'ID   SNP'
	
	poscount=1
	
	
	allaln=open(outfile+'_whole.aln','w')
	for pool in pools:
		print >> allaln, '>'+pool.name
		print >> allaln, ''.join(pool.sequence)
	
	allaln.close()
	
	snpsort=allsnps.keys()
	snpsort.sort()
				
	for i in snpsort:
		j=i-1
		numbases=0
		foundbases={}	
		for pool in pools:
			base=pool.sequence[j].upper()
			if base not in foundbases.keys() and base!='-':
				foundbases[base]=1
				numbases=numbases+1
			elif base!='-':
				foundbases[base]=foundbases[base]+1
		if numbases>0:
			snpbases.append(foundbases)
			
				
	
	for x, i in enumerate(snpsort):
		j=i-1
	#for x, j in enumerate(snplocations):
	
		if graphs=='y':
			while j>poscount:
				print >> snpgraph, '0'
				print >> covgraph, str(float(avcoverage[poscount])/len(pools))
				print >> covcountgraph, str(covcount[poscount])
				poscount=poscount+1
		
		
		outstring=''
		tabstring=''
		
		#if len(allsnps.keys())>1:
		#	outstring=outstring+key+'\t'+str(j)
		#else:
		outstring=outstring+str(j+1)
		snpcolour='1'
		if embl!='':
			
			if embldata[j]!='C':
				outstring=outstring+'\t'+embldata[j]+'\t-'
			else:
				outstring=outstring+'\t'+embldata[j]+'\t'+snptoCDSname[j]
			snptype='-'
			for pool in pools:
				if pool.name==ref:
					continue
				if pool.SNPtypes.has_key(j):
					if pool.SNPtypes[j]=='S':
						snpcolour='3'
					elif pool.SNPtypes[j]=='N':
						snpcolour='2'
					elif pool.SNPtypes[j]=='2':
						snpcolour='4'
					elif pool.SNPtypes[j]=='3':
						snpcolour='4'
					if snptype=='-':
						snptype=pool.SNPtypes[j]
					elif pool.SNPtypes[j] not in snptype:
						snptype=snptype+'/'+pool.SNPtypes[j]
				elif snptopseudoname.has_key(j):
					snpcolour='11'
					pool.SNPtypes[j]="P"
			outstring=outstring+'\t'+snptype
		
		
		outstring=outstring+'\t'+refseq[j]
		refbase=refseq[j]
		snpsequence[ref]=snpsequence[ref]+refseq[j]
			
		snpbase=''
		snpbasecount=0
		for base in snpbases[x].keys():
			if base!=refbase and snpbases[x][base]>snpbasecount:
				snpbasecount=snpbases[x][base]
				snpbase=base
			elif base!=refbase and  snpbases[x][base]==snpbasecount:
				snpbase=snpbase+','+base
	
		outstring=outstring+'\t'+snpbase+'\t'+str(snpbasecount)# need to think how to do this
		
		if tabfile=='y':
			tabstring=tabstring+'FT   SNP             '+str(j+1)+'\n'
			tabstring=tabstring+'FT                   /note="refAllele: '+refseq[j]
			#tabstring=tabstring+'FT                   /SNPAllele="'+snpbase+'"\n'
			tabstring=tabstring+' SNPstrains: '
		
		
	
		pattern=''
		numsnpbases=1
		sitesnpbases=[]
		for pool in pools:
			if pool.sequence[j]!='-':# and pool.SNPs[key][j][1]>8:# and pool.SNPs[key][j][2]>(pool.SNPs[key][j][3]*4):
				if pool.sequence[j]!=refseq[j]:
					
					if tabfile=='y':
					
						print >> straintabfiles[pool.name], 'FT   SNP             '+str(j+1)+'\nFT                   /colour='+snpcolour
					
						tabstring=tabstring+pool.name+'='+pool.sequence[j]+' '
						if pool.SNPtypes.has_key(j):
							if pool.SNPtypes[j]=='S':
								tabstring=tabstring+'(synonymous) '
							elif pool.SNPtypes[j]=='N':
								if pool.AAtypes.has_key(j):
									tabstring=tabstring+'(non-synonymous, AA='+pool.AAtypes[j]+') '
								else:
									tabstring=tabstring+'(non-synonymous) '
							elif pool.SNPtypes[j]=='P':
								tabstring=tabstring+'(SNP in pseudogene) '
							elif pool.SNPtypes[j]=='1':
								tabstring=tabstring+'(gap in SNP codon) '
							elif pool.SNPtypes[j]=='2':
								tabstring=tabstring+'(SNP codon is STOP) '
							elif pool.SNPtypes[j]=='3':
								tabstring=tabstring+'(ref codon is STOP) '
					pattern=pattern+'1'
					outstring=outstring+'\t'+pool.sequence[j]#+str(pool.mapped[key][j][1])
					
				else:
					pattern=pattern+'0'
					outstring=outstring+'\t'+pool.sequence[j]
				
			else:
				pattern=pattern+'-'
				#snpsequence[pool]=snpsequence[pool]+'-'
				outstring=outstring+'\t-'
			
			snpsequence[pool.name]=snpsequence[pool.name]+pool.sequence[j]
		
		
		if graphs=='y':
			print >> snpgraph, '1'
			snpgraph.flush()
		print >> output, outstring+'\t'+pattern
		output.flush()
		if tabfile=='y':
			print >> tabout, tabstring+'"\nFT                   /colour='+snpcolour
			tabout.flush()
	
		if graphs=='y':
			print >> covgraph, str(float(avcoverage[poscount])/len(pools))
			covgraph.flush()
			print >> covcountgraph, covcount[poscount]
			covcountgraph.flush()
			poscount=poscount+1
	
	
	if graphs=='y':
		while reflen>=poscount:
			print >> snpgraph, '0'
			print >> covgraph, str(float(avcoverage[poscount])/len(pools))
			print >> covcountgraph, str(covcount[poscount])
			poscount=poscount+1
	
	
	if align=='y':
		aln={}
		alnfile=open(outfile+'_NUC.aln', 'w')
		count=0
		inaln=[]
		for pool in snpsequence.keys():
			if len(snpsequence[pool].replace('-',''))>(len(snpsequence[ref])/2):
				inaln.append(pool)
				count=count+1
		
		
		
		print >> alnfile, count, len(snpsequence[ref])
		#print >> alnfile, ref+' '*(10-len(ref))+refseq
		for pool in inaln:
			if len(pool.replace('pool','').split('.')[0])<9:
				print >> alnfile, pool.replace('pool','').split('.')[0]+' '*(10-len(pool.replace('pool','').split('.')[0]))+snpsequence[pool]
			else:
				print >> alnfile, pool.replace('pool','').split('.')[0]+' '+snpsequence[pool]
	
					
		alnfile.close()
		
#		alnfileb=open(outfile+'_AA.aln', 'w')
#		count=0
#		inaln=[]
#		for pool in pools:
#			if len(pool.AAsequence.replace('-',''))>(refAAseqlen/2):
#				inaln.append(pool)
#				count=count+1
#		
#		print >> alnfileb, count, len(snpsequence[ref])
#		for pool in inaln:
#			print >> alnfile, pool.replace('pool','').split('.')[0]+' '*(10-len(pool.replace('pool','').split('.')[0]))+pool.AAsequence
#		
#		
	
	
	if tabfile=='y':
		tabout.close()
		for pool in pools:
			straintabfiles[pool.name].close()
	
	if graphs=='y':
		snpgraph.close()
		covgraph.close()
		covcountgraph.close()
	output.close()
	
	
	
	summaryfile=open(outfile+'_summary.out','w')
	
	print >> summaryfile, 'Strain\tPercent of Reference Mapped\tN\tS\tpN\tpS\tdN/dS',
	
	for x in ['A', 'C', 'G', 'T']:
		for y in ['A', 'C', 'G', 'T']:
			if x!=y:
				print >> summaryfile, '\t'+x+'->'+y,
	
	print >> summaryfile, '\tTotal SNPs',
	
	if embl!='':
		print >> summaryfile, '\t% Mapped Intragenic sites with SNP\t% Mapped Intergenic sites with SNP'
	else:
		print >> summaryfile, '\n',
	
	for pool in pools:
		intragenic=0.0
		intragenic_mapped=0.0
		if embl!='':
			if pool.dNdSstats['dS']!=0:
				print >> summaryfile, pool.name+'\t'+str(pool.percentmapped)+'\t'+str(pool.dNdSstats['N'])+'\t'+str(pool.dNdSstats['S'])+'\t'+str(pool.dNdSstats['pN'])+'\t'+str(pool.dNdSstats['pS'])+'\t'+str(pool.dNdSstats['dN']/pool.dNdSstats['dS']),
			else:
				print >> summaryfile, pool.name+'\t'+str(pool.percentmapped)+'\t'+str(pool.dNdSstats['N'])+'\t'+str(pool.dNdSstats['S'])+'\t'+str(pool.dNdSstats['pN'])+'\t'+str(pool.dNdSstats['pS'])+'\t-',
		else:
			print >> summaryfile, pool.name+'\t'+str(pool.percentmapped)+'\t-\t-\t-\t-\t-',
		for x in ['A', 'C', 'G', 'T']:
			for y in ['A', 'C', 'G', 'T']:
				if x!=y:
					print >> summaryfile, '\t'+str(pool.snpsummary[x][y]),
					count=count+pool.snpsummary[x][y]
		
		if embl!='':
			for CDS in CDSposns.keys():
				for y in CDSposns[CDS]:
					y=y
					if pool.sequence[y[0]-1]!=refseq[y[0]-1] and pool.sequence[y[0]-1]!='-':
						intragenic=intragenic+1
					if pool.sequence[y[0]-1]!='-':
						intragenic_mapped=intragenic_mapped+1
				#print pool.mapped[y[0]]
			#print pool.nummapped, intragenic_mapped, len(pool.mapped.keys())
			intergenic_mapped=pool.nummapped-intragenic_mapped
		
		#print intragenic_mapped, intergenic_mapped, pool.nummapped
		
		
			if intergenic_mapped>0:
				intergenic=(float(pool.SNPs-intragenic)/intergenic_mapped)*100
			else:
				intergenic=0
			if intragenic_mapped>0:
				intragenic=(float(intragenic)/intragenic_mapped)*100
			else:
				intragenic=0
		
		print >> summaryfile, '\t'+str(pool.SNPs),
		
		if embl!='':
			print >> summaryfile, '\t'+str(intragenic)+'\t'+str(intergenic)
		else:
			print >> summaryfile
	
	
	if embl!='':
		print >> summaryfile, '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(len(snpsort))
	else:
		print >> summaryfile, '\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(len(snpsort))
	
	
	print 'Done'
	
	summaryfile.close()	
	
	if raxml=='y':
		RAxMLphylogeny(outfile+'_NUC.aln', model, bootstrap, fastboot)
	
		
	
	
	
	
	print "All finished.\n"
	
	
