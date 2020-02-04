#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re, gzip
import os, sys, getopt, random, math


#########
# Usage #
#########

def Usage():
	os.system('clear')
	print '\nsummarise_snps_from_maq3.py Usage:\n'
	print 'summarise_snps_from_maq3.py [OPTIONS] <list of fastq files>'
	print "\nINPUT OPTIONS:"
	print "-r Reference dna sequence:\t<fasta file> (required)"
	print "-e Reference embl file:\t\t<embl or tab file>"
	print "\nMAPPING OPTIONS:"
	print "-R Run Maq:\t\t\t<y/n>"
	print "-n max no. of mismatches for mapping:\t<integer between 1 and 3>"
	print "-m max no. of snps per read for consensus:\t<integer between 1 and 100>"
	print "-p use paired-end reads:\t<y/n>"
	print "-i max insert size + 2*read length:\t<integer between 100 and 3000> [default=300]"
	print "-q Minimum mapping quality:\t<integer between 1 and 60> [default=30]"	
	print "-d Minimum mapping depth:\t<integer between 1 and 100,000> [default=5]"
	print "-v Assemble non-mapping reads:\t<y/n>"
	print "\nOUTPUT OPTIONS:"
	print "-o Prefix for output files:\t<file name prefix>"
	print "-a Create SNP alignment file:\t<y/n>"
	print "-g Create coverage plots:\t<y/n>"
	print "-t Create SNP tab file:\t\t<y/n>"
	print "-P Run phylogeny with RAxML:\t<y/n>"
	print "-M Model of evolution\t\t<GTR, GTRGAMMA, etc.> [default=GTRGAMMA]"
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
		opts, args = getopt.getopt(argv, "hd:e:r:o:tai:Rcpgb:M:q:s:vIm:Pn:", ["depth=", "model=", "bootstrap=", "help", "phylogeny" "embl=" "ref=", "out=", "align", "maq", "tabfile", "graphs", "indels", "circular", "quality=", "snptype=", "velvet=", "interactive", "maxsnps=", "pairedend", "insert=", "mismatches="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	ref=''
	outfile=''
	inputdirs=[]
	i=0
	runmaq='n'
	tabfile='n'
	indels='n'
	align='n'
	circular='n'
	embl=''
	dirty='n'
	graphs='n'
	raxml='n'
	bootstrap=100
	model='GTRGAMMA'
	quality=30
	depth=5
	snptype='maq'
	velvet='n'
	interactive='n'
	maxsnps=2
	pairedend='n'
	maxinsertsize=300
	mismatches=2

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-R", "--maq"):
			runmaq='n'
		elif opt in ("-t", "--tabfile"):
			tabfile='y'
		elif opt in ("-i", "--indels"):
			indels='y'
		elif opt in ("-c", "--circular"):
			circular='y'
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
		elif opt in ("-q", "--quality"):
			quality=int(arg)
		elif opt in ("-d", "--depth"):
			depth=int(arg)
		elif opt in ("-m", "--maxsnps"):
			maxsnps=int(arg)
		elif opt in ("-s", "--snptype"):
			snptype=arg.lower()
		elif opt in ("-v", "--velvet"):
			velvet='y'
		elif opt in ("-I", "--interactive"):
			interactive='y'
		elif opt in ("-p", "--pairedend"):
			pairedend='y'
		elif opt in ("-i", "--insert"):
			maxinsertsize=int(arg)
		elif opt in ("-n", "--mismatches"):
			mismatches=int(arg)

	inputdirs=args
	
	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=ref.split('.')[0]+"_q"+str(quality)+"_d"+str(depth)
	if align=='n' and raxml=='y':
		align='y'
	
	
	if interactive=='y':
		ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches=menusystem(ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches)

	if ref=='':
		print 'Error: No reference dna file selected!'
		Usage()
		sys.exit()
	
	
	return ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches



####################
# Interactive Menu #
####################

def menusystem(ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches):
	
	os.system('clear')
	
	if outfile==ref.split('.')[0]+"_q"+str(quality)+"_d"+str(depth):
		outorig='y'
	else:
		outorig='n'
	
	print "\nsummarise_snps_from_maq3.py: Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009"
	
	print "\nINPUT OPTIONS:"
	
	if ref=='':
		print "r: Reference dna sequence:\t\tNone selected (required)"
	else:
		print "r: Reference dna sequence:\t\t"+ref
		if embl=='':
			print "e: Reference embl file:\t\t\tNone selected"
		else:
			print "e: Reference embl file:\t\t\t"+embl
		print "\nMAPPING OPTIONS:"
		if runmaq=='n':
			print "R: Run Maq:\t\t\t\tno, use output files from previous run"
		else:
			print "R: Run Maq:\t\t\t\tyes"
		if pairedend=='n':
			print "p: Use paired-end reads:\t\tno"
		else:
			print "p: Use paired-end reads:\t\tyes"
			print "i: Maximum insert size (inc reads):\t"+str(maxinsertsize)
		print "n: Maximum mismatches:\t\t\t"+str(mismatches)
		print "m: Maximum SNPs per read:\t\t"+str(maxsnps)
		print "q: Minimum mapping quality:\t\t"+str(quality)
		print "d: Minimum mapping depth:\t\t"+str(depth)
		if velvet=='n':
			print "v: Assemble non-mapping reads:\t\tno"
		else:
			print "v: Assemble non-mapping reads:\t\tyes"
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
		if raxml=='n':
			print "P: Run phylogeny with RAxML:\t\tno"
		else:
			print "P: Run phylogeny with RAxML:\t\tyes"
		if raxml=='y':
			print "M:  Model of evolution:\t\t\t"+model
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
		inputlist=['r', 'e', 'R', 'q','d','m','v','o','a','g','t','P', 'y', 'Q', 'D', 'p']
		if pairedend=='y':
			inputlist=inputlist+['i']
		if raxml=='y':
			inputlist=inputlist+['M','b']
			
	ui=''
	while ui not in inputlist:
		ui=raw_input(message+' ')

	if ui=='y':
		os.system('clear')
		return ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches
		
	elif ui=='r':
		oldref=ref
		ref=''
		while not os.path.isfile(ref):
			ref=raw_input('Enter reference file name or Q to go back to the menu: ')
			if ref=='Q':
				ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches=menusystem(oldref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches)
			elif not os.path.isfile(ref):
				print "File not found"
				
	elif ui=='e':
		oldembl=embl
		embl=''
		while not os.path.isfile(embl):
			embl=raw_input('Enter embl file name or Q to go back to the menu: ')
			if embl=='Q':
				ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches=menusystem(ref, inputdirs, outfile, tabfile, align, oldembl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches)
			elif not os.path.isfile(embl):
				print "File not found"
	
	
	elif ui=='R':
		if runmaq=='n':
			runmaq='y'
		else:
			runmaq='n'
	
	elif ui=='n':
		mismatches=0
		while mismatches > 3 or mismatches < 1:
			mismatches=int(raw_input('Enter maximum mismatches (1-3): '))
		
	elif ui=='m':
		maxsnps=0
		while maxsnps > 100 or maxsnps < 1:
			maxsnps=int(raw_input('Enter maximum number of SNPs per read (1-100): '))
	
	elif ui=='q':
		quality=0
		while quality > 60 or quality < 1:
			quality=int(raw_input('Enter minimum mapping quality (1-60): '))
	
	elif ui=='d':
		depth=0
		while depth > 100000 or depth < 1:
			depth=int(raw_input('Enter minimum read depth for mapping and SNP calling (1-100,000): '))
	
	elif ui=='v':
		if velvet=='n':
			velvet='y'
		else:
			velvet='n'
	
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
		maxsnps=2
	
	elif ui=='p':
		if pairedend=='n':
			pairedend='y'
		else:
			pairedend='n'
	
	elif ui=='i':
		maxinsertsize=0
		while maxinsertsize > 10000 or maxinsertsize < 100:
			maxinsertsize=int(raw_input('Enter maximum insert size (including reads) (100-10000): '))
			
	elif ui=='Q':
		sys.exit()
	
	if outorig=='y':
		outfile=ref.split('.')[0]+"_q"+str(quality)+"_d"+str(depth)
		
	ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches=menusystem(ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches)

	return ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, depth, maxsnps, pairedend, maxinsertsize, mismatches


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
		
		if codon!=SNPcodon:
			
			for y,z in enumerate(codon):
			
				if SNPcodon[y]!=z and SNPcodon[y]!='-':
					newSNPcodon=codon[:y]+SNPcodon[y]+codon[y+1:]
			
					if geneticcode[newSNPcodon]=='*':
						SNPtype[codonposn[y]]='2'
					elif geneticcode[newSNPcodon]=='*':
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
					
					
		if '-' in codon or '-' in SNPcodon:
			gapcount=gapcount+3
			continue
		
		if geneticcode[codon]=='*' or geneticcode[SNPcodon]=='*':
			continue
			print codon, x
		
		
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
	N=len(CDS)-gapcount-S
	
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
		dS=(-3*(math.log(1-((pS*4)/3))))/4
		dN=(-3*(math.log(1-((pN*4)/3))))/4
		varianceS=(9 * pS * (1 -pS))/(((3 - 4 *pS) **2) * (len(CDS)-gapcount));
		varianceN=(9 * pN * (1 -pN))/(((3 - 4 *pN) **2) * (len(CDS)-gapcount));
		z=(dN - dS) / math.sqrt(varianceS + varianceN)
		
	else:
		#print "Too divergent for JC! Using pN/pS instead."
		dS=pS
		dN=pN

	
	return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}, SNPtype, AAtype


###############################################################
# Run a bootstrapped phyogeny of an alignment file with RAxML #
###############################################################

def RAxMLphylogeny(alignfile, model, bootstrap, force='y'):
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
	os.system("RAxML -f i -b "+str(random.randrange(1,99999))+" -# "+str(bootstrap)+" -m "+model+" -s "+alignfile+" -n "+alignfileprefix+"boot")
	os.system("RAxML -f d -m "+model+" -s "+alignfile+" -n "+alignfileprefix+"ml")
	os.system("RAxML -f b -t RAxML_result."+alignfileprefix+"ml -z RAxML_bootstrap."+alignfileprefix+"boot -m "+model+" -s "+alignfile+" -n "+alignfileprefix)


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







#####################
# SNPanalysis class #
#####################

class SNPanalysis:
	def __init__(self, fastq='', name='', mapped={}, runmaq='n', CDSseq=''):
		self.fastq=fastq
		self.name=name
		self.runname=''
		self.mapped=mapped
		self.runmaq=runmaq
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
	def runMaq(self, ref, quality, mindepth, maxsnps, pairedend, maxinsertsize, mismatches):
		print "\nRunning Maq on "+self.name+'...'
		
		#os.system("maq fasta2bfa "+ref+" "+self.runname+"/ref.bfa")
		
		
		if pairedend=='n':
			#print "maq fastq2bfq fastqs/"+self.name+".fastq "+self.runname+"/reads_1.bfq"
			#os.system("maq fastq2bfq fastqs/"+self.name+".fastq "+self.runname+"/reads_1.bfq")
			#os.system("maq map -n "+mismatches+" -u "+self.runname+"/unmap_1.txt "+self.runname+"/reads_1.map "+self.runname+"/ref.bfa  "+self.runname+"/reads_1.bfq")
			os.system("maq mapcheck -s -m "+maxsnps+" -q "+quality+" "+self.runname+"/ref.bfa "+pool.runname+"/reads_1.map >"+self.runname+"/mapcheck.txt")
			os.system("maq assemble -s -m "+maxsnps+" -q "+quality+" "+self.runname+"/consensus.cns "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map 2>"+self.runname+"/assemble.log")
		else:
			#os.system("maq fastq2bfq fastqs/"+self.name+"_1.fastq "+self.runname+"/reads_1.bfq")
			#os.system("maq fastq2bfq fastqs/"+self.name+"_2.fastq "+self.runname+"/reads_2.bfq")
			#os.system("maq map -n "+mismatches+" -a "+maxinsertsize+" -u "+self.runname+"/unmap_1.txt "+self.runname+"/reads_1.map "+self.runname+"/ref.bfa  "+self.runname+"/reads_1.bfq "+self.runname+"/reads_2.bfq ")
			os.system("maq mapcheck -m "+maxsnps+" -q "+quality+" "+self.runname+"/ref.bfa "+pool.runname+"/reads_1.map >"+self.runname+"/mapcheck.txt")
			os.system("maq assemble -p -m "+maxsnps+" -q "+quality+" "+self.runname+"/consensus.cns "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map 2>"+self.runname+"/assemble.log")
		os.system("maq cns2fq "+self.runname+"/consensus.cns >"+self.runname+"/cns.fq")
		os.system("maq cns2snp "+self.runname+"/consensus.cns >"+self.runname+"/cns.snp")
		if pairedend=='n':
			#os.system("maq.pl SNPfilter -a -d "+str(int(mindepth))+" -q "+quality+" "+self.runname+"/cns.snp >"+self.runname+"/cns.final.snp")
			os.system("maq.pl SNPfilter -a -d "+str(int(mindepth))+" -q 20 "+self.runname+"/cns.snp >"+self.runname+"/cns.final.snp")
		else:
			#os.system("maq.pl SNPfilter -d "+str(int(mindepth))+" -q "+quality+" "+self.runname+"/cns.snp >"+self.runname+"/cns.final.snp")
			os.system("maq.pl SNPfilter -d "+str(int(mindepth))+" -q 20 "+self.runname+"/cns.snp >"+self.runname+"/cns.final.snp")
	
		print "\nRunning Maq pileup on "+self.name+'...'
		#Maq pileup commands. change -m 2 -q 30 to change max snps per read and min read quality. Remove -s for paired-end
		if pairedend=='n':
			os.system("maq pileup -s -m "+maxsnps+" -q "+quality+" -v "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map > "+self.runname+"/all.pileup")
		else:
			os.system("maq pileup -p -m "+maxsnps+" -q "+quality+" -v "+self.runname+"/ref.bfa "+self.runname+"/reads_1.map > "+self.runname+"/all.pileup")


	def runVelvet(self, output, mindepth, pairedend, meaninsertsize):
		print "\nRunning velvet on unassembled reads from "+self.name+'...'
		#Creating single unmap file where there are multiple
		os.system("cat "+self.runname+"/unmap*_*.txt > "+self.runname+"/unmap.txt")
		
		#Changing unmap txt file to fastq format
		
		unmapout=open(self.runname+"/unmap.fastq", "w")
		lines=open(self.runname+"/unmap.txt", "rU").readlines()
		for line in lines:
			words=line.split()
			print >> unmapout, '@'+words[0]
			print >> unmapout, words[2]
			print >> unmapout, '+'
			print >> unmapout, words[3]
		unmapout.close()
		
		#Running Velvet.
		if pairedend=='n':
			os.system("bsub /nfs/pathogen/sh16_scripts/velvet_assembly.sh "+self.runname+"/unmap.fastq "+meaninsertsize)
			#os.system("velveth "+self.runname+"/velvet 21 -fastq -short "+self.runname+"/unmap.fastq")
			#os.system("velvetg "+self.runname+"/velvet -cov_cutoff "+mindepth+" -min_contig_lgth 100")
		else:
			os.system("bsub /nfs/pathogen/sh16_scripts/velvet_assembly.sh "+self.runname+"/unmap.fastq "+meaninsertsize)
			#os.system("velveth "+self.runname+"/velvet 21 -fastq -shortPaired "+self.runname+"/unmap.fastq")
			#os.system("velvetg "+self.runname+"/velvet -ins_length "+meaninsertsize+" -cov_cutoff "+mindepth+" -min_contig_lgth 100")
		
		
		
#		print "\nRunning blast on contigs from unassembled reads from "+self.name+'...'
#		
#		#Blasting Velvet contigs against Bacterial databse
#		os.system("blastall -p blastn -o "+self.runname+"/velvet/contigs.blast -d /data/blastdb/Bacteria_DB -v 1 -b 1 -m 8 -i "+self.runname+"/velvet/contigs.fa")
#				
#		#Reading blast output and extracting some hit info using fastacmd
#		lines=open(self.runname+"/velvet/contigs.blast", "rU").readlines()
#		lastcontig=''
#		
#		for line in lines:
#			words=line.strip().split()
#			if words[0]==lastcontig:
#				continue
#			else:
#				lastcontig=words[0]
#			
#			fastacmdout=os.popen("fastacmd -d /data/blastdb/Bacteria_DB -s "+words[1].split('|')[1])
#			fastacmdstr=fastacmdout.readlines()[0]
#			
#			print >> output, self.name+','+words[0]+','+words[1]+','+'_'.join(fastacmdstr[1:].split()[1:]).replace(',','')+','+','.join(words[2:])




########
# Main #
########


if __name__ == "__main__":
        argv=sys.argv[1:]
        ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet, mindepth, maxsnps, pairedend, maxinsertsize, mismatches=getOptions(argv)		

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
		
		if not os.path.isfile(pool) and not os.path.islink(pool):
			print "File "+pool+" not found! Skipping..."
			continue
	
		if pool[-1]=='/':
			pool=pool[:-1]
		pool=pool.split('/')[-1].split('.')[0]
		if pairedend=='y':
			if pool.split('_')[-1]=='1':
				pool='_'.join(pool.split('_')[:-1])
			else:
				continue
		name=pool
		pool=pool+'_maq'
		print pool+'...',
		sys.stdout.flush()
		if not os.path.isdir(pool):
		 	print "pool "+pool+" not found! Creating...",
			os.system("mkdir "+pool)
			runmaq='y'
	
		userinput='x'
		if runmaq=='n':
			if not os.path.isfile(pool+'/reads_1.map') and not os.path.isfile(pool+'/reads_1.map.gz'):
				print "\nError: No Maq map ("+pool+"/reads_1.map) found!"
				while userinput not in ['i','c','a']:
					userinput=raw_input("Would you like to ignore current folder, run maq on current folder or run maq on all folders? (i/c/a): ")
					userinput.lower()
				if userinput=='a':
					runmaq='y'
	
			elif not os.path.isfile(pool+'/all.pileup') and not os.path.isfile(pool+'/all.pileup.gz'):
				if runpileups=='n' and runmaq=='n':
					print "\nError: No pileup file ("+pool+"/all.pileup) found!"
					while userinput not in ['y','n','a']:
						userinput=raw_input("Would you like to create it on the fly? (y/n/a): ")
						userinput.lower()
					if userinput=='a':
						runpileups='y'
						userinput='y'
				else:
					userinput=='y'
			else:
				if os.path.isfile(pool+'/all.pileup'):
					templine=open(pool+"/all.pileup", "rU").readlines(1)[0].strip().split()
				else:
					templine=gzip.open(pool+"/all.pileup.gz", "r").readlines(1)[0].strip().split()
					
				if len(templine)!=7:
					if runpileups=='n' and runmaq=='n':
						print "\nError: Pileup file ("+pool+"/all.pileup) of wrong format!"
						while userinput not in ['y','n','a']:
							userinput=raw_input("Would you like to recreate it on the fly? (y/n/a): ")
							userinput.lower()
						if userinput=='a':
							runpileups='y'
							userinput='y'
					else:
						userinput='y'
				
		if runmaq=='y' or userinput in['c','y','x']:
			pools.append(SNPanalysis())
			pools[count].runname=pool
			pools[count].name=name
			poolsort.append(pool)
			if userinput=='c':
				pools[count].runmaq='m'
			elif userinput=='y':
				pools[count].runmaq='p'
			print 'ok'
			sys.stdout.flush()
			count=count+1
			
		else:
			print 'excluded'
			sys.stdout.flush()
		
		
	if len(pools)==0:
		print "\nError: No valid input folders!"
		sys.exit()
	
	if snptype=='ssaha':
		mindepth=float(mindepth)/5
	
	
	
	

	
	#Loading reference sequence
	if not os.path.isfile(ref) and not os.path.islink(ref):
		print "\nError: Cannot open reference fasta file "+ref
		sys.exit()
	if ref.split('.')[-1]=="gz":
		reflines=gzip.open(ref, "r").readlines()
	else:
		reflines=open(ref, "rU").readlines()
		
	refseq=''
	contigstarts={}
	for line in reflines:
		if len(line.strip())>0 and line.strip()[0]=='>':
			contigstarts[line.strip().split()[0][1:]]=len(refseq)
		elif len(line.strip())>0 and line.strip()[0]!='>':
			refseq=refseq+line.strip().upper()
	
	
	
	
	
	#Running Maq and Velvet where required
	
	if velvet=='y':
		output=open(outfile+'_unmapped_contig_hits.csv','w')
	
	for pool in pools:
		if runmaq=='y' or pool.runmaq=='m':
			pool.runMaq(ref, str(quality), str(mindepth), str(maxsnps), pairedend, str(maxinsertsize), str(mismatches))
			
		if velvet=='y':
			pool.runVelvet(output, str(mindepth), pairedend, str(maxinsertsize))
	
	if velvet=='y':
		output.close()
	if os.path.isfile(pools[0].runname+'/assemble.log'):
		lenstring=os.popen('grep "reference length:" '+pools[0].runname+'/assemble.log').read()
	elif os.path.isfile(pools[0].runname+'/assemble.log.gz'):
		os.system("gunzip "+pools[0].runname+'/assemble.log.gz')
		lenstring=os.popen('grep "reference length:" '+pools[0].runname+'/assemble.log').read()
		os.system("gzip "+pools[0].runname+'/assemble.log')
	else:
		print "Cannot find assembly log file"
		sys.exit()
		
	reflen=len(refseq)
#	if int(lenstring.strip().split(":")[1])!=reflen:
#		print "Error: reference sequence length is not equal to length in "+pools[0].runname+"/assemble.log"
#		sys.exit()
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
		for i in range(reflennogaps):
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
				while y<len(lines) and foundpool=='n' and lines[y].split()[1] not in ('CDS', 'rRNA', 'tRNA', 'gene', 'misc_feature', 'sig_peptide', 'repeat_unit', 'repeat_region', 'misc_RNA'):
					if lines[y].split()[1][0]!='/' and y==x+extralocationlines:
						location=location+lines[y].split()[1].strip()
						extralocationlines=extralocationlines+1
					if '/colour=11' in lines[y].split()[1]:
						pseudo='y'
					if '/systematic_id' in lines[y].split()[1] or '/locus_tag' in lines[y].split()[1] or '/gene' in lines[y].split()[1]:
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
	print reflen
	
	print "\nIdentifying SNPs and % of reference mapped..."
	sys.stdout.flush()
	for pool in pools:
		print pool.name+'...',
		sys.stdout.flush()
		#lines=open(pool.name+'/cns.final.snp', 'rU').readlines()
		
		if os.path.isfile(pool.runname+'/cns.final.snp'):
			lines=open(pool.runname+'/cns.final.snp', 'rU').readlines()
		elif os.path.isfile(pool.runname+'/cns.final.snp.gz'):
			lines=gzip.open(pool.runname+'/cns.final.snp.gz', 'r').readlines()
		poolcovout=open(pool.runname+'/'+pool.name+'_cov.txt', 'w')
		poolcovynout=open(pool.runname+'/'+pool.name+'_covyn.txt', 'w')
		
		snpcount=0
		tempsnps={}
		tempseq=['-']*reflen
		tempsnpsummary={'A':{'C':0, 'G':0, 'T':0}, 'C':{'A':0, 'G':0, 'T':0}, 'G':{'C':0, 'A':0, 'T':0}, 'T':{'C':0, 'G':0, 'A':0}}
		hetsnpcount=0
		for line in lines:
			words=line.split()
			contig=words[0]
			snplocation=int(words[1])+contigstarts[contig]
			refbase=words[2]
			cnsbase=words[3]
			if int(words[5])<mindepth:
				continue
			if cnsbase not in ['A','C','G','T']:
				hetsnpcount=hetsnpcount+1
				#continue#Change if you want heterozygous SNP calls
				cnsbase=words[9]
				if cnsbase not in ['A','C','G','T'] or int(words[4]<20) or cnsbase==refbase or words[11]==refbase:# or int(words[8])<quality
					continue
			
			
			if allsnps.has_key(snplocation):
				allsnps[snplocation][0]=allsnps[snplocation][0]+1
				allsnps[snplocation][2][cnsbase]=allsnps[snplocation][2][cnsbase]+1
			else:
				allsnps[snplocation]=[1,refbase,{'A':0, 'C':0, 'G':0, 'T':0}]
				allsnps[snplocation][2][cnsbase]=1
			tempsnpsummary[refbase][cnsbase]=tempsnpsummary[refbase][cnsbase]+1
			tempseq[snplocation-1]=cnsbase
			snpcount=snpcount+1
		
		print '%d SNPs found...' % (snpcount),
		print hetsnpcount,
		sys.stdout.flush()
		if os.path.isfile(pool.runname+'/all.pileup'):
			lines=open(pool.runname+'/all.pileup', 'rU').readlines()
		elif os.path.isfile(pool.runname+'/all.pileup.gz'):
			lines=gzip.open(pool.runname+'/all.pileup.gz', 'r').readlines()
		mapped=0.0
		
		#for x in hassnps.keys():
			#line=lines[x-1]
		hetcount=0
		for line in lines:
			
			words=line.split()
			snplocation=int(words[1])
			refbase=words[2]
			readdepth=int(words[3])
			reads=words[4][1:].upper()
			
			avcoverage[snplocation]=avcoverage[snplocation]+readdepth
	
			if readdepth>=mindepth:
				if tempseq[snplocation-1]=='-':
					tempseq[snplocation-1]=refbase
				covcount[snplocation]=covcount[snplocation]+1
				mapped=mapped+1
				
				
				reads=reads.replace(',','.')
				curread=reads[0]
				reads=reads.replace(reads[0],'')
				if len(reads)>0:
					hetcount=hetcount+1
				
				print >> poolcovynout, 1
			else:
				print >> poolcovynout, 0
			print >> poolcovout, readdepth
			
		
		poolcovout.close()
		poolcovynout.close()
		pool.SNPs=snpcount
		pool.sequence=''.join(tempseq)
		#pool.mapped=tempsnps
		pool.snpsummary=tempsnpsummary
		pool.nummapped=mapped
		pool.hetsnpcount=hetsnpcount
		pool.hetcount=hetcount
		pool.percentmapped=((mapped/reflen)*100)
		print '%.2f%% of reference mapped' % (pool.percentmapped),
		print hetcount
		sys.stdout.flush()
	
	
	
	
	#
	
	#testout=open("testing.out","w")
	#
	#print >> testout, ">ref\n"+refseq
	#
	#for strain in pools:
	#	print >> testout, ">"+strain.name+"\n"+strain.sequence
	#
	#testout.close()
	
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
		#dndsout=open(outfile+'_dnds.out','w')
		
		#dnbydsstats={}
		
		
	
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
		for pool in pools:
			if pool.name==ref:
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
			
			
			pool.dNdSstats, pool.SNPtypes, pool.AAtypes=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers)
			
			
			#N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd
			if pool.dNdSstats['dS']!=0:
				print "%.2f" % (pool.dNdSstats['dN']/pool.dNdSstats['dS'])
				#print >> dndsout, '\t'+str(pool.dNdSstats['dN']/pool.dNdSstats['dS']),
				sys.stdout.flush()
			else:
				print "-"
			
			
			
		#Testing doing each gene seperately
		
		#for sequence in sequences.keys():
		#		pool.sequence=''
		#		for i, j in enumerate(sequences[ref]):
		#			if j!='-':
		#				pool.sequence=pool.sequence+sequences[sequence][i]
	#	print "\nCalculating CDS dN/dS values...",
	#	for key in CDSposns.keys():
	#	
	#		refCDSseq=''
	#		
	#		for x in CDSposns[key]:
	#			if x[1]=='f':
	#				refCDSseq=refCDSseq+refseq[x[0]-1].upper()
	#			else:
	#				refCDSseq=refCDSseq+comp[refseq[x[0]-1].upper()]
	#		
	#		print >> dndsout, '\n'+CDSnames[key],
	#		#print CDSnames[key]
	#		sys.stdout.flush()
	#		
	#		for pool in pools:
	#			if pool.name==ref:
	#				continue
	#				
	#			tmpCDSseq=''
	#			
	#			for y in CDSposns[key]:
	#				if y[1]=='f' and pool.sequence[y[0]-1]!='-':
	#					tmpCDSseq=tmpCDSseq+pool.sequence[y[0]-1].upper()
	#				elif pool.sequence[y[0]-1]!='-':
	#					tmpCDSseq=tmpCDSseq+comp[pool.sequence[y[0]-1].upper()]
	#				else:
	#					tmpCDSseq=tmpCDSseq+'-'
	#			#print refCDSseq+'\n'+tmpCDSseq+'\n'
	#			dnbydsstats, gene_snptypes=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers)
	#			
	#			if dnbydsstats['dS']!=0:
	#				print >> dndsout, '\t'+str(dnbydsstats['dN']/dnbydsstats['dS']),
	#			else:
	#				print >> dndsout, '\t-',
	#			dndsout.flush()
	#			
	#	
	#	
	#	dndsout.close()
			
			
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
	
	if align=='y':
		aln={}
		alnfile=open(outfile+'.aln', 'w')
	if tabfile=='y':
	        tabout=open(outfile+'.tab','w')
	if graphs=='y':
		snpgraph=open(outfile+'_snps.txt','w')
		covgraph=open(outfile+'_cov.txt','w')
		covcountgraph=open(outfile+'_covcount.txt','w')
	
	#minquality=10
	#minqualityperdepth=4
	
	snpsequence={}
	for pool in pools:
		snpsequence[ref]=''
		snpsequence[pool.name]=''
	
	if ref=='':
		ref=pools[0].name
	
	if tabfile=='y':
		print >> tabout, 'ID   SNP'
	
	poscount=1
	
	
	
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
					if snptype=='-':
						snptype=pool.SNPtypes[j]
					elif pool.SNPtypes[j] not in snptype:
						snptype=snptype+'/'+pool.SNPtypes[j]
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
						tabstring=tabstring+pool.name+'='+pool.sequence[j]+' '
						if pool.SNPtypes.has_key(j):
							if pool.SNPtypes[j]=='S':
								tabstring=tabstring+'(synonymous) '
							elif pool.SNPtypes[j]=='N':
								if pool.AAtypes.has_key(j):
									tabstring=tabstring+'(non-synonymous, AA='+pool.AAtypes[j]+') '
								else:
									tabstring=tabstring+'(non-synonymous, codon='+pool.sequence[j-1]+pool.sequence[j]+pool.sequence[j+1]+') '
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
		count=0
		inaln=[]
		for pool in snpsequence.keys():
			if len(snpsequence[pool].replace('-',''))>(len(snpsequence[ref])/2):
				inaln.append(pool)
				count=count+1
		
		
		
		print >> alnfile, count, len(snpsequence[ref])
		#print >> alnfile, ref+' '*(10-len(ref))+refseq
		for pool in inaln:
			print >> alnfile, pool.replace('pool','').split('.')[0]+' '*(10-len(pool.replace('pool','').split('.')[0]))+snpsequence[pool]
	
					
		alnfile.close()
	
	
	summaryfile=open(outfile+'_summary.out','w')
	
	print >> summaryfile, 'Strain\tPercent of Reference Mapped\tdN/dS',
	
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
				print >> summaryfile, pool.name+'\t'+str(pool.percentmapped)+'\t'+str(pool.dNdSstats['dN']/pool.dNdSstats['dS']),
			else:
				print >> summaryfile, pool.name+'\t'+str(pool.percentmapped)+'\t-',
		else:
			print >> summaryfile, pool.name+'\t'+str(pool.percentmapped)+'\t-',
		for x in ['A', 'C', 'G', 'T']:
			for y in ['A', 'C', 'G', 'T']:
				if x!=y:
					print >> summaryfile, '\t'+str(pool.snpsummary[x][y]),
					count=count+pool.snpsummary[x][y]
		
		if embl!='':
			for CDS in CDSposns.keys():
				for y in CDSposns[CDS]:
					if pool.sequence[y[0]]!=refseq[y[0]] and pool.sequence[y[0]]!='-':
						intragenic=intragenic+1
					if pool.sequence[y[0]]!='-':
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
	
	
	
	if raxml=='y':
		RAxMLphylogeny(outfile.split('/')[-1]+".aln", model, bootstrap)
	
		
		
	if tabfile=='y':
		tabout.close()
	
	if graphs=='y':
		snpgraph.close()
		covgraph.close()
		covcountgraph.close()
	output.close()
	summaryfile.close()
	
	
	print "All finished.\n"
	
	
