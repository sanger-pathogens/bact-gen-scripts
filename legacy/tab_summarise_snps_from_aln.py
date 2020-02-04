#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math
from scipy.stats import chi2



def Usage():
	print 'summarise_snps_from_aln.py Usage:'
	print 'summarise_snps_from_aln.py [options] <input alignment file>'
	print 'Options:'
	print '-r\t\tname of reference strain'
	print '-e <filename>\tembl file of reference strain [optional]'
	print '-x\t\tproduce moving window snp plot and chi-squared test plot'
	print '-c\t\tproduce recombination hotspot tab files for each sequence'
	print '-o\t\tprefix for output file names'
	print '-t\t\tproduce tab files of snp locations'
	print '-a\t\tproduce alignment of snp locations'
	print '-g\t\tproduce mapping graph files'
	print '-p\t\trun phylogeny of snp locations using RAxML'
	print '-m <model>\tmodel of evolution for phylogeny [GTRGAMMA/GTRGAMMAI/GTRCAT]'
	print '-b <int>\tnumber of bootstrap replicates [0 = do not run bootstrap]'
	print '-h\t\tshow this help'
	print 'Written by Simon R Harris, Wellcome Trust Sanger Institute, UK. 2009'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "he:r:o:tapgb:m:xcT:k", ["model=", "bootstrap=", "help", "phylogeny" "embl=" "ref=", "out=", "align", "tabfile", "graphs", "recombination", "chisquare", "keep", "tabin="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	ref=''
	outfile=''
	inputfile=[]
	i=0
	tabfile='n'
	align='n'
	embl=''
	dirty='n'
	graphs='n'
	raxml='n'
	bootstrap=100
	model='GTRGAMMA'
	chisquare='n'
	recomb='n'
	tabin=""
	keep="n"
	
	
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
		elif opt in ("-p", "--phylogeny"):
			raxml='y'
		elif opt in ("-a", "--align"):
			align='y'
		elif opt in ("-x", "--chisquare"):
			chisquare='y'
		elif opt in ("-c", "--recombination"):
			recomb='y'
		elif opt in ("-g", "--graphs"):
			graphs='y'
		elif opt in ("-e", "--embl"):
			embl=arg
		elif opt in ("-b", "--bootstrap"):
			bootstrap=int(arg)
		elif opt in ("-m", "--model"):
			model=arg
		elif opt in ("-T", "--tabin"):
			tabin=arg
		elif opt in ("-k", "--keep"):
			keep="y"
	
	inputfile=args
	
	
	if inputfile==[]:
		print 'Error: No input file selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=ref.split('.')[0]
	if align=='n' and raxml=='y':
		align='y'

	return ref, inputfile[0], outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, chisquare, recomb, tabin, keep


def countcodonchanges(codon, SNPcodon, geneticcode, sd, nd, loopsd=0, loopnd=0, pathcount=0):
	
	for x in range(3):
		if codon[x]!=SNPcodon[x]:
			newSNPcodon=SNPcodon[:x]+codon[x]+SNPcodon[x+1:]
			
			#print  SNPcodon, newSNPcodon, geneticcode[SNPcodon], geneticcode[newSNPcodon]
			
			if geneticcode[newSNPcodon]=='STOP':
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
	
	
	
def getCDSseq(embldata, CDSstring, direction, CDScount, feature):
	
	if '(' in CDSstring:
		if CDSstring.split('(')[0]=='complement':
			if direction=='f':
				direction='r'
			else:
				direction='f'
		parts='('.join(CDSstring.split('(')[1:])
		parts=')'.join(parts.split(')')[:-1])
		getCDSseq(embldata, parts, direction, CDScount, feature)
		
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
				
				embldata[donelen+x]=words[1][0]
				if words[1]=='CDS':
					CDSposns[CDScount].append((x, direction))


	return embldata, CDSposns






class SNPanalysis:
	def __init__(self, fastq='', directory='', mapped={}, runmaq='n', CDSseq=''):
		self.fastq=fastq
		self.directory=directory
		self.mapped=mapped
		self.runmaq=runmaq
		self.CDSseq=CDSseq
		self.N=0.0
		self.S=0.0
		self.dN=0.0
		self.dS=0.0
		self.pN=0.0
		self.pS=0.0
		self.varianceN=0.0
		self.varianceS=0.0
		self.z=0.0
		self.Nd=0.0
		self.Ns=0.0
		self.goodlen=0
		self.CDSseq=''
		self.snpsummary={}
		self.nummapped=0
		self.percentmapped=0.0
		self.intragenic=0
		self.intergenic=0
		self.SNPs={}



if __name__ == "__main__":
        argv=sys.argv[1:]
        ref, inputfile, outfile, tabfile, align, embl, raxml, graphs, bootstrap, model, chisquare, recomb, tabin, keep=getOptions(argv)
		
	
	snps={}
	
	refbases={}
	bases={}#0=A 1=C 2=G 3=T
	nostates=0
	converter={}
	convertback={}
	snpbases={}
	
	snpstructs=[]
	count=0
	namesort=[]
	runpileups='n'
	
	print '\nReading input alignment...',
	sys.stdout.flush()
	
	sequences={}
	#for x in hassnps.keys():
		#line=lines[x-1]


	currseq=''

	if os.path.getsize(inputfile)<2000000000:
		lines=open(inputfile, "rU").read().split('>')[1:]
	else:
		lines=[]
		count=-1
		for linea in open(inputfile, "rU"):
			if linea[0]==">":
				count=count+1
				lines.append(linea.split()[0][1:]+'\n')
			else:	
				lines[count]=lines[count]+linea
		linesa=[]
		sequences={}
	
	
	for line in lines:
		words=line.strip().split('\n')
		sequences[words[0]]=''.join(words[1:])
	snpstructs.append(SNPanalysis())

	#print sequences.keys()

#	lines=open(inputfile, 'rU').readlines()
#	for line in lines:
#		line=line.strip()
#		if len(line)>0 and line[0]=='>':
#			currseq=line.split()[0][1:]
#			sequences[currseq]=''
#			snpstructs.append(SNPanalysis())
#		elif len(line)>0 and currseq!='':
#			sequences[currseq]=sequences[currseq]+line.upper().replace('N','-')

		

	reflen=len(sequences[sequences.keys()[0]])

	for sequence in sequences.keys():
                if len(sequences[sequence])!=reflen:
                        print "\nERROR!: sequences are not all of the same length!!!\n"
                        sys.exit()
                sequences[sequence]=sequences[sequence].upper().replace('N','-')

	
	#print sequences.keys()
	print "Found", len(sequences.keys()), "sequences of length", reflen
	sys.stdout.flush()
	
	if not ref in sequences.keys():
		print "Error!!!: Reference ("+ref+") not in alignment"
		sys.exit()

	

	if tabin!="":
		print '\nReading tab file...',
		sys.stdout.flush()
		regions=[]
		
		for line in open(tabin, 'rU'):
			if len(line.split())>2 and line.split()[0]=="FT" and line.split()[1]=="misc_feature":
				line=line.replace(")","").replace("complement(", "")
				regions.append([int(line.split()[2].split('..')[0]),int(line.split()[2].split('..')[1])])
		
		
		
		if keep=='y':
			regionlen=0
		else:
			regionlen=reflen
		for region in regions:
			if keep=='y':
				regionlen=regionlen+(region[1]-(region[0]-1))
			else:
				regionlen=regionlen-(region[1]-(region[0]-1))
		if keep=='y':
			print "Found", len(regions), "regions totalling", regionlen, "bases to keep"
		else:
			print "Found", len(regions), "regions totalling", reflen-regionlen, "bases to remove.", regionlen, "bases remaining"
		
		for i, sequence in enumerate(sequences.keys()):
			if sequence==ref:
				continue
			newsequence=""
			for x, region in enumerate(regions):
				if keep=='y':
					if x==0:
						newsequence=newsequence+"-"*(region[0]-1)+sequences[sequence][region[0]-1:region[1]]
					else:
						newsequence=newsequence+"-"*(region[0]-(regions[x-1][1]+1))+sequences[sequence][region[0]-1:region[1]]
					if x==len(regions)-1:
						newsequence=newsequence+"-"*(len(sequences[sequence])-region[1])
				else:
					if x==0:
						newsequence=newsequence+sequences[sequence][:region[0]-1]+"-"*(region[1]-(region[0]-1))
					else:
						newsequence=newsequence+sequences[sequence][regions[x-1][1]:region[0]-1]+"-"*(region[1]-(region[0]-1))
					if x==len(regions)-1:
						newsequence=newsequence+sequences[sequence][region[1]:]
			sequences[sequence]=newsequence
	
	else:
		regions=[[1,reflen-1]]
		keep="y"
		
	count=0
	if ref!='':
		refconvert={}
		refconvertb={}
		#print reflen
		reflennogaps=len(sequences[ref].replace('-',''))
		for x, y in enumerate(sequences[ref]):
			if y!='-':
				
				refconvert[x]=count#use refconvert to convert alignment positions to positions in ref sequence
				refconvertb[count]=x#opposite to refconvert
				count=count+1
	
	snplocations=[]
	snpbases=[]
	#basecounts=[]
	print "\nIdentifying SNPs...",
	sys.stdout.flush()
	if keep=="y":
		for region in regions:
			for x in range(region[0]-1,region[1]+1):
				numbases=0
				foundbases={}	
				for key in sequences.keys():
					base=sequences[key][x].upper()
					if base not in foundbases.keys() and base!='-':
						foundbases[base]=1
						numbases=numbases+1
					elif base!='-':
						foundbases[base]=foundbases[base]+1
				if numbases>1:
					snplocations.append(x)
					snpbases.append(foundbases)
					#basecounts.append(numbases)
	
	
	else:
		for x in range(reflen):
			numbases=0
			foundbases={}	
			for key in sequences.keys():
				base=sequences[key][x].upper()
				if base not in foundbases.keys() and base!='-':
					foundbases[base]=1
					numbases=numbases+1
				elif base!='-':
					foundbases[base]=foundbases[base]+1
			if numbases>1:
				snplocations.append(x)
				snpbases.append(foundbases)
				#basecounts.append(numbases)
	
	
	print "Found", len(snplocations), "sites with a SNP"
	sys.stdout.flush()
	
	
	if embl!='':
		print "\nReading EMBL file(s)...",
		sys.stdout.flush()
		CDSposns={}
		CDSnames={}
		embldata=[]
		embls=embl.split(',')
		donelen=0
		CDScount=0		
		for i in range(reflennogaps):
			embldata.append('I')
		for em in embls:
			lines=open(em,'rU').readlines()
			chromosomelen=0
			for x, line in enumerate(lines):
				words=line.strip().split()
				if len(words)>2 and words[1] in ('CDS', 'rRNA', 'tRNA') and '..' in words[2]:
					#print words
					foundname='n'
					y=x+1
					locustag=''
					location=words[2].replace('<','').replace('>','')
					extralocationlines=1
					while foundname=='n' and lines[y].split()[1] not in ('CDS', 'rRNA', 'tRNA', 'gene', 'misc_feature', 'sig_peptide', 'repeat_unit', 'repeat_region', 'misc_RNA'):
						if lines[y].split()[1][0]!='/' and y==x+extralocationlines:
							location=location+lines[y].split()[1].strip()
							extralocationlines=extralocationlines+1
						if '/systematic_id' in lines[y].split()[1] or '/locus_tag' in lines[y].split()[1]:
							locustag=lines[y].split()[1].split('"')[1]
							#print locustag
							foundname='y'
						y=y+1
					
					direction='f'
					feature=words[1]
					if feature=='CDS':
						CDScount=CDScount+1
						CDSposns[CDScount]=[]
						CDSnames[CDScount]=locustag
					embldata, CDSposns=getCDSseq(embldata, location, direction, CDScount, feature)
					
				donelen=donelen+chromosomelen
		
		snptoCDSname={}
		
		for i in CDSposns.keys():
			inregion="n"
			if keep=="y":
				for region in regions:
					if (CDSposns[i][0][0]>region[0] and CDSposns[i][0][0]<region[1]) or (CDSposns[i][-1][0]>region[0] and CDSposns[i][-1][0]<region[1]) or (CDSposns[i][0][0]<region[0] and CDSposns[i][-1][0]>region[1]):
						inregion="y"
						break
			elif keep=="n":
				for region in regions:
					if CDSposns[i][0][0]>region[0] and CDSposns[i][-1][0]<region[1]:
						inregion="y"
						break

			if (keep=="y" and inregion=="n") or (keep=="n" and inregion=="y"):
				del CDSposns[i]
		
		print "Found", len(CDSposns.keys()), "CDSs"
		
		for i in CDSposns.keys():
			for j in CDSposns[i]:
				snptoCDSname[refconvertb[j[0]]]=CDSnames[i]
			
		print "Done"
		sys.stdout.flush()

	
	#print CDSposns
	#sys.exit()
	namesort.sort()
	
	#avcoverage={}
	#covcount={}
	
	allsnps={}
	refbases={}
	
	alssnpsummary={}
	
				
	
	tempsnps={}
	tempsnpsummary={'A':{'C':0, 'G':0, 'T':0}, 'C':{'A':0, 'G':0, 'T':0}, 'G':{'C':0, 'A':0, 'T':0}, 'T':{'C':0, 'G':0, 'A':0}}
	
	olddone=-1
	if chisquare=='y':
	
		if not os.path.isdir("chisquared_plots"):
			os.system("mkdir chisquared_plots")
	
		print "\nCalculating chi-squared for snps..."
		sys.stdout.flush()
		
		#Chi-squared looking for snp clustering. Window size=1000
		window=100
		windowcounts=[0]*(reflennogaps)
		X2=[]
		pvalues=[]
		currwindowcount=0
		currposn=0
		X2converter={}
		
		
		expectedsnps=float(len(snplocations))/reflennogaps
		expectedsnps=expectedsnps*window
		expectednonsnps=window-expectedsnps
		
		while (expectedsnps<10 or expectednonsnps<10) and window<(reflen/1000):
			print "Expected frequencies too low for chi-squared test using a window size of "+str(window)+" (",expectedsnps,"and", expectednonsnps,"). Increasing window size to "+str(window*10)
			expectedsnps=expectedsnps*10
			expectednonsnps=expectednonsnps*10
			window=window*10
	
		
		print "\nCalculating moving window snp counts..."
		sys.stdout.flush()
		
		for x in snplocations:
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
	
		print "\nCalculating chi-squared values..."
		sys.stdout.flush()	
				
		for x in windowcounts:
			if not X2converter.has_key(x):
				X2converter[x]=[ float(((abs(x-expectedsnps)-0.5)**2)/expectedsnps) + float(((abs((window-x)-expectednonsnps)-0.5)**2)/expectednonsnps),0]
				if x<expectedsnps:
					X2converter[x][1]=1
				else:
					X2converter[x][1]=(1-chi2.cdf(X2converter[x][0],1))
			
			X2.append(X2converter[x][0])
			pvalues.append(X2converter[x][1])
			
			
	
		print "\nPrinting chi-squared files..."
		sys.stdout.flush()
		
		output=open("chisquared_plots/"+outfile+'_chi_squared_pvalue.plot','w')
		outputb=open("chisquared_plots/"+outfile+'_chi_squared.plot','w')
		outputc=open("chisquared_plots/"+outfile+'_moving_window_snps'+str(window)+'.plot','w')
		
		for x, pvalue in enumerate(pvalues):
			print >> output, pvalue
			print >> outputb, X2[x]
			print >> outputc, windowcounts[x]
			
		output.close()
		outputb.close()
		outputc.close()
		
		print "Done."
		
		
		keysort=sequences.keys()
		keysort.sort()
		
		newaln=open('chisquared_plots/'+outfile+'_norecomb.aln', 'w')
		print >> newaln, ">"+ref
		print >> newaln, sequences[ref]
		for strainnum, strain in enumerate(keysort):
			if strain==ref:
				continue
			print "Calculating chi-squared for "+strain+"...",
			sys.stdout.flush()
				
			#Chi-squared looking for snp clustering. Window size=1000
			window=100
			windowcounts=[0]*(reflennogaps)
			X2=[0]*(reflennogaps)
			pvalues=[1]*(reflennogaps)
			currwindowcount=0
			currposn=0
			X2converter={}
			
			snplocationstemp=[]
			for x in snplocations:
				if sequences[ref][x]!=sequences[strain][x] and sequences[strain][x]!='-' and sequences[ref][x]!='-':
					snplocationstemp.append(x)
			
			expectedsnps=float(len(snplocationstemp))/reflennogaps
			
			expectedsnps=expectedsnps*window
			expectednonsnps=window-expectedsnps
			
			while (expectedsnps<10 or expectednonsnps<10) and window<(reflen/100):
				expectedsnps=expectedsnps*10
				expectednonsnps=expectednonsnps*10
				window=window*10
				
			
			print window, expectedsnps, expectednonsnps,
			
			for x in snplocationstemp:
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
		
					
			for x in windowcounts:
				if not X2converter.has_key(x):
					if expectedsnps>=10 and expectednonsnps>=10:
						X2converter[x]=[ float(((abs(x-expectedsnps)-0.5)**2)/expectedsnps) + float(((abs((window-x)-expectednonsnps)-0.5)**2)/expectednonsnps),0]
					else:
						X2converter[x]=[0,1]
					if x<expectedsnps or expectedsnps<10 or expectednonsnps<10:
						X2converter[x][1]=1
					else:
						X2converter[x][1]=(1-chi2.cdf(X2converter[x][0],1))
				
				X2.append(X2converter[x][0])
				pvalues.append(X2converter[x][1])
				
			
			output=open("chisquared_plots/"+strain+'_chi_squared_pvalue.plot','w')
			outputb=open("chisquared_plots/"+strain+'_chi_squared.plot','w')
			outputc=open("chisquared_plots/"+strain+'_moving_window_snps'+str(window)+'.plot','w')
			
			newseq=sequences[strain]
			for x, pvalue in enumerate(pvalues):
				print >> output, pvalue
				print >> outputb, X2[x]
				print >> outputc, windowcounts[x]
				
			count=0
			for x in snplocations:
				if pvalues[x]<0.05:
					newseq=newseq[:x]+'-'+newseq[x+1:]
					count=count+1
				
			output.close()
			outputb.close()
			outputc.close()
			
			print >> newaln, ">"+strain
			print >> newaln, newseq
			
			print count, "bases removed as recombinations"
			
			print "Done."
		print "Done."
			
		newaln.close()
			
	
	
	if recomb=='y':
		print "\nCalculating moving window closest pairs...\n"
		sys.stdout.flush()
		currposn=0
		step=1
		window=10
		closest=[]
		pairwisesnps={}
		
		
		keysort=sequences.keys()
		keysort.sort()
	
		for strainnum, strain in enumerate(keysort):
			if not pairwisesnps.has_key(strain):
				pairwisesnps[strain]={}
			for strainb in keysort[strainnum+1:]:
				if not pairwisesnps.has_key(strainb):
					pairwisesnps[strainb]={}
				if not pairwisesnps[strain].has_key(strainb):
					pairwisesnps[strain][strainb]=[0]*len(snplocations)
				if not pairwisesnps[strainb].has_key(strain):
					pairwisesnps[strainb][strain]=[0]*len(snplocations)
				print '\r',strain, "vs", strainb, "        ",
				sys.stdout.flush()
				if strain==strainb:
					continue
				for base in range(0,len(snplocations),int(step)):
					start=base-(window/2)
					end=base+(window/2)
					if start>=0 and end<len(snplocations):
						
						
						for site in snplocations[start:end]:
							if sequences[strain][site]!=sequences[strainb][site] and sequences[strain][site]!='-' and sequences[strainb][site]!='-':
								pairwisesnps[strain][strainb][base]=pairwisesnps[strain][strainb][base]+1
								pairwisesnps[strainb][strain][base]=pairwisesnps[strainb][strain][base]+1
									
					elif start<0 and end<len(snplocations):
						
						for site in snplocations[len(snplocations)+start:]:
							if sequences[strain][site]!=sequences[strainb][site] and sequences[strain][site]!='-' and sequences[strainb][site]!='-':
								pairwisesnps[strain][strainb][base]=pairwisesnps[strain][strainb][base]+1
								pairwisesnps[strainb][strain][base]=pairwisesnps[strainb][strain][base]+1
						
						for site in snplocations[:end]:
							if sequences[strain][site]!=sequences[strainb][site] and sequences[strain][site]!='-' and sequences[strainb][site]!='-':
								pairwisesnps[strain][strainb][base]=pairwisesnps[strain][strainb][base]+1
								pairwisesnps[strainb][strain][base]=pairwisesnps[strainb][strain][base]+1
						
						
					elif start>=0 and end>=len(snplocations):
						
						for site in snplocations[start:]:
							if sequences[strain][site]!=sequences[strainb][site] and sequences[strain][site]!='-' and sequences[strainb][site]!='-':
								pairwisesnps[strain][strainb][base]=pairwisesnps[strain][strainb][base]+1
								pairwisesnps[strainb][strain][base]=pairwisesnps[strainb][strain][base]+1
						
						for site in snplocations[:end-len(snplocations)]:
							if sequences[strain][site]!=sequences[strainb][site] and sequences[strain][site]!='-' and sequences[strainb][site]!='-':
								pairwisesnps[strain][strainb][base]=pairwisesnps[strain][strainb][base]+1
								pairwisesnps[strainb][strain][base]=pairwisesnps[strainb][strain][base]+1
						
						
						
					else:
						print "Error! Your window size is larger than your sequence!"
						sys.exit()
			
			
		print "\rPrinting pairwise closest tab files..."
		sys.stdout.flush()
		
		
		strainblocks={}
		strainclosest={}
		
		for strain in keysort:
			
			strainblocks[strain]=[]
			
			pairwiseclosest=''
			sums=[]
			for strainb in keysort:
				if strain==strainb:
					continue
				sums.append([0,strainb])
				for base in pairwisesnps[strain][strainb]:
					sums[-1][0]=sums[-1][0]+base
			
			
			sums.sort()
			keysortb=[]
			for sum in sums:
			 keysortb.append(sum[1])
			pairwiseclosest=keysortb[0]
			strainclosest[strain]=pairwiseclosest
			

			lastclosest=''
			start=0
			snpsum=0
			snpsumclosest=0
			for base in range(0,len(snplocations),int(step)):
				closest=''
				numsnps=1000000000
				for strainb in keysortb:
					if strain==strainb:
						continue
					if strainb==pairwiseclosest:	
						numclosestsnps=pairwisesnps[strain][strainb][base]
					if pairwisesnps[strain][strainb][base]<numsnps:
						closest=strainb
						numsnps=pairwisesnps[strain][strainb][base]
					if numsnps==0:
						break
				
				if closest!=lastclosest and lastclosest!='':
					
					strainblocks[strain].append([start, base, lastclosest])
					lastclosest=closest
					start=base
				elif closest!=lastclosest:
					lastclosest=closest
					start=base
			strainblocks[strain].append([start, base, lastclosest])

		
		newstrainblocks={}
		
		for strain in strainblocks.keys():
			newstrainblocks[strain]=[]
			numblocks=0
			other=0
			total=0
			closest=0
			both=0
			
			for block in strainblocks[strain]:
			
				if strainclosest[strain]!=block[2]:
					both=0
					closest=0
					other=0
					total=0
					for site in snplocations[block[0]:block[1]-1]:
						snpinother='n'
						if sequences[strain][site]!=sequences[block[2]][site]:
							snpinother='y'
							other=other+1
							total=total+1
						if sequences[strain][site]!=sequences[strainclosest[strain]][site]:
							closest=closest+1
							if snpinother=='y':
								both=both+1
							else:
								total=total+1
						
				if closest<=other or strainclosest[strain]==block[2]:
					tempclosest=strainclosest[strain]
				else:
					tempclosest=block[2]
				
				if numblocks>0 and newstrainblocks[strain][numblocks-1][0]==tempclosest:
					newstrainblocks[strain][numblocks-1][2]=block[1]
				else:
					newstrainblocks[strain].append([tempclosest, block[0], block[1], other, closest, both, total])
					numblocks=numblocks+1
				
		
		strainblocks=[]
		
		alltabfile=open("recombinations.tab", 'w')
		print >> alltabfile, 'ID   Recombinations'
		for strain in newstrainblocks.keys():

			tabfile=open(strain+"_closest.tab", 'w')
			print >> tabfile, 'ID   Closest'
			
			for block in newstrainblocks[strain]:
				
				print >> tabfile, "FT   misc_feature    "+str(snplocations[block[1]])+'..'+str(snplocations[block[2]-1])
				print >> tabfile, 'FT                   /note='+block[0]
				if strainclosest[strain]==block[0]:
					print >> tabfile, 'FT                   /colour=0 0 255'
				else:
					#colour is based on 
					#colourcode=(float(block[3])/block[6])/(float(block[4])/block[6])*255
					#colourcode=(float(block[4])-block[3])/(float(block[6])-block[3])*255
					colourcode=(float(block[3]+10))/(float(block[4]+10))*255
					print >> tabfile, 'FT                   /note=', block[3],"snps in", block[0], block[4], "in", strainclosest[strain]
					print >> tabfile, 'FT                   /colour='+str(255-int(colourcode)), 0, str(int(colourcode))
					
					if (float(block[3]+10))/(float(block[4]+10))<0.75 and block[4]>10:
						print >> alltabfile, "FT   misc_feature    "+str(snplocations[block[1]])+'..'+str(snplocations[block[2]-1])
						print >> alltabfile, 'FT                   /note=', strain, "has", block[3],"snps in", block[0], block[4], "in", strainclosest[strain]
						print >> alltabfile, 'FT                   /colour='+str(255-int(colourcode)), 0, str(int(colourcode))
						
			tabfile.close()
		alltabfile.close()		

		sys.exit()
		
		
	
	

	if embl!='' and ref!='':
		print "\nCalculating dN/dS values...\nConcatenating reference CDSs...",
		sys.stdout.flush()
		dndsout=open(outfile+'_dNdS_by_CDS.out','w')
		
		sys.stdout.flush()
		dnbydsstats={}
		AAtypes={}
		
		comp={'A':'T','T':'A','G':'C','C':'G', 'N':'N'}

		
		refseq=sequences[ref].replace('-','')
		refCDSseq=''
		CDSbasenumbers=[]
		
		for key in CDSposns.keys():
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
		for sequence in sequences.keys():
			if sequence==ref:
				continue
			
			print sequence+'...',
			sys.stdout.flush()
			
			tmpseq=''
			for x in range (len(refseq)):
				tmpseq=tmpseq+sequences[sequence][refconvertb[x]]
				
			tmpseqs[sequence]=tmpseq
			
			tmpCDSseq=''
			
			for key in CDSposns.keys():
				for x in CDSposns[key]:
					if x[1]=='f' and tmpseqs[sequence][x[0]-1].upper() in ['A','C','G','T','N']:
						tmpCDSseq=tmpCDSseq+tmpseqs[sequence][x[0]-1].upper()
					elif tmpseqs[sequence][x[0]-1].upper() in ['A','C','G','T','N']:
						tmpCDSseq=tmpCDSseq+comp[tmpseqs[sequence][x[0]-1].upper()]
					else:
						tmpCDSseq=tmpCDSseq+'-'
			
			
			
			dnbydsstats[sequence], snptypes[sequence], AAtypes[sequence]=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers)
			
			#N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd
			if dnbydsstats[sequence]["dS"]!=0:
				print "%.2f" % (dnbydsstats[sequence]["dN"]/dnbydsstats[sequence]["dS"])
				print >> dndsout, '\t'+str(dnbydsstats[sequence]["dN"]/dnbydsstats[sequence]["dS"]),
				sys.stdout.flush()
			else:
				print "-"			
			
		#Testing doing each gene seperately
		
		#for sequence in sequences.keys():
		#		tmpseqs[sequence]=''
		#		for i, j in enumerate(sequences[ref]):
		#			if j!='-':
		#				tmpseqs[sequence]=tmpseqs[sequence]+sequences[sequence][i]
		
		for key in CDSposns.keys():
			refCDSseq=''
			
			for x in CDSposns[key]:
				if x[1]=='f':
					refCDSseq=refCDSseq+refseq[x[0]-1].upper()
				else:
					refCDSseq=refCDSseq+comp[refseq[x[0]-1].upper()]
			
			print >> dndsout, '\n'+CDSnames[key],
			#print CDSnames[key]
			sys.stdout.flush()
			
			for sequence in sequences.keys():
				if sequence==ref:
					continue
					
				tmpCDSseq=''
				
				for y in CDSposns[key]:
					if y[1]=='f' and tmpseqs[sequence][y[0]-1] in ['A','C','G','T','N']:
						tmpCDSseq=tmpCDSseq+tmpseqs[sequence][y[0]-1].upper()
					elif tmpseqs[sequence][y[0]-1] in ['A','C','G','T','N']:
						tmpCDSseq=tmpCDSseq+comp[tmpseqs[sequence][y[0]-1].upper()]
					else:
						tmpCDSseq=tmpCDSseq+'-'
				#print refCDSseq+'\n'+tmpCDSseq+'\n'
				genednbydsstats, gene_snptypes, geneAAtypes=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers)
				if genednbydsstats["dS"]!=0:
					print >> dndsout, '\t'+str(genednbydsstats["dN"]/genednbydsstats["dS"]),
				else:
					print >> dndsout, '\t-',
				
		
		
		dndsout.close()
			
		
		
		
		
		print "\nCalculating pairwise dN/dS values"
		sys.stdout.flush()
		dndsout=open(outfile+'_pairwise_snps.out','w')
		print >> dndsout, 'Sequence1\tSequence1\tPairwise_snps\tNonsynonymous\tSynonymous\tintergenic\tdN/dS'
		sys.stdout.flush()
		dnbydsstats={}
		AAtypes={}
		
		comp={'A':'T','T':'A','G':'C','C':'G', 'N':'N', "-":"-"}

		for i, tmpref in enumerate(sequences.keys()):
			refseq=sequences[tmpref]
			refCDSseq=''
			CDSbasenumbers=[]
			
			for key in CDSposns.keys():
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
			for j in range(i+1,len(sequences.keys())):
				sequence=sequences.keys()[j]
				#if sequence==tmpref:
				#	continue
				
				print tmpref, sequence+'...',
				sys.stdout.flush()
				
				tmpseq=''
#				for x in range (len(refseq)):
#					tmpseq=tmpseq+sequences[sequence][refconvertb[x]]
					
				tmpseqs[sequence]=sequences[sequence]
				
				tmpCDSseq=''
				
				for key in CDSposns.keys():
					for x in CDSposns[key]:
						if x[1]=='f' and tmpseqs[sequence][x[0]-1].upper() in ['A','C','G','T','N']:
							tmpCDSseq=tmpCDSseq+tmpseqs[sequence][x[0]-1].upper()
						elif tmpseqs[sequence][x[0]-1].upper() in ['A','C','G','T','N']:
							tmpCDSseq=tmpCDSseq+comp[tmpseqs[sequence][x[0]-1].upper()]
						else:
							tmpCDSseq=tmpCDSseq+'-'
				
				
				
#				dnbydsstats, snptypes, AAtypes=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers)
#				
#				
#				nsyn=0
#				nnon=0
#				for f in snptypes.keys():
#					if snptypes[f] in ["N", "2", "3"]:
#						nnon=nnon+1
#					elif snptypes[f]=="S":
#						nsyn=nsyn+1
						
				pairwisesnps=0
				pairwiseintergenic=0
				
				
				for x, base in enumerate(sequences[tmpref]):
					if base!= sequences[sequence][x] and base!="-" and sequences[sequence][x]!="-":
						pairwisesnps=pairwisesnps+1
						if embldata[x]=="I":
							pairwiseintergenic=pairwiseintergenic+1
				
				#N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd
#				if dnbydsstats["dS"]!=0 and dnbydsstats["dN"]!=0:
#					print "%.2f" % (dnbydsstats["dN"]/dnbydsstats["dS"])
#					print >> dndsout, tmpref+'\t'+sequence+'\t'+str(pairwisesnps)+'\t'+str(nnon)+'\t'+str(nsyn)+'\t'+str(pairwisesnps-(nnon+nsyn))+'\t'+str(dnbydsstats["dN"]/dnbydsstats["dS"])
#					sys.stdout.flush()
#				else:
#					print >> dndsout, tmpref+'\t'+sequence+'\t'+str(pairwisesnps)+'\t'+str(nnon)+'\t'+str(nsyn)+'\t'+str(pairwisesnps-(nnon+nsyn))+'\t'+"-"
#					print "-"

				print >>dndsout, tmpref+'\t'+sequence+'\t'+str(pairwisesnps)+str(pairwiseintergenic)
			
		
		
		
		
		
		
		
		
			
		print "Done"








	
	#print summary file and alignment file
	
	print "\nWriting output file(s)...",
	sys.stdout.flush()
	
	output=open(outfile+'.out','w')
	
	#if len(allsnps.keys())>1:
	#	print >> output, 'Ref_sequence\tPosition_in_ref',
	#else:
	print >> output, 'Position_in_alignment',
	if ref!='':
		 print >> output, '\tPosition_in_'+ref,
	if embl!='':
	        print >> output, '\tCDS/rRNA/tRNA/Intergenic\tCDS_name\tSynonymous/Non-synonymous',
	print >> output, '\tRef_base\tSNP_base\tTotal',
	
	for name in sequences.keys():
		print >> output, '\t'+name,
	print >> output, '\tPattern\n',
	
	if align=='y':
		aln={}
		alnfile=open(outfile+'_snps.aln', 'w')
	if tabfile=='y':
	        tabout=open(outfile+'_snps.tab','w')
	if graphs=='y':
		snpgraph=open(outfile+'_snps.txt','w')
		#covgraph=open(outfile+'_cov.txt','w')
		#covcountgraph=open(outfile+'_covcount.txt','w')
	
	#minquality=10
	#minqualityperdepth=4
	
	snpsequence={}
	snpsummary={}
	for name in sequences.keys():
		snpsequence[name]=''
		snpsummary[name]={'A':{'C':0, 'G':0, 'T':0}, 'C':{'A':0, 'G':0, 'T':0}, 'G':{'C':0, 'A':0, 'T':0}, 'T':{'C':0, 'G':0, 'A':0}}
	
	
	refsequence=''
	if ref=='':
		ref=sequences.keys()[0]
	
	if tabfile=='y':
		print >> tabout, 'ID   SNP'
	
	poscount=1
	
	
	
	for x, j in enumerate(snplocations):
	
		if graphs=='y':
			while j>poscount:
				print >> snpgraph, '0'
				#print >> covgraph, str(float(avcoverage[poscount])/len(snpfiles))
				#print >> covcountgraph, str(covcount[poscount])
				poscount=poscount+1
		
		
		outstring=''
		tabstring=''
		
		#if len(allsnps.keys())>1:
		#	outstring=outstring+key+'\t'+str(j)
		#else:
		outstring=outstring+str(j+1)
		snpcolour='1'
		if sequences[ref][j]!='-':
			outstring=outstring+'\t'+str(refconvert[j]+1)
			if embl!='':
				
				if embldata[refconvert[j]]!='C':
					outstring=outstring+'\t'+embldata[refconvert[j]]+'\t-'
				else:
					outstring=outstring+'\t'+embldata[refconvert[j]]+'\t'+snptoCDSname[j]
				snptype='-'
				for name in sequences.keys():
					if name==ref:
						continue
					if snptypes[name].has_key(refconvert[j]):
						if snptypes[name][refconvert[j]]=='S':
							snpcolour='3'
						elif snptypes[name][refconvert[j]]=='N':
							snpcolour='2'
						if snptype=='-':
							snptype=snptypes[name][refconvert[j]]
						elif snptypes[name][refconvert[j]] not in snptype and snptypes[name][refconvert[j]]!='1':
							snptype=snptype+'/'+snptypes[name][refconvert[j]]
				outstring=outstring+'\t'+snptype
		else:
			outstring=outstring+'\t-'
			snpcolour='1'
			if embl!='':
				outstring=outstring+'\t-\t-\t-'
		
		
		outstring=outstring+'\t'+sequences[ref][j]
		refbase=sequences[ref][j]
		refsequence=refsequence+refbase
		
			
		snpbase=''
		snpbasecount=0
		for base in snpbases[x].keys():
			if base!=refbase and snpbases[x][base]>snpbasecount:
				snpbasecount=snpbases[x][base]
				snpbase=base
			elif base!=refbase and  snpbases[x][base]==snpbasecount:
				snpbase=snpbase+','+base
	
		outstring=outstring+'\t'+snpbase+'\t'+str(snpbasecount)# need to think how to do this
		
		if tabfile=='y' and sequences[ref][j]!='-':
			tabstring=tabstring+'FT   SNP             '+str(refconvert[j]+1)+'\n'
			tabstring=tabstring+'FT                   /note="refAllele: '+sequences[ref][j]
			#tabstring=tabstring+'FT                   /SNPAllele="'+snpbase+'"\n'
			tabstring=tabstring+' SNPstrains: '

		

		pattern=''
		numsnpbases=1
		sitesnpbases=[]
		for name in sequences.keys():
			
			if sequences[name][j]!='-':# and name.SNPs[key][j][1]>8:# and name.SNPs[key][j][2]>(name.SNPs[key][j][3]*4):
				if sequences[name][j]!=sequences[ref][j]:
					
					if tabfile=='y' and sequences[ref][j]!='-':
						tabstring=tabstring+name+'='+sequences[name][j]+' '
						if embl!='':
							if snptypes[name].has_key(refconvert[j]):
								if snptypes[name][refconvert[j]]=='S':
									tabstring=tabstring+'(synonymous) '
								elif snptypes[name][refconvert[j]]=='N':
									tabstring=tabstring+'(non-synonymous) '
								elif snptypes[name][refconvert[j]]=='1':
									tabstring=tabstring+'(gap in SNP codon) '
								elif snptypes[name][refconvert[j]]=='2':
									tabstring=tabstring+'(SNP codon is STOP) '
								elif snptypes[name][refconvert[j]]=='3':
									tabstring=tabstring+'(ref codon is STOP) '
					pattern=pattern+'1'
					outstring=outstring+'\t'+sequences[name][j]#+str(name.mapped[key][j][1])
					#if sequences[name][j] not in sitesnpbases:
					#	sitesnpbases.append(sequences[name][j])
					#	numsnpbases=numsnpbases+1
					#name.SNPs[j]=''
					
				else:
					pattern=pattern+'0'
					#outstring=outstring+'\t'+sequences[name][j]
					outstring=outstring+'\t.'
				
			else:
				pattern=pattern+'-'
				#snpsequence[name]=snpsequence[name]+'-'
				outstring=outstring+'\t-'
			
			snpsequence[name]=snpsequence[name]+sequences[name][j]
			if sequences[ref][j]!=sequences[name][j] and sequences[ref][j] in ["A", "C", "G", "T"] and sequences[name][j] in ["A", "C", "G", "T"]:
				snpsummary[name][sequences[ref][j]][sequences[name][j]]=snpsummary[name][sequences[ref][j]][sequences[name][j]]+1
		
		
		if graphs=='y':
			print >> snpgraph, '1'
		print >> output, outstring+'\t'+pattern
		if tabfile=='y' and sequences[ref][j]!='-':
			print >> tabout, tabstring+'"\nFT                   /colour='+snpcolour

		if graphs=='y':
			#print >> covgraph, str(float(avcoverage[poscount])/len(snpfiles))
			#print >> covcountgraph, covcount[poscount]
			poscount=poscount+1

	
	if graphs=='y':
		while reflen>=poscount:
			print >> snpgraph, '0'
			#print >> covgraph, str(float(avcoverage[poscount])/len(snpfiles))
			#print >> covcountgraph, str(covcount[poscount])
			poscount=poscount+1
	
	
	if align=='y':
		
		count=0
		
		
		
		for name in snpsequence.keys():
			#print >> alnfile, name.replace('pool','').replace('.snps','')+' '*(10-len(name))+snpsequence[name]
			if len(snpsequence[name].replace("-",""))>(len(snpsequence[name])/3):
				count=count+1
		
		print >> alnfile, count, len(refsequence)
			
		for name in snpsequence.keys():
			#print >> alnfile, name.replace('pool','').replace('.snps','')+' '*(10-len(name))+snpsequence[name]
			if len(snpsequence[name].replace("-",""))>(len(snpsequence[name])/3):
				print >> alnfile, name.replace('pool','').replace('.snps','')+' '+snpsequence[name]	
					
		alnfile.close()
	
	
	
	#Create tab file for each sequence vs ref (extended to length of total alignment to allow insertions)
	#print refconvertb
	if tabfile=='y':
		extrefout=open('Padded_'+ref+'.fna', 'w')
		print >> extrefout, '>'+ref
		count=0
		outstring=''
		for x in range(len(sequences[ref])):
			outstring=outstring+sequences[ref][x]
			count=count+1
			if count==60:
				outstring=outstring+'\n'
				count=0
		print >> extrefout, outstring
		outstring=''
		#print reflen, len(sequences[ref]), len(sequences[ref].replace('-',''))
		extrefout.close()
		if embl!='':
			lines=open(embl,'rU').readlines()
			extrefout=open('Padded_'+ref+'.embl', 'w')
			for line in lines:
				if len(line.strip().split())==3 and line.strip().split()[1]=='source':
					#print '..'+str(len(sequences[ref])), '..'+str(len(sequences[ref].replace('-','')))
					print >> extrefout, line.replace( '..'+str(len(sequences[ref].replace('-',''))), '..'+str(len(sequences[ref]))),
				elif len(line.strip().split())==3 and len(line.strip().split()[2].split('..'))==2:
					positions=line.strip().split()[2]
					#print positions
					num_list = re.findall(r'.[0-9]+.', positions)
					#print num_list
					for numa in num_list:
					
						num = re.findall(r'[0-9]+', numa)[0]
						
						#print numa, num,
						
						numb=numa.replace(num,str(refconvertb[int(num)-1]+1))
						#print numb
						
						positions=positions.replace(numa, numb)
						 
						#print num, num[0]+str(refconvertb[int(num[1:-1])-1]+1)+num[-1], refconvertb[int(num[1:-1])-1]
					
					
					#print positions
					#start=refconvertb[int(line.strip().split()[2].split('..')[0])-1]
					#end=refconvertb[int(line.strip().split()[2].split('..')[1])-1]
					#posn=line[:22]+str(start+1)+'..'+str(end+1)
					print >> extrefout, line[:21]+positions
				elif len(line.strip().split())>1 and line.strip().split()[0]=='SQ':
					A=0
					C=0
					G=0
					T=0
					Others=0
					
					for x in sequences[ref]:
						if x=='A':
							A=A+1
						elif x=='C':
							C=C+1
						elif x=='G':
							G=G+1
						elif x=='T':
							T=T+1
						else:
							Others=Others+1
					print >> extrefout, "SQ   Sequence "+str(len(sequences[ref]))+" BP; "+str(A)+" A; "+str(C)+" C; "+str(G)+" G; "+str(T)+" T; "+str(Others)+" other;"
					
					count=0
					countb=0
					totalsofar=0
					outstring=''
					for x in range(len(sequences[ref])):
						if count==0:
							outstring=outstring+'     '
						outstring=outstring+sequences[ref][x].lower()
						count=count+1
						totalsofar=totalsofar+1
						countb=countb+1
						
						if count==60:
							outstring=outstring+'       '+str(totalsofar)+'\n'
							count=0
							countb=0
						if countb==10:
							outstring=outstring+' '
							countb=0
					
					if count!=0:
						while count!=0:
							outstring=outstring+' '
							count=count+1
							countb=countb+1
							if count==60:
								outstring=outstring+'       '+str(totalsofar)+'\n//'
								count=0
								countb=0
							if countb==10:
								outstring=outstring+' '
								countb=0
	
					
					print >> extrefout, outstring,
					outstring=''
					
					break
					
					
				else:
					print >> extrefout, line.strip()
			extrefout.close()
		for name in sequences.keys():
			indtabfile=open(name+'.tab', 'w')
			print >> indtabfile, 'ID   SNP'
			delstart=0
			insstart=0
			for i,j in enumerate(sequences[name]):
				if j=='-':
					if i==0 or sequences[name][i-1]!='-':
						delstart=i
					if (i+1)==len(sequences[ref]) or sequences[name][i+1]!='-':
						if delstart!=i:
							print >> indtabfile, 'FT   misc_feature    '+str(delstart+1)+'..'+str(i+1)
						else:
							print >> indtabfile, 'FT   misc_feature    '+str(i+1)
						print >> indtabfile, 'FT                   /colour=3'#Deletion = grey
						
#				elif sequences[ref][i]=='-' and j!=sequences[ref][i]:
#					if i==0 or sequences[ref][i-1]!='-' or sequences[sequence][i-1]=='-':
#						insstart=i
					#if (i+1)==len(sequences[ref]) or sequences[ref][i+1]!='-' or sequences[sequence][i+1]=='-':
						#if insstart!=i:
						#	print >> indtabfile, 'FT   misc_feature    '+str(insstart+1)+'..'+str(i+1)
						#else:
						#	print >> indtabfile, 'FT   misc_feature    '+str(i+1)
						#print >> indtabfile, 'FT                   /colour=10'#Insertion = orange
						
				elif j!=sequences[ref][i] and sequences[ref][i]!='-':
					print >> indtabfile, 'FT   misc_feature    '+str(i+1)
					print >> indtabfile, 'FT                   /refBase="'+sequences[ref][i]+'"'
					#tabstring=tabstring+'FT                   /SNPAllele="'+snpbase+'"\n'
					print >> indtabfile, 'FT                   /SNPBase="'+j+'"'
					if embl!='' and snptypes[name].has_key(refconvert[i]):
						if snptypes[name][refconvert[i]]=='S':
							print >> indtabfile, 'FT                   /colour=2'#Synonymous=red
						elif snptypes[name][refconvert[i]]=='N':
							print >> indtabfile, 'FT                   /colour=4'#Non-synonymous=blue
						elif snptypes[name][refconvert[i]]=='1':
							print >> indtabfile, 'FT                   /colour=1'#add synonymous stuff here
						else:
							print >> indtabfile, 'FT                   /colour=1'#add synonymous stuff here
					elif embl!='':
						print >> indtabfile, 'FT                   /colour=1'#add synonymous stuff here
					else:
						print >> indtabfile, 'FT                   /colour=14'#add synonymous stuff here

			indtabfile.close()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	#sys.exit()
	summary='y'
	if summary=='y':
		summaryfile=open(outfile+'_summary.out','w')
		
		print >> summaryfile, 'Strain\tPercent non-gap\tdN/dS',
		
		for x in ['A', 'C', 'G', 'T']:
			for y in ['A', 'C', 'G', 'T']:
				if x!=y:
					print >> summaryfile, '\t'+x+'->'+y,
		
		print >> summaryfile, '\tTotal SNPs'#\t% Mapped Intragenic sites with SNP\t% Mapped Intergenic sites with SNP'
		
		for name in sequences.keys():
			if name==ref:
				continue
			intragenic=0
			intragenic_mapped=0
			totalSNPs=0
			print float(len(sequences[name].replace("-","")))/regionlen, regionlen
			if embl!="" and dnbydsstats[name]["dS"]!=0:
				print >> summaryfile, name+'\t'+str(100*(float(len(sequences[name].replace("-","")))/regionlen))+"\t"+str(dnbydsstats[name]["dN"]/dnbydsstats[name]["dS"]),
			else:
				print >> summaryfile, name+'\t'+str(100*(float(len(sequences[name].replace("-","")))/regionlen))+'\t-',
			for x in ['A', 'C', 'G', 'T']:
				for y in ['A', 'C', 'G', 'T']:
					if x!=y:
						print >> summaryfile, '\t'+str(snpsummary[name][x][y]),
						totalSNPs=totalSNPs+snpsummary[name][x][y]
			
#			for CDS in CDSposns.keys():
#				for y in CDSposns[CDS]:
#					if name.SNPs.has_key(y[0]):
#						intragenic=intragenic+1
#					for chromosome in name.mapped.keys():
#						if name.mapped[chromosome].has_key(y[0]):
#							intragenic_mapped=intragenic_mapped+1
#			intergenic_mapped=name.nummapped-intragenic_mapped
			
			#print intragenic_mapped, intergenic_mapped, name.nummapped
			
#			if intergenic_mapped>0:
#				intergenic=(float(totalSNPs-intragenic)/intergenic_mapped)*100
#			else:
#				intergenic=0
#			if intragenic_mapped>0:
#				intragenic=(float(intragenic)/intragenic_mapped)*100
#			else:
#				intragenic=0
			
			
			print >> summaryfile, '\t'+str(totalSNPs)#+'\t'+str(intragenic)+'\t'+str(intergenic)
		
		print >> summaryfile, '\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(len(snpsequence[ref]))
		summaryfile.close()
		
		print 'Done'
	
	if tabfile=='y':
		tabout.close()
	
	if graphs=='y':
		snpgraph.close()
		#covgraph.close()
		#covcountgraph.close()
	output.close()
	
	#sys.exit()
	if raxml=='y':
		outlen=0
		userinput='x'
		if '/' in outfile:
			outlen=-1*len(outfile.split('/')[-1])
		if os.path.isfile('RAxML_info.'+outfile.split('/')[-1]):
			print '\nRAxML files with extension '+outfile.split('/')[-1]+' already exist!'
			while userinput not in ['y','n']:
				userinput=raw_input('Overwrite? (y/n): ')
				userinput.lower()
			if userinput=='y':
				os.system('rm RAxML_*'+outfile.split('/')[-1])
		print "Running RAxML phylogeny with "+model+" model of evolution and "+str(bootstrap)+" bootstrap replicates..."
		os.system("RAxML -f a -x "+str(random.randrange(1,99999))+" -p "+str(random.randrange(1,99999))+" -# "+str(bootstrap)+" -m "+model+" -s "+outfile+"_snps.aln -n "+outfile.split('/')[-1])
		
	
	
	
	
	print "All finished.\n"
