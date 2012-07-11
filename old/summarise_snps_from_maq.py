#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math


def Usage():
	print 'compare_mummer_out.py Usage:'
	print 'compare_mummer_out.py -r=[reference sequence] [input alignment(s)] > [output file name] {-h}'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "he:r:o:taimcpgb:M:q:s:v", ["model=", "bootstrap=", "help", "phylogeny" "embl=" "ref=", "out=", "align", "maq", "tabfile", "graphs", "indels", "circular", "quality=", "snptype=", "velvet="])
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
	model='GTRGAMMAI'
	quality=8
	snptype='maq'
	velvet='n'

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-m", "--maq"):
			runmaq='y'
		elif opt in ("-t", "--tabfile"):
			tabfile='y'
		elif opt in ("-i", "--indels"):
			indels='y'
		elif opt in ("-c", "--circular"):
			circular='y'
		elif opt in ("-p", "--phylogeny"):
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
		elif opt in ("-s", "--snptype"):
			snptype=arg
		elif opt in ("-v", "--velvet"):
			velvet='y'

	inputdirs=args
	


	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if outfile=='':
		outfile=ref.split('.')[0]
	if align=='n' and raxml=='y':
		print "Can't create phylogeny without alignment!"
		raxml='n'

	return ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet



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

def dnbyds(CDS, SNPseq):
	
	geneticcode={'TTT':'Phe', 'TTC':'Phe', 'TTA':'Leu', 'TTG':'Leu', 'TCT': 'Ser', 'TCC': 'Ser','TCA': 'Ser','TCG': 'Ser', 'TAT': 'Tyr','TAC': 'Tyr', 'TAA': 'STOP', 'TAG': 'STOP', 'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp', 'CTT': 'Leu','CTC': 'Leu','CTA': 'Leu','CTG': 'Leu', 'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'CAT': 'His', 'CAC': 'His', 'CAA': 'Gin', 'CAG': 'Gin', 'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met', 'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val', 'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu', 'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'}
	
	#codonsynonyms={'TTT':1, 'TTC':1, 'TTA':2, 'TTG':2, 'TCT': 3, 'TCC': 3,'TCA': 3,'TCG': 3, 'TAT': 1,'TAC': 1, 'TAA': 2, 'TAG': 1, 'TGT': 1, 'TGC': 1, 'TGA': 1, 'TGG': 0, 'CTT': 3,'CTC': 3,'CTA': 4,'CTG': 4, 'CCT': 3, 'CCC': 3, 'CCA': 3, 'CCG': 3, 'CAT': 1, 'CAC': 1, 'CAA': 1, 'CAG': 1, 'CGT': 3, 'CGC': 3, 'CGA': 5, 'CGG': 5, 'ATT': 2, 'ATC': 2, 'ATA': 2, 'ATG': 0, 'ACT': 3, 'ACC': 3, 'ACA': 3, 'ACG': 3, 'AAT': 1, 'AAC': 1, 'AAA': 1, 'AAG': 1, 'AGT': 1, 'AGC': 1, 'AGA': 1, 'AGG': 1, 'GTT': 3, 'GTC': 3, 'GTA': 3, 'GTG': 3, 'GCT': 3, 'GCC': 3, 'GCA': 3, 'GCG': 3, 'GAT': 1, 'GAC': 1, 'GAA': 1, 'GAG': 1, 'GGT': 3, 'GGC': 3, 'GGA': 3, 'GGG': 3}
	
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
					elif geneticcode[newcodon]=='STOP':
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
	varS=0.0
	varN=0.0
	z=0.0
	dN=0.0
	dS=0.0
	
	if len(CDS)!=len(SNPseq):
		print "Error: sequences must be the same length to calculate dN/dS!"
		return N, S, dN, dS, pN, pS, varS, varN, z, (len(CDS)-gapcount), Nd, Sd
	
	
	
	
	for x in range(0,len(CDS),3):
		numcodons=numcodons+1
		codon=CDS[x:x+3]
		SNPcodon=SNPseq[x:x+3]
		
		
		if '-' in codon or '-' in SNPcodon:
			gapcount=gapcount+3
			continue
		
		if geneticcode[codon]=='STOP' or geneticcode[SNPcodon]=='STOP':
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
		
		#if nochanges==2:
		#	sd=float(sd)/2
		#	nd=float(nd)/2
		#elif nochanges==3:
		#	sd=float(sd)/6
		#	nd=float(nd)/6
		if pathcount>0:
			sd=float(sd)/pathcount
			nd=float(nd)/pathcount
		#elif codon!= SNPcodon:
		#	print codon, SNPcodon, "All paths lead to early stop codon in SNP sequence."
			#return
		
		
		Sd=Sd+sd
		Nd=Nd+nd
		#pS=pS+(sd/s)
		#pN=pN+(nd/n)
		
		#print codon, SNPcodon
	
	S=(S1+S2)/2
	N=len(CDS)-gapcount-S
	
	#pNb=pN/numcodons
	#pSb=pS/numcodons
	
	if N!=0:
		pN=Nd/N
	if S!=0:
		pS=Sd/S
	
	if pS==0:
		print "No sites are synonymous."
		return N, S, dN, dS, pN, pS, varS, varN, z, (len(CDS)-gapcount), Nd, Sd 
	if pS<0.75 and pN<0.75:
		dS=(-3*(math.log(1-((pS*4)/3))))/4
		dN=(-3*(math.log(1-((pN*4)/3))))/4
		varS=(9 * pS * (1 -pS))/(((3 - 4 *pS) **2) * (len(CDS)-gapcount));
		varN=(9 * pN * (1 -pN))/(((3 - 4 *pN) **2) * (len(CDS)-gapcount));
		z=(dN - dS) / math.sqrt(varS + varN)
		
	else:
		print "Too divergent for JC! Using pN/pS instead."
		dS=pS
		dN=pN
		
	
	
	#print dN, dS, S, N, S+N, Sd, Nd, pS, pN#, pSb, pNb, dS, dN, pN/pS, dN/dS
	#print "N =", N
	#print "S =", S
	#print "dN/dS =", dN/dS
	
	return N, S, dN, dS, pN, pS, varS, varN, z, (len(CDS)-gapcount), Nd, Sd
	



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
		self.varN=0.0
		self.varS=0.0
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
        ref, inputdirs, outfile, tabfile, align, embl, raxml, runmaq, graphs, bootstrap, model, quality, snptype, velvet=getOptions(argv)

		

snps={}

refbases={}
bases={}#0=A 1=C 2=G 3=T
nostates=0
converter={}
convertback={}
snpbases={}

snpfiles=[]
count=0
namesort=[]
runpileups='n'

print '\nChecking input files...'
sys.stdout.flush()

for inputfile in inputdirs:
	if inputfile[-1]=='/':
		inputfile=inputfile[:-1]
	print inputfile+'...',
	sys.stdout.flush()
	if os.path.isdir(inputfile):
		userinput='x'
		if runmaq=='n':
			if not os.path.isfile(inputfile+'/all.map'):
				print "\nError: No Maq map ("+inputfile+"/all.map) found!"
				while userinput not in ['i','c','a']:
					userinput=raw_input("Would you like to ignore current folder, run maq on current folder or run maq on all folders? (i/c/a): ")
					userinput.lower()
				if userinput=='a':
					runmaq='y'

			elif not os.path.isfile(inputfile+'/all.pileup'):
				if runpileups=='n':
					print "\nError: No pileup file ("+inputfile+"/all.pileup) found!"
					while userinput not in ['y','n','a']:
						userinput=raw_input("Would you like to create it on the fly? (y/n/a): ")
						userinput.lower()
					if userinput=='a':
						runpileups='y'
						userinput='y'
				else:
					userinput=='y'
			else:
				templine=open(inputfile+"/all.pileup", "rU").readlines(1)[0].strip().split()
				if len(templine)!=7:
					if runpileups=='n':
						print "\nError: Pileup file ("+inputfile+"/all.pileup) of wrong format!"
						while userinput not in ['y','n','a']:
							userinput=raw_input("Would you like to recreate it on the fly? (y/n/a): ")
							userinput.lower()
						if userinput=='a':
							runpileups='y'
							userinput='y'
					else:
						userinput='y'
				
		if runmaq=='y' or userinput in['c','y','x']:
			snpfiles.append(SNPanalysis())
			snpfiles[count].directory=inputfile
			namesort.append(inputfile)
			if userinput=='c':
				snpfiles[count].runmaq='m'
			elif userinput=='y':
				snpfiles[count].runmaq='p'
			print 'ok'
			sys.stdout.flush()
			count=count+1
			
		else:
			print 'excluded'
			sys.stdout.flush()
	else:
		print "\nError: Directory "+inputfile+" not found!"

if len(snpfiles)==0:
	print "\nError: No valid input folders!"
	sys.exit()

if snptype=='saha':
	mindepth=float(quality)/5
else:
	mindepth=float(quality)

if velvet=='y':
	output=open(outfile+'_unmapped_contig_hits.csv','w')

for inputfile in snpfiles:
	if runmaq=='y' or inputfile.runmaq=='m':
		print "\nRunning Maq on "+inputfile.directory+'...'
		#print "/nfs/team81/sh16/scripts/maq.pl easyrun -u -d "+inputfile.directory+" "+ref+" fastqs/"+inputfile.directory+".fastq"
		#os.system("/nfs/team81/sh16/scripts/maq.pl easyrun -u "+inputfile.directory+"/unmap.fasta -d "+inputfile.directory+" "+ref+" fastqs/"+inputfile.directory+".fastq")
		print "maq.pl easyrun -e "+str(int(mindepth))+" -d "+inputfile.directory+"test "+ref+" fastqs/"+inputfile.directory+".fastq"
		os.system("maq.pl easyrun -e "+str(int(mindepth))+" -d "+inputfile.directory+"test "+ref+" fastqs/"+inputfile.directory+".fastq")
		
	if velvet=='y':
		print "\nRunning velvet on unassembled reads from "+inputfile.directory+'...'
		
		os.system("cat "+inputfile.directory+"test/unmap*@*.txt > "+inputfile.directory+"test/unmap.txt")
		unmapout=open(inputfile.directory+"test/unmap.fastq", "w")
		lines=open(inputfile.directory+"test/unmap.txt", "rU").readlines()
		for line in lines:
			words=line.split()
			print >> unmapout, '@'+words[0]
			print >> unmapout, words[2]
			print >> unmapout, '+'
			print >> unmapout, words[3]
		unmapout.close()
		
		os.system("/nfs/team81/tdo/bin/velvet/0.7.26/velveth "+inputfile.directory+"test/velvet 21 -fastq -short "+inputfile.directory+"test/unmap.fastq")
		os.system("/nfs/team81/tdo/bin/velvet/0.7.26/velvetg "+inputfile.directory+"test/velvet -cov_cutoff "+str(mindepth)+" -min_contig_lgth 100")
		print "\nRunning blast on contigs from unassembled reads from "+inputfile.directory+'...'
		
		os.system("blastall -p blastn -o "+inputfile.directory+"test/velvet/contigs.blast -d /data/blastdb/Bacteria_DB -v 1 -b 1 -m 8 -i "+inputfile.directory+"test/velvet/contigs.fa")
		
		lines=open(inputfile.directory+"test/velvet/contigs.blast", "rU").readlines()
		lastcontig=''
		
		for line in lines:
			words=line.strip().split()
			if words[0]==lastcontig:
				continue
			else:
				lastcontig=words[0]
			
			fastacmdout=os.popen("fastacmd -d /data/blastdb/Bacteria_DB -s "+words[1].split('|')[1])
			name=fastacmdout.readlines()[0]
			print >> output, inputfile.directory+','+words[0]+','+words[1]+','+'_'.join(name[1:].split()[1:]).replace(',','')+','+','.join(words[2:])

if velvet=='y':
	output.close()
sys.exit()
		

lenstring=os.popen('grep "reference length:" '+snpfiles[0].directory+'/assemble.log').read()
reflen=int(lenstring.strip().split(":")[1])





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
				while foundname=='n' and lines[y].split()[1] not in ('CDS', 'rRNA', 'tRNA', 'gene', 'misc_feature', 'sig_peptide', 'repeat_unit', 'repeat_region', 'misc_RNA'):
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
				embldata, CDSposns=getCDSseq(embldata, words[2].replace('<','').replace('>',''), direction, CDScount, feature)
				
			donelen=donelen+chromosomelen
	
	snptoCDSname={}
	for i in CDSposns.keys():
		for j in CDSposns[i]:
			snptoCDSname[refconvertb[j[0]]]=CDSnames[i]
		
	print "Done"
	sys.stdout.flush()


namesort.sort()

avcoverage={}
covcount={}

#hassnps={}
allsnps={}
refbases={}

#print "\nIdentifying SNPs in pileup files..."
#for inputfile in snpfiles:
#	print inputfile.directory+'...',
#	sys.stdout.flush()
#	if runmaq=='y' or inputfile.runmaq in ['m','p']:
#		print 'Running Maq pileup...',
#		sys.stdout.flush()
#		os.system('maq pileup -v '+inputfile.directory+'/ref.bfa '+inputfile.directory+'/all.map > '+inputfile.directory+'/all.pileup')
	
#	lines=open(inputfile.directory+'/all.pileup', 'rU').readlines()
#	
#	for line in lines:
#		words=line.split()
#		snplocation=int(words[1])
#		readdepth=int(words[3])
#		reads=words[4][1:]
#		
#		if avcoverage.has_key(snplocation):
#			avcoverage[snplocation]=avcoverage[snplocation]+readdepth
#			if readdepth>0:
#				covcount[snplocation]=covcount[snplocation]+1
#		else:
#			avcoverage[snplocation]=readdepth
#			if readdepth>0:
#				covcount[snplocation]=1
#			else:
#				covcount[snplocation]=0
#		if readdepth>1:
#			if len(reads.replace(',','').replace('.',''))>1:
#				hassnps[snplocation]=1
#
#	print 'Done'
#	sys.stdout.flush()
	


alssnpsummary={}



if snptype=='saha':
	print "\nCalculating Ssaha-like mapping qualities and identifying SNPs..."
	sys.stdout.flush()
	for inputfile in snpfiles:
					
		lines=open(inputfile.directory+'/all.pileup', 'rU').readlines()
		tempsnps={}
		tempsnpsummary={'A':{'C':0, 'G':0, 'T':0}, 'C':{'A':0, 'G':0, 'T':0}, 'G':{'C':0, 'A':0, 'T':0}, 'T':{'C':0, 'G':0, 'A':0}}
		mapped=0.0
		olddone=-1
		#for x in hassnps.keys():
			#line=lines[x-1]
		for line in lines:
			
			words=line.split()
			chromosome=words[0]
			snplocation=int(words[1])
			refbase=words[2]
			readdepth=int(words[3])
			reads=words[4][1:]
			
			if not allsnps.has_key(chromosome):
				allsnps[chromosome]={}
			if not tempsnps.has_key(chromosome):
				tempsnps[chromosome]={}
			
			if avcoverage.has_key(snplocation):
				avcoverage[snplocation]=avcoverage[snplocation]+readdepth
			else:
				avcoverage[snplocation]=readdepth
				covcount[snplocation]=0
	
	
			if readdepth>mindepth:
				phred=words[5][1:]
				mapscore=words[6][1:]
				basescores={}
				
				for i, base in enumerate(reads):
					if base in [',','.']:
						base=refbase
					else:
						base=base.upper()
		
					if ord(phred[i])>=30:
						F=1
					else:
						F=0.5
						
					currmapscore=ord(mapscore[i])
					if currmapscore>50:
						currmapscore=50
					iscore=(currmapscore*F)/10
			
					if basescores.has_key(base):
						basescores[base]=basescores[base]+iscore
					else:
						basescores[base]=iscore
			
		
				bestscore=0
				otherscores=0
				bestSNP=''
				for base in basescores.keys():
					if basescores[base]>bestscore:
						otherscores=otherscores+bestscore
						bestscore=basescores[base]
						bestSNP=base
					else:
						otherscores=otherscores+basescores[base]
				
				SNPscore=bestscore-otherscores
				
				if SNPscore>quality:
					if bestSNP!=refbase:
						percentdone=int((float(snplocation)/reflen)*100)
						if percentdone!=olddone:
							print inputfile.directory+'... %d%%\r' % percentdone,
							sys.stdout.flush()
							olddone=percentdone
						
						if allsnps[chromosome].has_key(snplocation):
							allsnps[chromosome][snplocation][0]=allsnps[chromosome][snplocation][0]+1
							allsnps[chromosome][snplocation][2][bestSNP]=allsnps[chromosome][snplocation][2][bestSNP]+1
						else:
							allsnps[chromosome][snplocation]=[1,refbase,{'A':0, 'C':0, 'G':0, 'T':0}]
							allsnps[chromosome][snplocation][2][bestSNP]=1
						tempsnpsummary[refbase][bestSNP]=tempsnpsummary[refbase][bestSNP]+1
						
					tempsnps[chromosome][snplocation]= bestSNP#[bestSNP, SNPscore]
					covcount[snplocation]=covcount[snplocation]+1
					mapped=mapped+1
				
				
					
		#inputfile.mapped=tempsnps
		inputfile.snpsummary=tempsnpsummary
		inputfile.nummapped=mapped
		inputfile.percentmapped=((mapped/reflen)*100)
		print inputfile.directory+'... %.2f%% mapped' % inputfile.percentmapped
		sys.stdout.flush()

else:
	print "\nIdentifying SNPs..."
	sys.stdout.flush()
	for inputfile in snpfiles:
		print inputfile.directory+'...',
		sys.stdout.flush()
		lines=open(inputfile.directory+'/cns.final.snp', 'rU').readlines()
		
		tempsnps={}
		tempsnpsummary={'A':{'C':0, 'G':0, 'T':0}, 'C':{'A':0, 'G':0, 'T':0}, 'G':{'C':0, 'A':0, 'T':0}, 'T':{'C':0, 'G':0, 'A':0}}
		for line in lines:
			words=line.split()
			chromosome=words[0]
			snplocation=int(words[1])
			refbase=words[2]
			cnsbase=words[3]
			
			if cnsbase not in ['A','C','G','T']:
				cnsbase=words[9]
				if cnsbase not in ['A','C','G','T'] or int(words[8])<30 or cnsbase==refbase:
					continue
			
			if not allsnps.has_key(chromosome):
				allsnps[chromosome]={}
			if not tempsnps.has_key(chromosome):
				tempsnps[chromosome]={}
			
			if allsnps[chromosome].has_key(snplocation):
				allsnps[chromosome][snplocation][0]=allsnps[chromosome][snplocation][0]+1
				allsnps[chromosome][snplocation][2][cnsbase]=allsnps[chromosome][snplocation][2][cnsbase]+1
			else:
				allsnps[chromosome][snplocation]=[1,refbase,{'A':0, 'C':0, 'G':0, 'T':0}]
				allsnps[chromosome][snplocation][2][cnsbase]=1
			tempsnpsummary[refbase][cnsbase]=tempsnpsummary[refbase][cnsbase]+1
			tempsnps[chromosome][snplocation]=cnsbase
		
		
		lines=open(inputfile.directory+'/all.pileup', 'rU').readlines()
		mapped=0.0
		
		#for x in hassnps.keys():
			#line=lines[x-1]
		for line in lines:
			
			words=line.split()
			chromosome=words[0]
			snplocation=int(words[1])
			refbase=words[2]
			readdepth=int(words[3])
			reads=words[4][1:]
			
			if not allsnps.has_key(chromosome):
				allsnps[chromosome]={}
			if not tempsnps.has_key(chromosome):
				tempsnps[chromosome]={}
			
			if avcoverage.has_key(snplocation):
				avcoverage[snplocation]=avcoverage[snplocation]+readdepth
			else:
				avcoverage[snplocation]=readdepth
				covcount[snplocation]=0
	
	
			if readdepth>mindepth:
				covcount[snplocation]=covcount[snplocation]+1
				mapped=mapped+1
				if not tempsnps[chromosome].has_key(snplocation):
					tempsnps[chromosome][snplocation]=refbase
				
					
		inputfile.mapped=tempsnps
		inputfile.snpsummary=tempsnpsummary
		inputfile.nummapped=mapped
		inputfile.percentmapped=((mapped/reflen)*100)
		print '%.2f%% mapped' % inputfile.percentmapped
		sys.stdout.flush()
		



if embl!='':
	print "\nCalculating dN/dS values...\n"+ref+'...',
	sys.stdout.flush()
	
	comp={'A':'T','T':'A','G':'C','C':'G'}
	
	reflines=open(ref, "rU").readlines()
	refseq=''
	for line in reflines:
		if line.strip()[0]!='>':
			refseq=refseq+line.strip()
	refCDSseq=''
	for key in CDSposns.keys():
		for x in CDSposns[key]:
			if x[1]=='f':
				refCDSseq=refCDSseq+refseq[x[0]-1].upper()
			else:
				refCDSseq=refCDSseq+comp[refseq[x[0]-1].upper()]
	
	print "Done"
	sys.stdout.flush()
	for inputfile in snpfiles:
		print inputfile.directory+'...',
		sys.stdout.flush()
		lines=open(inputfile.directory+'/all.pileup', 'rU').readlines()
		tempseq=''
		for key in CDSposns.keys():
			for x in CDSposns[key]:
				CDSmapped='n'
				#print x[0],
				
				for chromosome in inputfile.mapped.keys():#NOTE this won't work with more than one chromosome!
					if inputfile.mapped[chromosome].has_key(x[0]):
						CDSmapped=inputfile.mapped[chromosome][x[0]]#[0]
						#print inputfile.mapped[chromosome][x[0]][0]
					
				if CDSmapped!='n':
					#print CDSmapped
					if x[1]=='f':
						tempseq=tempseq+CDSmapped.upper()
					else:
						tempseq=tempseq+comp[CDSmapped.upper()]
					
				else:
					tempseq=tempseq+'-'
#			
			
#			for x in CDSposns[key]:
#				line=lines[x[0]-1]
#				words=line.split()
#				chromosome=words[0]
#				snplocation=int(words[1])
#				refbase=words[2]
#				readdepth=int(words[3])
#				reads=words[4][1:]
#				phred=words[5][1:]
#				mapscore=words[6][1:]
#				basescores={}			
#				for i, base in enumerate(reads):
#					if base in [',','.']:
#						base=refbase
#					else:
#						base=base.upper()
		
#					if ord(phred[i])>=30:
#						F=1
#					else:
#						F=0.5
						
#					currmapscore=ord(mapscore[i])
#					if currmapscore>50:
#						currmapscore=50
#					iscore=(currmapscore*F)/10
			
#					if basescores.has_key(base):
#						basescores[base]=basescores[base]+iscore
#					else:
#						basescores[base]=iscore
			
		
#				bestscore=0
#				otherscores=0
#				bestSNP=''
#				for base in basescores.keys():
#					if basescores[base]>bestscore:
#						otherscores=otherscores+bestscore
#						bestscore=basescores[base]
#						bestSNP=base
#					else:
#						otherscores=otherscores+basescores[base]
#				
#				SNPscore=bestscore-otherscores
#				
#				if SNPscore>8:
#					if x[1]=='f':
#						tempseq=tempseq+bestSNP.upper()
#					else:
#						tempseq=tempseq+comp[bestSNP.upper()]
					
#				else:
#					tempseq=tempseq+'-'
			
			#print key, tempseq
			
		#inputfile.CDSseq=tempsnps
		
		#print refCDSseq, tempseq
		#print len(refCDSseq), len(tempseq)
		
		inputfile.N, inputfile.S, inputfile.dN, inputfile.dS, inputfile.pN, inputfile.pS, inputfile.varS, inputfile.varN, inputfile.z, inputfile.goodlen, inputfile.Nd, inputfile.Ns=dnbyds(refCDSseq, tempseq)
		
		inputfile.CDSseq=tempseq
		
		if inputfile.dS!=0:
			print "%.2f" % (inputfile.dN/inputfile.dS)
		sys.stdout.flush()


#print summary file and alignment file

print "\nWriting output file(s)...",
sys.stdout.flush()

output=open(outfile+'.out','w')

if len(allsnps.keys())>1:
	print >> output, 'Ref_sequence\tPosition_in_ref',
else:
	print >> output, 'Position_in_ref',
if embl!='':
        print >> output, '\tCDS/rRNA/tRNA/Intergenic',
print >> output, '\tRef_base\tSNP_base\tTotal',

for name in namesort:
	print >> output, '\t'+name,
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

minquality=10
minqualityperdepth=4

sequence={}
for name in snpfiles:
	sequence[name.directory]=''

refsequence=''



if tabfile=='y':
	print >> tabout, 'ID   SNP'

poscount=1

for key in allsnps.keys():
	
	

	snpsort=allsnps[key].keys()
	snpsort.sort()

	for j in snpsort:
	
		if graphs=='y':
			while j>poscount:
				print >> snpgraph, '0'
				print >> covgraph, str(float(avcoverage[poscount])/len(snpfiles))
				print >> covcountgraph, str(covcount[poscount])
				poscount=poscount+1
		
		
		outstring=''
		tabstring=''
		
		if len(allsnps.keys())>1:
			outstring=outstring+key+'\t'+str(j)
		else:
			outstring=outstring+str(j)
		if embl!='':
                        outstring=outstring+'\t'+embldata[j]

		outstring=outstring+'\t'+allsnps[key][j][1]

		refsequence=refsequence+allsnps[key][j][1]
		snpbase=''
		snpbasecount=0
		for base in allsnps[key][j][2].keys():
			if allsnps[key][j][2][base]>snpbasecount:
				snpbasecount=allsnps[key][j][2][base]
				snpbase=base
			elif allsnps[key][j][2][base]==snpbasecount:
				snpbase=snpbase+','+base
	
		outstring=outstring+'\t'+snpbase+'\t'+str(allsnps[key][j][0])
		
		if tabfile=='y':
                	tabstring=tabstring+'FT   SNP             '+str(j)+'\n'
			tabstring=tabstring+'FT                   /refAllele="'+allsnps[key][j][1]+'"\n'
			#tabstring=tabstring+'FT                   /SNPAllele="'+snpbase+'"\n'
			tabstring=tabstring+'FT                   /SNPstrains="'

		

		pattern=''
		numsnpbases=1
		snpbases=[]
		for name in snpfiles:
			
			if name.mapped[key].has_key(j):# and name.SNPs[key][j][1]>8:# and name.SNPs[key][j][2]>(name.SNPs[key][j][3]*4):
				if name.mapped[key][j]!=allsnps[key][j][1]:
					
					if tabfile=='y':
						tabstring=tabstring+name.directory+'='+name.mapped[key][j]+' '
					pattern=pattern+'1'
					outstring=outstring+'\t1'#+str(name.mapped[key][j][1])
					if name.mapped[key][j] not in snpbases:
						snpbases.append(name.mapped[key][j])
						numsnpbases=numsnpbases+1
					name.SNPs[j]=''
					
				else:
					pattern=pattern+'0'
					outstring=outstring+'\t0'
				sequence[name.directory]=sequence[name.directory]+name.mapped[key][j]
			else:
				pattern=pattern+'-'
				sequence[name.directory]=sequence[name.directory]+'?'
				outstring=outstring+'\t-'
		
		if numsnpbases>1:
			if graphs=='y':
				print >> snpgraph, '1'
			print >> output, outstring+'\t'+pattern
			if tabfile=='y':
                        	print >> tabout, tabstring+'"\nFT                   /colour='+str(numsnpbases)

		else:
			if graphs=='y':
				print >> snpgraph, '0'
			refsequence=refsequence[:-1]
			for name in snpfiles:
				sequence[name.directory]=sequence[name.directory][:-1]
		if graphs=='y':
			print >> covgraph, str(float(avcoverage[poscount])/len(snpfiles))
			print >> covcountgraph, covcount[poscount]
			poscount=poscount+1


if graphs=='y':
	while reflen>=poscount:
		print >> snpgraph, '0'
		print >> covgraph, str(float(avcoverage[poscount])/len(snpfiles))
		print >> covcountgraph, str(covcount[poscount])
		poscount=poscount+1


if align=='y':
	
	count=0
	for name in snpfiles:
		if name.percentmapped>50:
			count=count+1
	
	print >> alnfile, count+1, len(refsequence)
	print >> alnfile, 'Reference\t'+refsequence
		
	for name in snpfiles:
		if name.percentmapped>50:
			print >> alnfile, name.directory.replace('pool','').replace('.snps','')+'\t'+sequence[name.directory]

				
	alnfile.close()




summaryfile=open(outfile+'_summary.out','w')

print >> summaryfile, 'Strain\tPercent of Reference Mapped\tdN/dS',

for x in ['A', 'C', 'G', 'T']:
	for y in ['A', 'C', 'G', 'T']:
		if x!=y:
			print >> summaryfile, '\t'+x+'->'+y,

print >> summaryfile, '\tTotal SNPs'

if embl!='':
	print >> summaryfile, '\t% Mapped Intragenic sites with SNP\t% Mapped Intergenic sites with SNP'

for name in snpfiles:
	intragenic=0
	intragenic_mapped=0
	if name.dS!=0:
		print >> summaryfile, name.directory+'\t'+str(name.percentmapped)+'\t'+str(name.dN/name.dS),
	else:
		print >> summaryfile, name.directory+'\t'+str(name.percentmapped)+'\t-',
	for x in ['A', 'C', 'G', 'T']:
		for y in ['A', 'C', 'G', 'T']:
			if x!=y:
				print >> summaryfile, '\t'+str(name.snpsummary[x][y]),
				count=count+name.snpsummary[x][y]
	
	totalSNPs=len(name.SNPs.keys())
	
	if embl!='':
		for CDS in CDSposns.keys():
			for y in CDSposns[CDS]:
				if name.SNPs.has_key(y[0]):
					intragenic=intragenic+1
				for chromosome in name.mapped.keys():
					if name.mapped[chromosome].has_key(y[0]):
						intragenic_mapped=intragenic_mapped+1
		intergenic_mapped=name.nummapped-intragenic_mapped
	
	#print intragenic_mapped, intergenic_mapped, name.nummapped
	
	
		if intergenic_mapped>0:
			intergenic=(float(totalSNPs-intragenic)/intergenic_mapped)*100
		else:
			intergenic=0
		if intragenic_mapped>0:
			intragenic=(float(intragenic)/intragenic_mapped)*100
		else:
			intragenic=0
	
	print >> summaryfile, '\t'+str(totalSNPs),
	
	if embl!='':
		print >> summaryfile, '\t'+str(intragenic)+'\t'+str(intergenic)

if embl!='':
	print >> summaryfile, '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(len(snpsort))
else:
	print >> summaryfile, '\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(len(snpsort))


print 'Done'


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
	os.system("RAxML -f a -x "+str(random.randrange(1,99999))+" -p "+str(random.randrange(1,99999))+" -# "+str(bootstrap)+" -m "+model+" -s "+outfile+".aln -n "+outfile.split('/')[-1])
	
if tabfile=='y':
	tabout.close()

if graphs=='y':
	snpgraph.close()
	covgraph.close()
	covcountgraph.close()
output.close()
summaryfile.close()


print "All finished.\n"
