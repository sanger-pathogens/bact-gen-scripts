#!/usr/bin/env python

#Compares mummer snp output files.  NOTE: must be from the same reference sequence



#------------------------------------------------------------------------------------
# Import modules
#------------------------------------------------------------------------------------

import string, re
import os, sys, getopt, random, math


def Usage():
	print 'compare_mummer_out.py Usage:'
	print 'compare_mummer_out.py -r=[reference sequence] [input alignment(s)] > [output file pool] {-h}'
	print 'Copyright Simon R Harris, Wellcome Trust Sanger Institute, UK. 2008'

	

#------------------------------------------------------------------------------------
# Get command line arguments
#------------------------------------------------------------------------------------


def getOptions(arg):

	try:
		opts, args = getopt.getopt(argv, "he:r:o:taimcpgb:M:q:s:vR:d:", ["model=", "bootstrap=", "help", "phylogeny" "embl=" "ref=", "out=", "align", "maq", "tabfile", "graphs", "indels", "circular", "quality=", "snptype=", "velvet=", "rtype=", "depth="])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)
	
	ref=''
	outfile=''
	inputdirs=[]
	i=0
	runssaha='n'
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
	quality=60
	snptype='maq'
	velvet='n'
	rtype='454'
	depth=10

	for opt, arg in opts:
	
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-r", "--ref"):
			ref=arg
		elif opt in ("-o", "--out"):
			outfile=arg
		elif opt in ("-m", "--maq"):
			runssaha='y'
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
			snptype=arg.lower()
		elif opt in ("-v", "--velvet"):
			velvet='y'
		elif opt in ("-R", "--rtype"):
			rtype=arg.lower()
		elif opt in ("-d", "--depth"):
			depth=int(arg)

	inputdirs=args
	


	if inputdirs==[]:
		print 'Error: No input files selected!'
		Usage()
		sys.exit()
	if rtype not in ["abi", "solexa", "454"]:
		print "rtype must be abi, solexa or 454"
		Usage()
		sys.exit()
	if outfile=='':
		outfile=ref.split('.')[0]
	if align=='n' and raxml=='y':
		print "Can't create phylogeny without alignment!"
		raxml='n'

	return ref, inputdirs, outfile, tabfile, align, embl, raxml, runssaha, graphs, bootstrap, model, quality, snptype, velvet, rtype, depth



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
	varianceS=0.0
	varianceN=0.0
	z=0.0
	dN=0.0
	dS=0.0
	SNPtype={}
	
	if len(CDS)!=len(SNPseq):
		print "Error: sequences must be the same length to calculate dN/dS!"
		return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}, SNPtype
	
	
	
	
	for x in range(0,len(CDS),3):
		numcodons=numcodons+1
		codon=CDS[x:x+3]
		SNPcodon=SNPseq[x:x+3]
		codonposn=[CDSbasenumbers[x], CDSbasenumbers[x+1], CDSbasenumbers[x+2]]
		
		if codon!=SNPcodon:
			
			for y,z in enumerate(codon):
				if SNPcodon[y]!=z:
					newSNPcodon=codon[:y]+SNPcodon[y]+codon[y+1:]
					if '-' in SNPcodon:
						SNPtype[codonposn[y]]='1'
					elif geneticcode[newSNPcodon]=='STOP':
						SNPtype[codonposn[y]]='2'
					elif geneticcode[newSNPcodon]=='STOP':
						SNPtype[codonposn[y]]='3'
					elif geneticcode[newSNPcodon] == geneticcode[codon]:
						SNPtype[codonposn[y]]='S'
					else:
						SNPtype[codonposn[y]]='N'
		
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
		#print "No sites are synonymous."
		return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}, SNPtype
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
		
	
	
	#print dN, dS, S, N, S+N, Sd, Nd, pS, pN#, pSb, pNb, dS, dN, pN/pS, dN/dS
	#print "N =", N
	#print "S =", S
	#print "dN/dS =", dN/dS
	
	return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}, SNPtype
	



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
	def __init__(self, fastq='', name='', mapped={}, runssaha='n', CDSseq=''):
		self.fastq=fastq
		self.name=name
		self.mapped=mapped
		self.runssaha=runssaha
		self.CDSseq=CDSseq
		self.dNdSstats={}
		self.SNPtypes={}
		self.goodlen=0
		self.CDSseq=''
		self.snpsummary={}
		self.nummapped=0
		self.percentmapped=0.0
		self.intragenic=0
		self.intergenic=0
		self.SNPs=0
		self.sequence=""



if __name__ == "__main__":
        argv=sys.argv[1:]
        ref, inputdirs, outfile, tabfile, align, embl, raxml, runssaha, graphs, bootstrap, model, quality, snptype, velvet, rtype, depth=getOptions(argv)

		

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

	if pool[-1]=='/':
		pool=pool[:-1]
	pool=pool.split('/')[-1].split('.')[0]
	print pool+'...',
	sys.stdout.flush()
	if not os.path.isdir(pool):
	 	print "Directory "+pool+" not found! Creating...",
		os.system("mkdir "+pool)
		runssaha='y'

	userinput='x'
	if runssaha=='n':
		if not os.path.isfile(pool+'/all.snp'):
			print "\nError: No Ssaha snp file ("+pool+"/all.snp) found!"
			while userinput not in ['i','c','a']:
				userinput=raw_input("Would you like to ignore current folder, run Ssaha on current folder or run Ssaha on all folders? (i/c/a): ")
				userinput.lower()
			if userinput=='a':
				runssaha='y'

		elif not os.path.isfile(pool+'/all.pileup'):
			if runpileups=='n' and runssaha=='n':
				print "\nError: No Ssaha snp file ("+pool+"/all.snp) found!"
				while userinput not in ['i','c','a']:
					userinput=raw_input("Would you like to ignore current folder, run Ssaha on current folder or run Ssaha on all folders? (i/c/a): ")
					userinput.lower()
				if userinput=='a':
					runssaha='y'

			
	if runssaha=='y' or userinput in['c','y','x']:
		pools.append(SNPanalysis())
		pools[count].name=pool
		poolsort.append(pool)
		if userinput=='c':
			pools[count].runssaha='m'
		elif userinput=='y':
			pools[count].runssaha='p'
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
	mindepth=float(quality)/5
else:
	mindepth=float(quality)



#Loading reference sequence
if not os.path.isfile(ref):
	print "\nError: Cannot open reference fasta file "+ref
	sys.exit()
reflines=open(ref, "rU").readlines()
refseq=''
for line in reflines:
	if line.strip()[0]!='>':
		refseq=refseq+line.strip().upper()


#Running Ssaha where required

for pool in pools:
	if runssaha=='y' or pool.runssaha=='m':
		print "\nRunning Ssaha on "+pool.name+'...'
		sys.stdout.flush()
		#Ssaha commands.
		print "ssaha2 -rtype "+rtype+" "+ref+" fastqs/"+pool.name+".fastq > "+pool.name+"/cigar.tmp"
		print "here1"
		sys.stdout.flush()
		os.system("ssaha2 -rtype "+rtype+" -output cigar "+ref+" fastqs/"+pool.name+".fastq > "+pool.name+"/cigar.tmp")
		print "here2"
		sys.stdout.flush()
		os.system('grep "^cigar" '+pool.name+"/cigar.tmp > "+pool.name+"/cigar2.tmp")
		print "here3"
		sys.stdout.flush()
		os.system("mv "+pool.name+"/cigar2.tmp "+pool.name+"/cigar.tmp")
		print "here4"
		sys.stdout.flush()
		os.system("/nfs/team81/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_cigar "+pool.name+"/cigar.tmp "+pool.name+"/cigar2.tmp")
		print "here5"
		sys.stdout.flush()
		os.system("awk '{print $2}' "+pool.name+"/cigar2.tmp > "+pool.name+"/readnames.tmp")
		print "here6"
		sys.stdout.flush()
		os.system("/nfs/team81/tdo/bin/pileup_v0.4/ssaha_pileup/other_codes/get_seqreads/get_seqreads "+pool.name+"/readnames.tmp fastqs/"+pool.name+".fastq "+pool.name+"/fastq.tmp")
		print "here7"
		sys.stdout.flush()
		if rtype=='454':
			os.system("/nfs/team81/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 0 -trans 0 "+pool.name+"/cigar2.tmp "+ref+" "+pool.name+"/fastq.tmp > "+pool.name+"/snp.tmp")
			os.system("egrep ^SNP "+pool.name+"/snp.tmp > "+pool.name+"/all.snp")
			os.system("/nfs/team81/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 0 -trans 0 "+pool.name+"/cigar2.tmp "+ref+" "+pool.name+"/fastq.tmp > "+pool.name+"/pileup.tmp")
		elif rtype=='solexa':
			os.system("/nfs/team81/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -solexa 1 -trans 0 "+pool.name+"/cigar2.tmp "+ref+" "+pool.name+"/fastq.tmp > "+pool.name+"/snp.tmp")
			print "here8"
			sys.stdout.flush()
			os.system("egrep ^SNP "+pool.name+"/snp.tmp > "+pool.name+"/all.snp")
			print "here9"
			sys.stdout.flush()
			os.system("/nfs/team81/tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_pileup -cons 1 -solexa 1 -trans 0 "+pool.name+"/cigar2.tmp "+ref+" "+pool.name+"/fastq.tmp > "+pool.name+"/pileup.tmp")
		os.system("egrep ^cons "+pool.name+"/pileup.tmp > "+pool.name+"/all.pileup")
		#os.system("rm -f "+pool.name+"/pileup.tmp "+pool.name+"/snp.tmp "+pool.name+"/cigar.tmp "+pool.name+"/cigar2.tmp "+pool.name+"/fastq.tmp "+pool.name+"/readnames.tmp")

sys.exit()

reflen=len(refseq)
reflennogaps=len(refseq.replace('-',''))

if embl!='':
	print "\nReading EMBL file...",
	sys.stdout.flush()
	CDSposns={}
	CDSpools={}
	embldata=[]
	donelen=0
	CDScount=0		
	for i in range(reflennogaps):
		embldata.append('I')

	lines=open(embl,'rU').readlines()
	
	for x, line in enumerate(lines):
		words=line.strip().split()
		if len(words)>2 and words[1] in ('CDS', 'rRNA', 'tRNA') and '..' in words[2]:
			#print words
			foundpool='n'
			y=x+1
			locustag=''
			while foundpool=='n' and lines[y].split()[1] not in ('CDS', 'rRNA', 'tRNA', 'gene', 'misc_feature', 'sig_peptide', 'repeat_unit', 'repeat_region', 'misc_RNA'):
				if '/systematic_id' in lines[y].split()[1] or '/locus_tag' in lines[y].split()[1]:
					locustag=lines[y].split()[1].split('"')[1]
					#print locustag
					foundpool='y'
				y=y+1
			
			direction='f'
			feature=words[1]
			if feature=='CDS':
				CDScount=CDScount+1
				CDSposns[CDScount]=[]
				CDSpools[CDScount]=locustag
			embldata, CDSposns=getCDSseq(embldata, words[2].replace('<','').replace('>',''), direction, CDScount, feature)
			

	
	snptoCDSpool={}
	for i in CDSposns.keys():
		for j in CDSposns[i]:
			snptoCDSpool[j[0]]=CDSpools[i]
		
	print "Done"
	sys.stdout.flush()


poolsort.sort()

avcoverage={}
covcount={}

#hassnps={}
allsnps={}
refbases={}
	


alssnpsummary={}


print "\nIdentifying SNPs..."
sys.stdout.flush()
for pool in pools:
	print pool.name+'...',
	sys.stdout.flush()
	lines=open(pool.name+'/all.snp', 'rU').readlines()
	
	tempsnps={}
	tempseq=""
	tempsnpsummary={'A':{'C':0, 'G':0, 'T':0}, 'C':{'A':0, 'G':0, 'T':0}, 'G':{'C':0, 'A':0, 'T':0}, 'T':{'C':0, 'G':0, 'A':0}}
	
	for line in lines:
		words=line.split()
		het_hez=words[0].split('_')[1].replace(':','')
		snplocation=int(words[3])
		refbase=words[5]
		cnsbase=words[6]
		snpscore=int(words[2])
		snpdepth=int(words[4])
		Big={'A':int(words[7])-int(words[13]), 'C':int(words[8])-int(words[13]),'G':int(words[9])-int(words[13]),'T':int(words[10])-int(words[13])}
		Small={'A':int(words[7]), 'C':int(words[8]), 'G':int(words[9]), 'T':int(words[10])}
		
		if snpscore<quality or snpdepth<depth or het_hez=='hez':
			continue
		
#		if cnsbase not in ['A','C','G','T']:
#			continue
#			y=0
#			bestsnp=[]
#			for x in Big.keys():
#				if Big[x]>y:
#					bestsnp=[x]
#					y=Big[x]
#				elif Big[x]==y:
#					bestsnp.append(x)
#			
#			if len(bestsnp)>1:
#				#continue
#				y=0
#				bestsnpb=[]
#				for x in bestsnp:
#					if Small[x]>y:
#						bestsnpb=[x]
#						y=Small[x]
#					elif Small[x]==y:
#						bestsnpb.append(x)
#				cnsbase=bestsnpb[0]
#			else:
#				cnsbase=bestsnp[0]
			
			
			
			
			
				
#			if cnsbase not in ['A','C','G','T'] or cnsbase==refbase:
#				continue
		
		if allsnps.has_key(snplocation):
			allsnps[snplocation][0]=allsnps[snplocation][0]+1
			allsnps[snplocation][2][cnsbase]=allsnps[snplocation][2][cnsbase]+1
		else:
			allsnps[snplocation]=[1,refbase,{'A':0, 'C':0, 'G':0, 'T':0}]
			allsnps[snplocation][2][cnsbase]=1
		tempsnpsummary[refbase][cnsbase]=tempsnpsummary[refbase][cnsbase]+1
		tempsnps[snplocation]=cnsbase
	
	
	snpcount=len(tempsnps.keys())
	
	lines=open(pool.name+'/all.pileup', 'rU').readlines()
	mapped=0.0
	
	#for x in hassnps.keys():
		#line=lines[x-1]
	for line in lines:
		
		words=line.split()
		snplocation=int(words[2])
		refbase=words[4]
		readdepth=int(words[3])
		#reads=words[4][1:]
		
		if avcoverage.has_key(snplocation):
			avcoverage[snplocation]=avcoverage[snplocation]+readdepth
		else:
			avcoverage[snplocation]=readdepth
			covcount[snplocation]=0

		if readdepth>depth:
			covcount[snplocation]=covcount[snplocation]+1
			mapped=mapped+1
			if not tempsnps.has_key(snplocation):
				tempsnps[snplocation]=refbase
				tempseq=tempseq+refbase
			else:
				tempseq=tempseq+tempsnps[snplocation]
		else:
			tempseq=tempseq+"-"
			if tempsnps.has_key(snplocation):
				snpcount=snpcount-1
				tempsnpsummary[refbase][tempsnps[snplocation]]=tempsnpsummary[refbase][tempsnps[snplocation]]-1
			tempsnps[snplocation]="-"
			if allsnps.has_key(snplocation):
				if allsnps[snplocation][0]==1:
					del allsnps[snplocation]
				else:
					allsnps[snplocation][0]=allsnps[snplocation][0]-1
					allsnps[snplocation][2][cnsbase]=allsnps[snplocation][2][cnsbase]-1
	
	pool.SNPs=snpcount
	pool.sequence=tempseq
	pool.mapped=tempsnps
	pool.snpsummary=tempsnpsummary
	pool.nummapped=mapped
	pool.percentmapped=((mapped/reflen)*100)
	print '%.2f%% of reference mapped, %d SNPs found' % (pool.percentmapped, pool.SNPs)
	sys.stdout.flush()
		



#testout=open("testing.out","w")
#
#print >> testout, ">ref\n"+refseq
#
#for strain in pools:
#	print >> testout, ">"+strain.name+"\n"+strain.sequence
#
#testout.close()



if embl!='' and ref!='':
	print "\nCalculating dN/dS values...\nConcatenating reference CDSs...",
	sys.stdout.flush()
	dndsout=open(outfile+'_dnds.out','w')
	
	sys.stdout.flush()
	#dnbydsstats={}
	
	comp={'A':'T','T':'A','G':'C','C':'G'}

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
	for pool in pools:
		if pool.name==ref:
			continue
		
		print pool.name+'...',
		sys.stdout.flush()
		
		tmpCDSseq=''
		
		for key in CDSposns.keys():
			for x in CDSposns[key]:
				if x[1]=='f' and pool.sequence[x[0]-1]!='-':
					tmpCDSseq=tmpCDSseq+pool.sequence[x[0]-1].upper()
				elif pool.sequence[x[0]-1]!='-':
					tmpCDSseq=tmpCDSseq+comp[pool.sequence[x[0]-1].upper()]
				else:
					tmpCDSseq=tmpCDSseq+'-'
		
		
		
		pool.dNdSstats, pool.snptypes=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers)
		
		#N, S, dN, dS, pN, pS, varianceS, varianceN, z, (len(CDS)-gapcount), Nd, Sd
		if pool.dNdSstats['dS']!=0:
			print "%.2f" % (pool.dNdSstats['dN']/pool.dNdSstats['dS'])
			print >> dndsout, '\t'+str(pool.dNdSstats['dN']/pool.dNdSstats['dS']),
			sys.stdout.flush()
		
		
		
	#Testing doing each gene seperately
	
	#for sequence in sequences.keys():
	#		pool.sequence=''
	#		for i, j in enumerate(sequences[ref]):
	#			if j!='-':
	#				pool.sequence=pool.sequence+sequences[sequence][i]
	print "\nCalculating CDS dN/dS values...",
	for key in CDSposns.keys():
	
		refCDSseq=''
		
		for x in CDSposns[key]:
			if x[1]=='f':
				refCDSseq=refCDSseq+refseq[x[0]-1].upper()
			else:
				refCDSseq=refCDSseq+comp[refseq[x[0]-1].upper()]
		
		print >> dndsout, '\n'+CDSpools[key],
		#print CDSpools[key]
		sys.stdout.flush()
		
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
			dnbydsstats, gene_snptypes=dnbyds(refCDSseq, tmpCDSseq, CDSbasenumbers)
			
			if dnbydsstats['dS']!=0:
				print >> dndsout, '\t'+str(dnbydsstats['dN']/dnbydsstats['dS']),
			else:
				print >> dndsout, '\t-',
			
	
	
	dndsout.close()
		
		
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
		if base not in foundbases.keys():
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
			outstring=outstring+'\t'+embldata[j]+'\t'+snptoCDSpool[j]
		snptype='-'
		for pool in pools:
			if pool.name==ref:
				continue
			if pool.snptypes.has_key(j):
				if pool.snptypes[j]=='S':
					snpcolour='3'
				elif pool.snptypes[j]=='N':
					snpcolour='2'
				if snptype=='-':
					snptype=pool.snptypes[j]
				elif pool.snptypes[j] not in snptype:
					snptype=snptype+'/'+pool.snptypes[j]
		outstring=outstring+'\t'+snptype
	
	
	outstring=outstring+'\t'+refseq[j]
	refbase=refseq[j]
		
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
					if pool.snptypes.has_key(j):
						if pool.snptypes[j]=='S':
							tabstring=tabstring+'(synonymous) '
						elif pool.snptypes[j]=='N':
							tabstring=tabstring+'(non-synonymous) '
						elif pool.snptypes[j]=='1':
							tabstring=tabstring+'(gap in SNP codon) '
						elif pool.snptypes[j]=='2':
							tabstring=tabstring+'(SNP codon is STOP) '
						elif pool.snptypes[j]=='3':
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
	print >> output, outstring+'\t'+pattern
	if tabfile=='y':
		print >> tabout, tabstring+'"\nFT                   /colour='+snpcolour

	if graphs=='y':
		print >> covgraph, str(float(avcoverage[poscount])/len(pools))
		print >> covcountgraph, covcount[poscount]
		poscount=poscount+1


if graphs=='y':
	while reflen>=poscount:
		print >> snpgraph, '0'
		print >> covgraph, str(float(avcoverage[poscount])/len(pools))
		print >> covcountgraph, str(covcount[poscount])
		poscount=poscount+1


if align=='y':
	
	inaln=[]
	for pool in snpsequence.keys():
		if len(snpsequence[pool].replace('-',''))>len(refseq)/2:
			inaln.append(pool)
	
	count=0
	
	print >> alnfile, len(snpsequence.keys()), len(refseq)
	print >> alnfile, ref+' '*(10-len(ref))+refseq
	for pool in inaln:
		print >> alnfile, pool.replace('pool','').replace('.snps','')+' '*(10-len(pool))+snpsequence[pool]

				
	alnfile.close()



#Create tab file for each sequence vs ref (extended to length of total alignment to allow insertions)
#print refconvertb

#if tabfile=='y':
#	extrefout=open('Extended_'+ref+'.fna', 'w')
#	print >> extrefout, '>'+ref
#	count=0
#	outstring=''
#	for x in range(len(refseq)):
#		outstring=outstring+refseq[x]
#		count=count+1
#		if count==60:
#			outstring=outstring+'\n'
#			count=0
#	print >> extrefout, outstring
#	outstring=''
#	#print reflen, len(sequences[ref]), len(sequences[ref].replace('-',''))
#	extrefout.close()
#	if embl!='':
#		lines=open(embl,'rU').readlines()
#		extrefout=open('Extended_'+ref+'.embl', 'w')
#		for line in lines:
#			if len(line.strip().split())==3 and line.strip().split()[1]=='source':
#				#print '..'+str(len(sequences[ref])), '..'+str(len(sequences[ref].replace('-','')))
#				print >> extrefout, line.replace( '..'+str(len(refseq.replace('-',''))), '..'+str(len(refseq[ref]))),
#			elif len(line.strip().split())==3 and len(line.strip().split()[2].split('..'))==2:
#				positions=line.strip().split()[2]
#				#print positions
#				num_list = re.findall(r'.[0-9]+.', positions)
#				#print num_list
#				for numa in num_list:
#				
#					num = re.findall(r'[0-9]+', numa)[0]
#					
#					#print numa, num,
#					
#					numb=numa.replace(num,str(int(num)))
#					#print numb
#					
#					positions=positions.replace(numa, numb)
#					 
#					#print num, num[0]+str(refconvertb[int(num[1:-1])-1]+1)+num[-1], refconvertb[int(num[1:-1])-1]
#				
#				
#				#print positions
#				#start=refconvertb[int(line.strip().split()[2].split('..')[0])-1]
#				#end=refconvertb[int(line.strip().split()[2].split('..')[1])-1]
#				#posn=line[:22]+str(start+1)+'..'+str(end+1)
#				print >> extrefout, line[:21]+positions
#			elif len(line.strip().split())>1 and line.strip().split()[0]=='SQ':
#				A=0
#				C=0
#				G=0
#				T=0
#				Others=0
#				
#				for x in sequences[ref]:
#					if x=='A':
#						A=A+1
#					elif x=='C':
#						C=C+1
#					elif x=='G':
#						G=G+1
#					elif x=='T':
#						T=T+1
#					else:
#						Others=Others+1
#				print >> extrefout, "SQ   Sequence "+str(len(sequences[ref]))+" BP; "+str(A)+" A; "+str(C)+" C; "+str(G)+" G; "+str(T)+" T; "+str(Others)+" other;"
#				
#				count=0
#				countb=0
#				totalsofar=0
#				outstring=''
#				for x in range(len(sequences[ref])):
#					if count==0:
#						outstring=outstring+'     '
#					outstring=outstring+sequences[ref][x].lower()
#					count=count+1
#					totalsofar=totalsofar+1
#					countb=countb+1
#					
#					if count==60:
#						outstring=outstring+'       '+str(totalsofar)+'\n'
#						count=0
#						countb=0
#					if countb==10:
#						outstring=outstring+' '
#						countb=0
#				
#				if count!=0:
#					while count!=0:
#						outstring=outstring+' '
#						count=count+1
#						countb=countb+1
#						if count==60:
#							outstring=outstring+'       '+str(totalsofar)+'\n//'
#							count=0
#							countb=0
#						if countb==10:
#							outstring=outstring+' '
#							countb=0
#
#				
#				print >> extrefout, outstring,
#				outstring=''
#				
#				break
#				
#				
#			else:
#				print >> extrefout, line.strip()
#		extrefout.close()
#	for pool in pools:
#		if pool.name==ref:
#			continue
#		indtabfile=open(pool.name+'.tab', 'w')
#		print >> indtabfile, 'ID   SNP'
#		delstart=0
#		insstart=0
#		for i,j in enumerate(pool.sequence):
#			if j=='-' and j!=refseq[i]:
#				if i==0 or pool.sequence[i-1]!='-' or refseq[i-1]=='-':
#					delstart=i
#				#if (i+1)==len(sequences[ref]) or sequences[pool][i+1]!='-' or sequences[ref][i+1]=='-':
#					#if delstart!=i:
#					#	print >> indtabfile, 'FT   misc_feature    '+str(delstart+1)+'..'+str(i+1)
#					#else:
#					#	print >> indtabfile, 'FT   misc_feature    '+str(i+1)
#				#	print >> indtabfile, 'FT                   /colour=3'#Deletion = green
#					
#			elif j!=refseq[i]:
#				if i==0 or pool.sequence[i-1]=='-':
#					insstart=i
#				#if (i+1)==len(sequences[ref]) or sequences[ref][i+1]!='-' or sequences[sequence][i+1]=='-':
#					#if insstart!=i:
#					#	print >> indtabfile, 'FT   misc_feature    '+str(insstart+1)+'..'+str(i+1)
#					#else:
#					#	print >> indtabfile, 'FT   misc_feature    '+str(i+1)
#					#print >> indtabfile, 'FT                   /colour=10'#Insertion = orange
#					
#			elif j!=refseq[i]:
#				print >> indtabfile, 'FT   misc_feature    '+str(i+1)+'..'+str(i+2)
#				print >> indtabfile, 'FT                   /refBase="'+refseq[i]+'"'
#				#tabstring=tabstring+'FT                   /SNPAllele="'+snpbase+'"\n'
#				print >> indtabfile, 'FT                   /SNPBase="'+j+'"'
#				if pool.snptypes.has_key(i):
#					if pool.snptypes[i]=='S':
#						print >> indtabfile, 'FT                   /colour=2'#Synonymous=red
#					elif pool.snptypes[refconvert[i]]=='N':
#						print >> indtabfile, 'FT                   /colour=4'#Non-synonymous=blue
#					elif pool.snptypes[refconvert[i]]=='1':
#						print >> indtabfile, 'FT                   /colour=13'#add synonymous stuff here
#					else:
#						print >> indtabfile, 'FT                   /colour=13'#add synonymous stuff here
#				else:
#					print >> indtabfile, 'FT                   /colour=13'#add synonymous stuff here
#		
#		indtabfile.close()
#
#
#
#
#
#
#
#






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


#
#
#
#sys.exit()
#summary='y'
#if summary=='y':
#	summaryfile=open(outfile+'_summary.out','w')
#	
#	print >> summaryfile, 'Strain\tdN/dS',
#	
#	for x in ['A', 'C', 'G', 'T']:
#		for y in ['A', 'C', 'G', 'T']:
#			if x!=y:
#				print >> summaryfile, '\t'+x+'->'+y,
#	
#	print >> summaryfile, '\tTotal SNPs\t% Mapped Intragenic sites with SNP\t% Mapped Intergenic sites with SNP'
#	
#	for pool in sequences.keys():
#		intragenic=0
#		intragenic_mapped=0
#		if pool.dS!=0:
#			print >> summaryfile, pool+'\t'+str(pool.dN/pool.dS),
#		else:
#			print >> summaryfile, pool+'\t'+str(pool.percentmapped)+'\t-',
#		for x in ['A', 'C', 'G', 'T']:
#			for y in ['A', 'C', 'G', 'T']:
#				if x!=y:
#					print >> summaryfile, '\t'+str(pool.snpsummary[x][y]),
#					count=count+pool.snpsummary[x][y]
#		
#		for CDS in CDSposns.keys():
#			for y in CDSposns[CDS]:
#				if pool.SNPs.has_key(y[0]):
#					intragenic=intragenic+1
#				for chromosome in pool.mapped.keys():
#					if pool.mapped[chromosome].has_key(y[0]):
#						intragenic_mapped=intragenic_mapped+1
#		intergenic_mapped=pool.nummapped-intragenic_mapped
#		totalSNPs=len(pool.SNPs.keys())
#		print intragenic_mapped, intergenic_mapped, pool.nummapped
#		
#		if intergenic_mapped>0:
#			intergenic=(float(totalSNPs-intragenic)/intergenic_mapped)*100
#		else:
#			intergenic=0
#		if intragenic_mapped>0:
#			intragenic=(float(intragenic)/intragenic_mapped)*100
#		else:
#			intragenic=0
#		
#		
#		print >> summaryfile, '\t'+str(totalSNPs)+'\t'+str(intragenic)+'\t'+str(intergenic)
#	
#	print >> summaryfile, '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(len(snpsort))
#	
#	summaryfile.close()
#	print 'Done'

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















#
#
#
#
#
#
#
#
#poscount=1
#
#for key in allsnps.keys():
#	
#	
#
#	snpsort=allsnps[key].keys()
#	snpsort.sort()
#
#	for j in snpsort:
#	
#		if graphs=='y':
#			while j>poscount:
#				print >> snpgraph, '0'
#				print >> covgraph, str(float(avcoverage[poscount])/len(pools))
#				print >> covcountgraph, str(covcount[poscount])
#				poscount=poscount+1
#		
#		
#		outstring=''
#		tabstring=''
#		
#		if len(allsnps.keys())>1:
#			outstring=outstring+key+'\t'+str(j)
#		else:
#			outstring=outstring+str(j)
#		if embl!='':
#                        outstring=outstring+'\t'+embldata[j]
#
#		outstring=outstring+'\t'+allsnps[key][j][1]
#
#		refsequence=refsequence+allsnps[key][j][1]
#		snpbase=''
#		snpbasecount=0
#		for base in allsnps[key][j][2].keys():
#			if allsnps[key][j][2][base]>snpbasecount:
#				snpbasecount=allsnps[key][j][2][base]
#				snpbase=base
#			elif allsnps[key][j][2][base]==snpbasecount:
#				snpbase=snpbase+','+base
#	
#		outstring=outstring+'\t'+snpbase+'\t'+str(allsnps[key][j][0])
#		
#		if tabfile=='y':
#                	tabstring=tabstring+'FT   SNP             '+str(j)+'\n'
#			tabstring=tabstring+'FT                   /refAllele="'+allsnps[key][j][1]+'"\n'
#			#tabstring=tabstring+'FT                   /SNPAllele="'+snpbase+'"\n'
#			tabstring=tabstring+'FT                   /SNPstrains="'
#
#		
#
#		pattern=''
#		numsnpbases=1
#		snpbases=[]
#		for pool in pools:
#			
#			if pool.mapped[key].has_key(j):# and pool.SNPs[key][j][1]>8:# and pool.SNPs[key][j][2]>(pool.SNPs[key][j][3]*4):
#				if pool.mapped[key][j]!=allsnps[key][j][1]:
#					
#					if tabfile=='y':
#						tabstring=tabstring+pool.name+'='+pool.mapped[key][j]+' '
#					pattern=pattern+'1'
#					outstring=outstring+'\t1'#+str(pool.mapped[key][j][1])
#					if pool.mapped[key][j] not in snpbases:
#						snpbases.append(pool.mapped[key][j])
#						numsnpbases=numsnpbases+1
#					pool.SNPs[j]=''
#					
#				else:
#					pattern=pattern+'0'
#					outstring=outstring+'\t0'
#				sequence[pool.name]=sequence[pool.name]+pool.mapped[key][j]
#			else:
#				pattern=pattern+'-'
#				sequence[pool.name]=sequence[pool.name]+'?'
#				outstring=outstring+'\t-'
#		
#		if numsnpbases>1:
#			if graphs=='y':
#				print >> snpgraph, '1'
#			print >> output, outstring+'\t'+pattern
#			if tabfile=='y':
#                        	print >> tabout, tabstring+'"\nFT                   /colour='+str(numsnpbases)
#
#		else:
#			if graphs=='y':
#				print >> snpgraph, '0'
#			refsequence=refsequence[:-1]
#			for pool in pools:
#				sequence[pool.name]=sequence[pool.name][:-1]
#		if graphs=='y':
#			print >> covgraph, str(float(avcoverage[poscount])/len(pools))
#			print >> covcountgraph, covcount[poscount]
#			poscount=poscount+1
#
#
#if graphs=='y':
#	while reflen>=poscount:
#		print >> snpgraph, '0'
#		print >> covgraph, str(float(avcoverage[poscount])/len(pools))
#		print >> covcountgraph, str(covcount[poscount])
#		poscount=poscount+1
#
#
#if align=='y':
#	
#	count=0
#	for pool in pools:
#		if pool.percentmapped>50:
#			count=count+1
#	
#	print >> alnfile, count+1, len(refsequence)
#	print >> alnfile, 'Reference\t'+refsequence
#		
#	for pool in pools:
#		if pool.percentmapped>50:
#			print >> alnfile, pool.name.replace('pool','').replace('.snps','')+'\t'+sequence[pool.name]
#
#				
#	alnfile.close()
#
#
#
#
#summaryfile=open(outfile+'_summary.out','w')
#
#print >> summaryfile, 'Strain\tPercent of Reference Mapped\tdN/dS',
#
#for x in ['A', 'C', 'G', 'T']:
#	for y in ['A', 'C', 'G', 'T']:
#		if x!=y:
#			print >> summaryfile, '\t'+x+'->'+y,
#
#print >> summaryfile, '\tTotal SNPs',
#
#if embl!='':
#	print >> summaryfile, '\t% Mapped Intragenic sites with SNP\t% Mapped Intergenic sites with SNP'
#else:
#	print >> summaryfile, '\n',
#
#for pool in pools:
#	intragenic=0
#	intragenic_mapped=0
#	if pool.dS!=0:
#		print >> summaryfile, pool.name+'\t'+str(pool.percentmapped)+'\t'+str(pool.dN/pool.dS),
#	else:
#		print >> summaryfile, pool.name+'\t'+str(pool.percentmapped)+'\t-',
#	for x in ['A', 'C', 'G', 'T']:
#		for y in ['A', 'C', 'G', 'T']:
#			if x!=y:
#				print >> summaryfile, '\t'+str(pool.snpsummary[x][y]),
#				count=count+pool.snpsummary[x][y]
#	
#	totalSNPs=len(pool.SNPs.keys())
#	
#	if embl!='':
#		for CDS in CDSposns.keys():
#			for y in CDSposns[CDS]:
#				if pool.SNPs.has_key(y[0]):
#					intragenic=intragenic+1
#				for chromosome in pool.mapped.keys():
#					if pool.mapped[chromosome].has_key(y[0]):
#						intragenic_mapped=intragenic_mapped+1
#		intergenic_mapped=pool.nummapped-intragenic_mapped
#	
#	#print intragenic_mapped, intergenic_mapped, pool.nummapped
#	
#	
#		if intergenic_mapped>0:
#			intergenic=(float(totalSNPs-intragenic)/intergenic_mapped)*100
#		else:
#			intergenic=0
#		if intragenic_mapped>0:
#			intragenic=(float(intragenic)/intragenic_mapped)*100
#		else:
#			intragenic=0
#	
#	print >> summaryfile, '\t'+str(totalSNPs),
#	
#	if embl!='':
#		print >> summaryfile, '\t'+str(intragenic)+'\t'+str(intergenic)
#
#if embl!='':
#	print >> summaryfile, '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(len(snpsort))
#else:
#	print >> summaryfile, '\t\t\t\t\t\t\t\t\t\t\t\t\t'+str(len(snpsort))
#
#
#print 'Done'
#
#
#if raxml=='y':
#	outlen=0
#	userinput='x'
#	if '/' in outfile:
#		outlen=-1*len(outfile.split('/')[-1])
#	if os.path.isfile('RAxML_info.'+outfile.split('/')[-1]):
#		print '\nRAxML files with extension '+outfile.split('/')[-1]+' already exist!'
#		while userinput not in ['y','n']:
#			userinput=raw_input('Overwrite? (y/n): ')
#			userinput.lower()
#		if userinput=='y':
#			os.system('rm RAxML_*'+outfile.split('/')[-1])
#	print "Running RAxML phylogeny with "+model+" model of evolution and "+str(bootstrap)+" bootstrap replicates..."
#	os.system("RAxML -f a -x "+str(random.randrange(1,99999))+" -p "+str(random.randrange(1,99999))+" -# "+str(bootstrap)+" -m "+model+" -s "+outfile+".aln -n "+outfile.split('/')[-1])
#	
#if tabfile=='y':
#	tabout.close()
#
#if graphs=='y':
#	snpgraph.close()
#	covgraph.close()
#	covcountgraph.close()
#output.close()
#summaryfile.close()
#
#
#print "All finished.\n"
