from Bio.Align import AlignInfo
from Bio.Align.Generic import Alignment
import sys, os
from Bio.Data import CodonTable


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






######################################
# Function to create a SNP alignment #
######################################

def Create_SNP_alignment(alignment, SNPlocations):
	
	alphabet = Gapped(IUPAC.unambiguous_dna)

	SNPalignment = Alignment(alphabet)

	
	for record in alignment:
		SNPsequenceObject=""
		
		for base in SNPlocations:
			SNPsequenceObject=SNPsequenceObject+record.seq[base].replace("-","?")
		
		SNPalignment.add_sequence(record.id, SNPsequenceObject)
		
	
	return SNPalignment


#############################################################################
# Function to return a list of snp locations when given an alignment object #
#############################################################################


def snps_locations_from_alignment(alignmentObject, unknowns=["-", "?", "N", "X"]):
	
	SNPlocations=[]
	
	#summary_align = AlignInfo.SummaryInfo(alignmentObject)
	
	#print summary_align
	
	#return
	
	print "Identifying SNP locations..."
	sys.stdout.flush()
	
	count=0
	total=0.0
	
	
	
	for x in range(alignmentObject.get_alignment_length()):
#
#		count=count+1
#		if count==10000:
#			total=total+count
#			count=0
#			print "%.2f%% complete\r" % (100*(total/alignmentObject.get_alignment_length())),
#			sys.stdout.flush()

		
		foundbases=[]
		for record in alignmentObject:
			base=record.seq[x].upper()
			
			if base not in unknowns and base not in foundbases:# 
				foundbases.append(base)
			if len(foundbases)>1:
				SNPlocations.append(x)
				break
	

#	print "100.00% complete"#Found %d SNP locations" % len(SNPlocations),
#	sys.stdout.flush()
	return SNPlocations
	




#############################################################################
# Function to return a list of snp locations when given an alignment object #
#############################################################################


def snps_locations_from_alignment_fast(alignmentObject, unknowns=["-", "?", "N", "X"]):
	
	SNPlocations=[]
	
	#summary_align = AlignInfo.SummaryInfo(alignmentObject)
	
	#print summary_align
	
	#return
	
	print "Identifying SNP locations..."
	sys.stdout.flush()
	
	count=0
	total=0.0
	
	
	
	for x in range(alignmentObject.get_alignment_length()):
		column=alignmentObject.get_column(columnnumber)
		print column
		sys.exit()

		
		foundbases=[]
		for record in alignmentObject:
			base=record.seq[x].upper()
			
			if base not in unknowns and base not in foundbases:# 
				foundbases.append(base)
			if len(foundbases)>1:
				SNPlocations.append(x)
				break
	

#	print "100.00% complete"#Found %d SNP locations" % len(SNPlocations),
#	sys.stdout.flush()
	return SNPlocations











######################################################
# Function to identify if SNPs are synonymous or not #
######################################################


def SNP_types_using_reference(SNPlocations, reference, SeqRecordObject):
	
	geneticcode={'TTT':'Phe', 'TTC':'Phe', 'TTA':'Leu', 'TTG':'Leu', 'TCT': 'Ser', 'TCC': 'Ser','TCA': 'Ser','TCG': 'Ser', 'TAT': 'Tyr','TAC': 'Tyr', 'TAA': 'STOP', 'TAG': 'STOP', 'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp', 'CTT': 'Leu','CTC': 'Leu','CTA': 'Leu','CTG': 'Leu', 'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'CAT': 'His', 'CAC': 'His', 'CAA': 'Gin', 'CAG': 'Gin', 'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met', 'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val', 'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu', 'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'}
	
	


	print CodonTable




	
	
	
class SNP():
	def __init__(self):
	
		self.position_in_DNA_alignment=position
		self.intergenic=False
		self.position_in_reference
		self.CDSid
		self.position_in_CDS
		self.refbase
		self.SNPbase
		self.codon
		self.aminoacid
		self.SNP_type #synonymous, nonsynonymous, intergenic
		self.codontype #synonymous, nonsynonymous, intergenic, from_stop, to_stop, pseudogene
		
		
		
	
	
	
	
	