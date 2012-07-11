from Bio.Align import AlignInfo
from Bio.Align.Generic import Alignment
import sys, os
from Bio.Data import CodonTable
import math

geneticcode={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

start_codons=["ATG", "GTG", "TTG", "ATT", "CTG", "ATC", "ATA"]
stop_codons=["TGA", "TAA", "TAG"]

snptypelong={"S":"Synonymous", "N":"Nonsynonymous", "P":"Pseudogene", "M":"STIP (stop codon to non-stop codon)", "O":"SNOP (non-stop codon to stop condon)", "*":"STOP (stop codon to stop codon)", "I":"Intergenic"}

snptype={"S":"Synonymous", "N":"Nonsynonymous", "P":"Pseudogene", "M":"STIP", "O":"SNOP", "*":"STOP", "I":"Intergenic"}

snptype_to_colours={"S":"3", "N":"2", "P":"11", "M":"4", "O":"5", "*":"6", "I":"1"}


def translate(sequence):
	
	aa_sequence=""
	
	for i in range(0,len(sequence),3):
		if geneticcode.has_key(sequence[i:i+3].upper()):
			aa_sequence=aa_sequence+geneticcode[sequence[i:i+3].upper()]
		else:
			aa_sequence=aa_sequence+"X"
	
	
	return aa_sequence



#############################################
# Function to reverse complement a sequence #
#############################################

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



############################################
# Function to identify snps from alignment #
############################################

#def snp_locations_from_alignment(alignmentObject, unknowns=[ "?", "N", "X"]):
#
#	SNPlocations=[]
#	
#	#summary_align = AlignInfo.SummaryInfo(alignmentObject)
#	
#	print "Identifying SNP locations"
#	sys.stdout.flush()
#	
#	
#	for columnnumber in range(alignmentObject.get_alignment_length()):
#		#when we get the new version of biopython, this can be changed to slice notation alignmentObject[:,base]
#		column=alignmentObject.get_column(columnnumber).upper()
#		
#		#print columnnumber
#		#print column
#		column=column.replace("N","").replace("?","").replace("X","")
#
#		
#		firstbase=column[0]
#		for base in column[1:]:
#			if base!=firstbase:
#				SNPlocations.append(columnnumber)
#				break
#	
#
#	print "Found", len(SNPlocations)
#	return SNPlocations






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

def dnbyds(CDS, SNPseq):
	
	#geneticcode={'TTT':'Phe', 'TTC':'Phe', 'TTA':'Leu', 'TTG':'Leu', 'TCT': 'Ser', 'TCC': 'Ser','TCA': 'Ser','TCG': 'Ser', 'TAT': 'Tyr','TAC': 'Tyr', 'TAA': 'STOP', 'TAG': 'STOP', 'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp', 'CTT': 'Leu','CTC': 'Leu','CTA': 'Leu','CTG': 'Leu', 'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'CAT': 'His', 'CAC': 'His', 'CAA': 'Gin', 'CAG': 'Gin', 'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met', 'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val', 'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu', 'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'}
	#geneticcode={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
	
	codonsynonyms={}
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
	
	#print CDS, SNPseq
	

	if ((float(len(SNPseq.replace("-",""))))/3)-(math.floor(float(len(SNPseq.replace("-","")))/3))!=0:
		#print "frameshift in daughter sequence"
		raise ValueError()
	if ((float(len(CDS.replace("-",""))))/3)-(math.floor(float(len(CDS.replace("-","")))/3))!=0:
		#print "frameshift in parent sequence"
		raise ZeroDivisionError()


	if len(CDS)!=len(SNPseq):
		#print "Error: sequences must be the same length to calculate dN/dS!"
		raise FloatingPointError
		return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}, SNPtype, AAtype

	
	#create a codonsynonyms dictionary
	
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
	


	

	
	for x in range(0,len(CDS),3):
		if x+3>len(CDS):
			break
		numcodons=numcodons+1
		codon=CDS[x:x+3]
		SNPcodon=SNPseq[x:x+3]
		
#		if codon!=SNPcodon:
#			print codon, SNPcodon
		
#		codonposn=[CDSbasenumbers[x], CDSbasenumbers[x+1], CDSbasenumbers[x+2]]
#		newSNPcodon='---'
		
		
		
	#	if "-" in codon and len(codon.replace("-","").replace("?","").replace("N",""))>0:
			#print "frameshift in parent sequence"
		#	raise ZeroDivisionError()
#		else:
#			gapcount=gapcount+3
#			continue
		#if "-" in SNPcodon and len(SNPcodon.replace("-","").replace("-","").replace("N",""))>0:
			#print "frameshift in daughter sequence"
			#raise ValueError()

		
		if '-' in codon or '-' in SNPcodon or 'N' in codon or 'N' in SNPcodon or '?' in codon or '?' in SNPcodon or geneticcode[codon]=='*' or geneticcode[SNPcodon]=='*':
		
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
		return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}
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

	return {'N':N, 'S':S, 'dN':dN, 'dS':dS, 'pN':pN, 'pS':pS, 'varS':varianceS, 'varN':varianceN, 'z':z, 'CDSlen':(len(CDS)-gapcount), 'Nd':Nd, 'Sd':Sd}






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



#################################################################################
# Function to return a strict consensus sequence when given an alignment object #
#################################################################################

def consensus_from_alignment(alignmentObject, unknowns=["?", "N", "X"]):
	
	consensus_sequence=""
	
	print "Creating consensus..."
	sys.stdout.flush()
	
	count=0
	total=0.0
	hundredth=alignmentObject.get_alignment_length()/100
	
	
	for x in range(alignmentObject.get_alignment_length()):

		count=count+1
		if count>=hundredth:
			total=total+count
			count=0
			print "%.0f%% complete\r" % (100*(total/alignmentObject.get_alignment_length())),
			sys.stdout.flush()

		
		foundbases=[]
		for record in alignmentObject:
			base=record.seq[x].upper()
			if base not in foundbases and base not in unknowns:# 
				foundbases.append(base)
			if len(foundbases)>1:
				consensus_sequence=consensus_sequence+"N"
				break
		if len(foundbases)==1:
			consensus_sequence=consensus_sequence+foundbases[0]
		elif len(foundbases)==0:
			consensus_sequence=consensus_sequence+"N"
	

	#print "100.00% complete"#Found %d SNP locations" % len(SNPlocations),
	
	
	print "100% complete"
	sys.stdout.flush()
	
	return consensus_sequence




#############################################################################
# Function to return a list of snp locations when given an alignment object #
#############################################################################


def snp_locations_from_alignment(alignmentObject, unknowns=set(["?", "N", "X"]), incgaps=True):
	
	SNPlocations=[]
	consensus_sequence=""
	
	#summary_align = AlignInfo.SummaryInfo(alignmentObject)
	
	#print summary_align
	
	
	print "Identifying SNP locations..."
	sys.stdout.flush()
	
	count=0
	total=0.0
	hundredth=alignmentObject.get_alignment_length()/100
	
	
	for x in xrange(alignmentObject.get_alignment_length()):

		count=count+1
		if count>=hundredth:
			total=total+count
			count=0
			print "%.0f%% complete\r" % (100*(total/alignmentObject.get_alignment_length())),
			sys.stdout.flush()

		
		foundbases=[]
		for record in alignmentObject:
			base=record.seq[x].upper()
			
			if base not in foundbases and base not in unknowns:# 
				foundbases.append(base)
			if len(foundbases)>1:
				SNPlocations.append(x)
				consensus_sequence=consensus_sequence+"N"
				break
		if len(foundbases)==1:
			if foundbases[0]=="-" and incgaps:
				SNPlocations.append(x)
			consensus_sequence=consensus_sequence+foundbases[0]
		elif len(foundbases)==0:
			consensus_sequence=consensus_sequence+"N"
	
		
	#print "100.00% complete"#Found %d SNP locations" % len(SNPlocations),
	
	
	if incgaps:
		outtext="locations with a SNP or gap"
	else:
		outtext="locations with a SNP"
	
	print "Found", len(SNPlocations), outtext
	sys.stdout.flush()
	
	return SNPlocations, consensus_sequence
	





#############################################################################
# Function to return a list of snp locations when given an alignment object #
#############################################################################


def snp_locations_from_alignment_fast(alignmentObject, unknowns=set(["?", "N", "X"]), incgaps=True):
	
	SNPlocations=[]
	consensus_sequence=""
	consensus_list=[]
	
	#summary_align = AlignInfo.SummaryInfo(alignmentObject)
	
	#print summary_align
	
	append=consensus_list.append
	print "Identifying SNP locations..."
	sys.stdout.flush()
	
	count=0
	total=0.0
	hundredth=alignmentObject.get_alignment_length()/100
	baseset=set(["A","C","G","T"])
	if incgaps:
		baseset.add("-")
	
	for x in xrange(alignmentObject.get_alignment_length()):
	
		count=count+1
		if count>=hundredth:
			total=total+count
			count=0
			print "%.0f%% complete\r" % (100*(total/alignmentObject.get_alignment_length())),
			sys.stdout.flush()
	
		column=set(alignmentObject.get_column(x).upper())
		foundbases=column.intersection(baseset)
		
		if len(foundbases)==1:
			if "-" in foundbases:
				SNPlocations.append(x)
			append(foundbases.pop())
		elif len(foundbases)==0:
			append("N")
		else:
			SNPlocations.append(x)
			append("N")
	
		
	print "100.00% complete"
	
	
	if incgaps:
		outtext="locations with a SNP or gap"
	else:
		outtext="locations with a SNP"
	
	print "Found", len(SNPlocations), len(consensus_list), alignmentObject.get_alignment_length(), outtext
	
	return SNPlocations, ''.join(consensus_list)





#############################################################################
# Function to return a list of snp locations when given an alignment object #
#############################################################################


def snp_locations_from_alignment_tmp(alignmentObject, unknowns=["?", "N", "X"]):
	
	SNPlocations=[]
	consensus_sequence=""
	
	#summary_align = AlignInfo.SummaryInfo(alignmentObject)
	
	#print summary_align
	
	
	print "Identifying SNP locations..."
	sys.stdout.flush()
	
	count=0
	total=0.0
	hundredth=50000/100
	
	
	
	
	for x in range(50000):

		count=count+1
		if count>=hundredth:
			total=total+count
			count=0
			print "%.0f%% complete\r" % (100*(total/50000)),
			sys.stdout.flush()

		
		foundbases=[]
		for record in alignmentObject:
			base=record.seq[x].upper()
			
			if base not in foundbases and base not in unknowns:# 
				foundbases.append(base)
			if len(foundbases)>1:
				SNPlocations.append(x)
				consensus_sequence=consensus_sequence+"N"
				break
		if len(foundbases)==1:
			consensus_sequence=consensus_sequence+foundbases[0]
		elif len(foundbases)==0:
			consensus_sequence=consensus_sequence+"N"
	

	#print "100.00% complete"#Found %d SNP locations" % len(SNPlocations),
	
	
	print "Found", len(SNPlocations), "locations with a SNP or gap"
	sys.stdout.flush()
	
	return SNPlocations, consensus_sequence


	
	
	
class SNP():
	def __init__(self):
	
		self.position=-1
		self.codon_position=-1
		#self.intergenic=False
		#self.position_in_reference=0
		self.CDSid=-1
		self.CDSname=""
		self.position_in_CDS=-1
		self.strand=0
		self.parent_base=""
		self.daughter_base=""
		self.parent_codon=""
		self.daughter_codon=""
		self.parent_aminoacid=""
		self.daughter_aminoacid=""
		self.homoplasy=False
		self.originalSNP=False
		self.oldhomoplasy=False
		self.homoplasies=[]
		self.oldhomoplasies=[]
		self.ambiguous=False
		self.parent=-1
		self.daughter=-1
		self.colour=0
		self.SNP_type="" #synonymous, nonsynonymous, intergenic
		self.codon_type="" #synonymous, nonsynonymous, intergenic, from_stop, to_stop, pseudogene
		self.recombination=False
		self.addrecombination=False
		
		
		
	def get_reference_position(self, reference_name, SeqRecordObject):
		
		print reference_name
	
	
	def write_tab_format(self, handle, strain_list=[], colourby="synonymous", rootbase="N", nodebase="N"):
		basecolours={"A": "2", "C": "3", "G": "10", "T":"4"}
		print >> handle, "FT   SNP             "+str(self.position+1)
		if self.strand==-1:
			print >> handle, 'FT                   /strand="reverse"'
		if self.parent!=-1 and self.daughter!=-1:
			print >>handle, 'FT                   /node="'+str(self.parent)+'->'+str(self.daughter)+'"'
		print >>handle, 'FT                   /SNP="'+', '.join(self.parent_base)+'->'+', '.join(self.daughter_base)+'"'
		if self.homoplasy:
		
			homoplasyline=[]
			for htype in self.homoplasies:
				if htype[0]=="c":
					homoplasyline.append("convergence with branch leading to "+str(htype[2]))
				elif htype[0]=="r":
					homoplasyline.append("reversal from branch leading to "+str(htype[2]))
				elif htype[0]=="d":
					homoplasyline.append("reversed in branch leading to "+str(htype[2]))
					
		
			print >>handle, 'FT                   /homoplasy="'+', '.join(homoplasyline)+'"'
		
		if self.recombination and len(self.recombination)>0:
			recombline=[]
			for key in self.recombination:
				recombline.append(str(key)+": "+str(self.recombination[key]))
			print >>handle, 'FT                   /recombination="'+', '.join(recombline)+'"'
		
		if self.codon_type!="":
	
			print >>handle, 'FT                   /codon_type="'+snptypelong[self.codon_type]+'"'
			
			if colourby=="synonymous":
				print >>handle, 'FT                   /colour='+snptype_to_colours[self.codon_type]
		if colourby=="homoplasy":
			if self.homoplasy:
				print >>handle, 'FT                   /colour=2'
			else:
				print >>handle, 'FT                   /colour=4'
		elif colourby=="base":
			
			if rootbase!=nodebase and nodebase in ["A", "C", "G", "T"] and rootbase in ["A", "C", "G", "T"]:
				print >> handle, 'FT                   /colour='+basecolours[nodebase]
		elif colourby=="homoplasy_bases":
			if self.homoplasy:
				
			
				if rootbase!=nodebase and nodebase in ["A", "C", "G", "T"] and rootbase in ["A", "C", "G", "T"]:
					print >> handle, 'FT                   /colour='+basecolours[nodebase]
			else:
				print >>handle, 'FT                   /colour=13'
		
			
		
		
#		else:
#			print >>handle, 'FT                   /colour=1'
		if len(strain_list)>0:
			print >> handle, 'FT                   /taxa="'+', '.join(strain_list)+'"'
		
		
	
	def display(self):
		print "Alignment position:", self.position
		if self.ambiguous:
			print "Ambiguous",
		print "SNP: from", ', '.join(self.parent_base), "to", ', '.join(self.daughter_base)
		if self.SNP_type!="":
			print "SNP type:", self.SNP_type
		if self.position_in_CDS!=-1:
			print "CDS:",	self.CDSname
			print "Strand:", self.strand
			print "CDS position:", self.position_in_CDS+1
		if self.codon_position!=-1:
			print "Codon position:", self.codon_position
		if self.parent_codon!="" and self.daughter_codon!="":
			print "Codon change: from", self.parent_codon, "("+self.parent_aminoacid+") to", self.daughter_codon, "("+self.daughter_aminoacid+")"
		if self.codon_type!="":
			print "Codon change type:", snptypelong[self.codon_type]
		
	

	def get_annotation_info(self, parent_annotation, daughter_annotation, parent_seq, daughter_seq, parent=-1, daughter=-1):
	
		if self.position==-1:
			return
		
		self.parent=str(parent)
		self.daughter=str(daughter)
		
		for x, feature in enumerate(daughter_annotation):
			if feature['location'][1]>self.position and feature['location'][0]<self.position:
				
				self.CDSname=feature['name']
				
				
				self.CDSid=x
				
				self.strand=feature['strand']
				
				if feature['strand']==1:
					if parent_annotation[x]['location'][1]>self.position and parent_annotation[x]['location'][0]<self.position:
						parent_CDS_seq=parent_seq[parent_annotation[x]['location'][0]:parent_annotation[x]['location'][1]]
					else:
						self.codon_type="P"
						break
					daughter_CDS_seq=daughter_seq[feature['location'][0]:feature['location'][1]]
				
					self.position_in_CDS=self.position-feature['location'][0]
					
					
					
					
				else:
					if parent_annotation[x]['location'][1]>self.position and parent_annotation[x]['location'][0]<self.position:
						parent_CDS_seq=revcomp(parent_seq[parent_annotation[x]['location'][0]:parent_annotation[x]['location'][1]])
					else:
						self.codon_type="P"
						break
					daughter_CDS_seq=revcomp(daughter_seq[feature['location'][0]:feature['location'][1]])
					
				
					self.position_in_CDS=((feature['location'][1]-1)-self.position)
					

				
				self.codon_position=(((float(self.position_in_CDS)/3)-math.floor(float(self.position_in_CDS)/3))*3)+1
				#print float(self.position_in_CDS)/3, self.codon_position
				
				#self.codon_position=int((((float(self.position_in_CDS)/3)-math.floor(float(self.position_in_CDS)/3))*3))+1
				
				self.parent_codon=str(parent_CDS_seq[int(math.floor(float(self.position_in_CDS))/3)*3:int((math.floor(float(self.position_in_CDS)/3)*3)+3)])
				
				if geneticcode.has_key(self.parent_codon):
					self.parent_aminoacid=geneticcode[self.parent_codon]
				else:
					self.parent_aminoacid="X"
					
				
				self.daughter_codon=str(daughter_CDS_seq[int(math.floor(float(self.position_in_CDS)/3)*3):int((math.floor(float(self.position_in_CDS)/3)*3)+3)])
				
				if geneticcode.has_key(self.daughter_codon):
					self.daughter_aminoacid=geneticcode[self.daughter_codon]
				else:
					self.daughter_aminoacid="X"
				
				
				if self.parent_aminoacid=="*" and self.daughter_aminoacid=="*":
					self.codon_type="*"
				elif self.parent_aminoacid=="*" and self.daughter_aminoacid!="*":
					self.codon_type="M"
				elif self.parent_aminoacid!="*" and self.daughter_aminoacid=="*":
					self.codon_type="O"
				elif self.parent_aminoacid==self.daughter_aminoacid:
					self.codon_type="S"
				else:
					self.codon_type="N"
				
				break
		
	
		if self.CDSid==-1:
			self.codon_type="I"
		self.colour=snptype_to_colours[self.codon_type]
	
		#self.display()









def find_gene_limits(sequence, predicted_start, predicted_end):
	#print sequence[:predicted_end]
	gap_list=set(["-"])
	missing_list=set(["?","X","N"])
	#codon=sequence[predicted_end-3:predicted_end]
	
	#print codon, translate(codon)

#	if codon in stop_codons:
#		return predicted_start, predicted_end
#		
#	may_be_stop=False
#	for stop_codon in stop_codons:
#		may_be_stop=True
#		for x, base in enumerate(codon):
#			if  base!=stop_codon[x] and base not in missing_list:
#				may_be_stop=False
#				break
#		
#		if may_be_stop:
#			#print "may be stop", codon
#			return predicted_start, predicted_end
	

	
	
	
	#first find start codon
	
#	if not sequence[predicted_start:predicted_start+3] in start_codons:
#		print "start codon has changed", sequence[predicted_start:predicted_start+3]
	
	start=predicted_start
	
	#print predicted_end
	
	#print predicted_end-predicted_start
	
	codon=""
	
	codon_bases=["","",""]
	end=start
	
	
	while codon not in stop_codons:
		#print codon, end, predicted_end-1
		if end>=predicted_end-1:
			for stop_codon in stop_codons:
				may_be_stop=True
				for x, base in enumerate(codon):
					if  base!=stop_codon[x] and base not in missing_list:
						may_be_stop=False
						#print codon, geneticcode[codon]
						break
				
				if may_be_stop:
					print "may be stop", codon, end, predicted_end
					return start, end
						
							
		codon_base=0
		
		while codon_base<3 and (end+(3-codon_base))<len(sequence):
		
			if sequence[end] not in gap_list:
				codon_bases[codon_base]=sequence[end]
				codon_base+=1
		

			end+=1
			
		codon=''.join(codon_bases)
	
		
		if (end+3)>=len(sequence):
			return start, end
	
	return start, end
	
	
	






def find_gene_limits_new(sequence, predicted_start, predicted_end):
	#print sequence[:predicted_end]
	gap_list=set(["-"])
	missing_list=set(["?","X","N"])
	
	start_codon_posn=-1
	stop_codon_posn=-1
	
	alternative_start_posns=[]
	alternative_end_posns=[]
	
	#codon=sequence[predicted_end-3:predicted_end]
	
	#print codon, translate(codon)

	
	start=predicted_start
	
	
	
	
	codon=""
	
	while codon not in stop_codons:
		codon=""
	
		codon_bases=["","",""]
		end=start
		codon_base=0
		while codon_base<3 and (end+(3-codon_base))<len(sequence):
			if sequence[end] not in gap_list:
				codon_bases[codon_base]=sequence[end]
				codon_base+=1
			end+=1
		
		codon=''.join(codon_bases)
		
		if len(codon)<3 or codon in stop_codons:
			break
#		print codon
		if start_codon_posn==-1:
		
			if codon in start_codons:
				start_codon_posn=start
			
			elif len(set(codon_bases).intersection(set(missing_list)))>0:
			#else:
				for start_codon in start_codons:
						may_be_start=True
						for x, base in enumerate(codon_bases):
							if  base!=start_codon[x] and base not in missing_list:
								may_be_start=False
#								print set(codon_bases).intersection(set(missing_list))
#								print codon_bases, codon, geneticcode[codon]
								break
						
						if may_be_start:
							alternative_start_posns.append(start)
							break
		
		
		if len(set(codon_bases).intersection(set(missing_list)))>0:
			for stop_codon in stop_codons:
					may_be_stop=True
					for x, base in enumerate(codon_bases):
						if  base!=stop_codon[x] and base not in missing_list:
							may_be_stop=False
#							print set(codon_bases).intersection(set(missing_list))
#							print codon_bases, codon, geneticcode[codon]
							break
					
					if may_be_stop:
						alternative_end_posns.append(end)
						break
		
		start=end
	
	if codon in stop_codons:
		stop_codon_posn=end
	
	
	#if start_codon_posn!=predicted_start or stop_codon_posn!=predicted_end:
	#if len(alternative_start_posns)>0 or len(alternative_end_posns)>0:
	#	print start_codon_posn, alternative_start_posns, stop_codon_posn, alternative_end_posns, predicted_start, predicted_end
	#	sys.stdout.flush()
	
	return start_codon_posn, alternative_start_posns, stop_codon_posn, alternative_end_posns
	
	
