import sys
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import SeqIO
from Si_general import *
from Si_SeqIO import *
from Bio.Alphabet import IUPAC, Gapped
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Align.Generic import Alignment
from Bio.Align.Applications import MuscleCommandline, ClustalwCommandline
import subprocess
from random import *
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
import Si_SeqIO



def muscle_gene_realign(alignmentObject, annotationObject):

	newseqs={}
	
	for record in alignmentObject:
		newseqs[record.id]=""

	chars = string.ascii_letters + string.digits
	tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
	prevend=0

	for featurenumber, feature in enumerate(annotationObject.features):
		if feature.type=="CDS":
			
			start=feature.location.nofuzzy_start
			end=feature.location.nofuzzy_end
			
			
			
			if feature.strand==1:
				strand=1
				end = end-3
			else:
				strand=-1
				start=start+3
				
			while start<prevend:
				start=start+3
				
			
			
			for record in alignmentObject:
				newseqs[record.id]=newseqs[record.id]+str(record.seq.tomutable())[prevend:start]
				
			
			
			if feature.qualifiers.has_key("locus_tag"):
				print feature.qualifiers["locus_tag"][0]
				sys.stdout.flush()
			
			alphabet = Gapped(IUPAC.protein)
			
			#alphabet = Gapped(generic_protein)

			geneAlignment = Alignment(alphabet)
			
			nualigns={}
			
			pseudogene=False
			needsaligning=False
			
			for record in alignmentObject:
				
				if strand==1:
					nucseq=Seq(str(record.seq)[start:end])
				else:
					nucseq=Seq(str(Seq.reverse_complement(Seq(str(record.seq)[start:end]))))
				
				
				if "-" in str(nucseq):
					needsaligning=True
					nucseq=Seq(str(nucseq).replace("-",""))
				
				nualigns[record.id]=str(nucseq)
					
				aaseq=str(Seq.translate(nucseq))

				if "*" in aaseq:
					pseudogene=True
				geneAlignment.add_sequence(record.id, str(Seq(aaseq, alphabet)))
			
			if needsaligning:
				print "realigning"
				if pseudogene:
					print feature.qualifiers["locus_tag"][0], "contains stop codons in one or more sequences. Aligning as nucleotides"
					alphabet = Gapped(IUPAC.unambiguous_dna)
					geneAlignment = Alignment(alphabet)
					for newseq in nualigns.keys():
						#print newseq, len(nualigns[newseq])
						geneAlignment.add_sequence(newseq, nualigns[newseq])
				#print geneAlignment
				
				handle=open(tmpname+".fasta", "w")
				AlignIO.write([geneAlignment], handle, "fasta")
				handle.close()
				
				cline = MuscleCommandline(input=tmpname+".fasta", out=tmpname+".aln")
				
				#cline = ClustalwCommandline(infile=tmpname+".fasta", outfile=tmpname+".aln")
				#child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
				os.system(str(cline)+" >  /dev/null 2>&1")
				
				handle=open(tmpname+".aln", "rU")
				align = AlignIO.read(handle, "fasta")
				handle.close()
				
				if pseudogene:
					for record in align:
						if strand==1:
							newseqs[record.id]=newseqs[record.id]+str(record.seq)
						else:
							newseqs[record.id]=newseqs[record.id]+str(Seq.reverse_complement(record.seq))
				else:
					for record in align:
		#				print nualigns[record.id]
						newseq=""
						x=0
						for base in str(record.seq):
							if base=="-":
								newseq=newseq+"---"
							else:
								newseq=newseq+nualigns[record.id][x:x+3]
								x=x+3
						
						if strand==1:
							newseqs[record.id]=newseqs[record.id]+newseq
						else:
							newseqs[record.id]=newseqs[record.id]+str(Seq.reverse_complement(Seq(newseq)))
					
					
					
					
					#print record.seq
					#print newseq
					#newseqs[record.id]=newseqs[record.id]+str(record.seq.tomutable())[prevend:start]
				
				
				os.system("rm "+tmpname+".*")
				
				#print align
				
	#			if len(newseqs[record.id])>end:
	#				print end, len(newseqs[record.id]), record.id
				
				prevend=end
			else:
				prevend=start
			
#			if featurenumber>18:
#				break
			
	for record in alignmentObject:
		newseqs[record.id]=newseqs[record.id]+str(record.seq.tomutable())[prevend:]
	
	alphabet = Gapped(IUPAC.unambiguous_dna)
			
			#alphabet = Gapped(generic_protein)

	newAlignment = Alignment(alphabet)
	
	
	for newseq in newseqs.keys():
		print newseq, len(newseqs[newseq])
		newAlignment.add_sequence(newseq, newseqs[newseq])
	
	return newAlignment
	
	
	
	
			#sys.exit()
			
			
			
			
			
def back_translate_protein_alignment(protein_alignment_file, nucleotide_sequence_file):
	geneticcode={'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT': 'S', 'TCC': 'S','TCA': 'S','TCG': 'S', 'TAT': 'Y','TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
	try:
		palignment=Si_SeqIO.read_alignment(protein_alignment_file, quiet=True)
	except StandardError:
		DoError("Cannot open alignment file")
	
	
	try:
		nalignment=Si_SeqIO.read_seq_file(nucleotide_sequence_file, quiet=True)
	except StandardError:
		DoError("Cannot open alignment file")
	
	
	pnames={}
	for num, pseq in enumerate(palignment):
		pnames[pseq.id.split("_")[0]]=num
		
	
	#output=open("Si.aln","w")
	for nseq in nalignment:
		
		if not nseq.id.split("_")[0] in pnames:
			print nseq.id.split("_")[0], "not in palign"
			sys.exit()
		
		nsequence=str(nseq.seq).upper().replace("-","")
		
		psequence=palignment[pnames[nseq.id.split("_")[0]]].seq.upper()
		

		ppos=0
		npos=0
		count=0
		codon=[]
		newseq=[]
		inx=False
		while npos<len( nsequence)-2 and ppos<len(psequence):
			if psequence[ppos]=="-":
				newseq.append("---")
				ppos+=1
				#print "-"
				continue
			if psequence[ppos]=="X":
				ppos+=1
				inx=True
				continue
			#print nsequence[npos:npos+3], geneticcode[nsequence[npos:npos+3]], psequence[ppos]
			if geneticcode[nsequence[npos:npos+3]]==psequence[ppos]:
				if count>0 and inx:
					#print codon
					if count==1:
						codon.append("--")
						newseq.append(''.join(codon))
					else:
						codon.append("-")
						newseq.append(''.join(codon))
					count=0
					inx=False
					codon=[]
						
				#print nsequence[npos:npos+3], geneticcode[nsequence[npos:npos+3]], psequence[ppos]
				newseq.append(nsequence[npos:npos+3])
				npos+=3
				ppos+=1
			else:
				codon.append(nsequence[npos])
				npos+=1
				count+=1
				if count==3:
					#newseq.append(''.join(codon))
					count=0
					codon=[]
					inx=False
					
				
		if ppos<len(psequence):
			sys.exit()
		print ">"+nseq.id
		print ''.join(newseq)

	#output.close()
	
