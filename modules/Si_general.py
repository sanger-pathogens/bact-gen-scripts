import sys
from Bio.Seq import Seq



#################################
# Simple Error Printing Funtion #
#################################

def DoError(ErrorString):
	print "!!!Error:", ErrorString,"!!!"
	sys.exit()



def remove_gaps_from_sequence(sequence):
	
	mutable_seq=str(sequence)
	
	mutable_seq=mutable_seq.replace("-","")
	
	seq=Seq(mutable_seq)
	
	return seq
	

def kyte_doolittle(sequence):
	
	aaseq=Seq.translate(sequence)
	
	kytescores={'A':1.800, 'R':-4.500, 'N':-3.500, 'D':-3.500, 'C':2.500, 'Q':-3.500, 'E':-3.500, 'G':-0.400, 'H':-3.200, 'I':4.500, 'L':3.800, 'K':-3.900, 'M':1.900, 'F':2.800, 'P':-1.600, 'S':-0.800, 'T':-0.700, 'W':-0.900, 'Y':-1.300, 'V':4.200, '-':0, 'X':0, '*':0}

	kyte=[]

	for aa in aaseq:
		kyte.append(kytescores[aa])
		kyte.append(kytescores[aa])
		kyte.append(kytescores[aa])
		
	return kyte


eisscores={'A':0.620, 'R':-2.530, 'N':-0.780, 'D':-0.900, 'C':0.290, 'Q':-0.850, 'E':-0.740, 'G':0.480, 'H':-0.400, 'I':1.380, 'L':1.060, 'K':-1.500, 'M':0.640, 'F':1.190, 'P':0.120, 'S':-0.180, 'T':-0.050, 'W':0.810, 'Y':0.260, 'V':1.080, '-':0, 'X':0, '*':0}


 