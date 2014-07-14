#!/usr/bin/env python

##################
# Import modules #
##################

import os, sys, math
from optparse import OptionParser, OptionGroup


##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] <list of kraken reports>"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-l", "--level", action="store", dest="level", choices=["D", "K", "P", "C", "O", "F", "G", "S", "T"], help="Taxonomic level to output. Choose from D (Domain), K (Kingdom), P (Phylum), C (Class), O (Order), F (Family), G (Genus), S (Species), T (Strain)", default="S", metavar="LEVEL", type="choice")
	parser.add_option("-c", "--cutoff", action="store", dest="cutoff", help="Cutoff minimum % match in at least one report to include in output", default=0.0, metavar="FLOAT", type="float")
	parser.add_option("-t", "--transpose", action="store_true", dest="transpose", help="Transpose output to have files in rows and matches in columns. Defautl is to have matches in rows and files in columns", default=False, metavar="BOOL")

	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.cutoff<0:
		DoError('Cutoff % match must be 0 or greater')
	elif options.cutoff>=100:
		DoError('Cutoff % match must be below 100')
		
	for arg in args:
		if not os.path.isfile(arg):
			DoError("Cannot find file "+arg)
	if len(args)==0:
		DoError("No kraken report files specified")
	
	return


Match_type={'D':'Domain', 'K':'Kingdom', 'P':'Phylum', 'C':'Class', 'O':'Order', 'F':'Family', 'G':'Genus', 'S':'Species', 'T':'Strain'}



################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	check_input_validity(options, args)
	
	values={}
	all_values={}
	UC={}
	
	files=args
	files.sort()
	
	for file in files:
		values[file]={}
		UC[file]={}
		
		species_level=float("Inf")
		fh=open(file, "rU")
		for line in fh:
			line=line.strip()
			words=line.split("\t")
			name=words[5].strip()
			PC=float(words[0])
			
			if options.level=="T" and words[3]=="-" and (len(words[5]) - len(words[5].lstrip(' ')))==species_level+2:
				values[file][name]=PC
				if not name in all_values:
					all_values[name]=[]
				all_values[name].append(PC)
			if options.level=="T" and words[3]=="S":
				species_level=(len(words[5]) - len(words[5].lstrip(' ')))
			elif words[3]==options.level:
				values[file][name]=PC
				if not name in all_values:
					all_values[name]=[]
				all_values[name].append(PC)
			elif words[3]=="U":
				UC[file]=PC
		fh.close()
	
	
	
		
	vl=[]	
	for v in all_values:
		vl.append([sum(all_values[v]), v])
	
	vl.sort()
	vl.reverse()
		
	
	
	if options.transpose:
		
		headers=[Match_type[options.level], "Unclassified"]
		
		for v in vl:
			if max(all_values[v[1]])>=options.cutoff:
				headers.append(v[1])
		print '\t'.join(map(str,headers))
		
		for file in files:
			outvalues=[file]
			if file in UC:
				outvalues.append(UC[file])
			else:
				outvalues.append(0.0)
			for v in vl:
			
				if max(all_values[v[1]])>=options.cutoff:
					
					if v[1] in values[file]:
						outvalues.append(values[file][v[1]])
					else:
						outvalues.append(0.0)
					
			print '\t'.join(map(str,outvalues))
	
	else:
	
		headers=[Match_type[options.level]]
		
		for file in files:
			headers.append(file)
		print '\t'.join(map(str,headers))
	
		outvalues=["Unclassified"]
		for file in files:
			if file in UC:
				outvalues.append(UC[file])
			else:
				outvalues.append(0.0)
		print '\t'.join(map(str,outvalues))
	
		for v in vl:
			
			if max(all_values[v[1]])>=options.cutoff:
				outvalues=[v[1]]
				for file in files:
					if file in UC:
						outvalues.append(UC[file])
					else:
						outvalues.append(0.0)
					if v[1] in values[file]:
						outvalues.append(values[file][v[1]])
					else:
						outvalues.append(0.0)
					
				print '\t'.join(map(str,outvalues))
	
	
	
			
