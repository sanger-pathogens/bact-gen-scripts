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
	parser.add_option("-c", "--cutoff", action="store", dest="cutoff", help="Cutoff minimum value in at least one report to include in output.", default=0.0, metavar="FLOAT", type="float")
	parser.add_option("-t", "--transpose", action="store_true", dest="transpose", help="Transpose output to have files in rows and matches in columns. Default is to have matches in rows and files in columns", default=False, metavar="BOOL")
	parser.add_option("-s", "--scale", action="store_true", dest="scale", help="Scale output values by the total percentage of unique matches at the chosen level for each file", default=False, metavar="BOOL")
	parser.add_option("-v", "--value", action="store", dest="value", choices=['1', '2', '3'], help="Value to output. 1 = Percentage of reads covered by the clade rooted at this taxon, 2 = Number of reads covered by the clade rooted at this taxon, 3 = Number of reads assigned directly to this taxon", default='1', metavar="VALUE", type="choice")

	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	
	if options.cutoff<0:
		DoError('Cutoff value must be 0 or greater')
	if options.value=="1":
		if options.cutoff>100:
			DoError('Cutoff value must be 100 or below when percentage of unique matches is chosen as the value')
		
	for arg in args:
		if not os.path.isfile(arg):
			DoError("Cannot find file "+arg)
	if len(args)==0:
		DoError("No kraken report files specified")
	
	return


Match_type={'D':'Domain', 'K':'Kingdom', 'P':'Phylum', 'C':'Class', 'O':'Order', 'F':'Family', 'G':'Genus', 'S':'Species', 'T':'Strain'}



##############################################################################
# Print output to stdout, but catch IOErrors that happen if you pipe to head #
##############################################################################

def print_list_to_stdout(my_list):
	try:
		sys.stdout.write('\t'.join(map(str,my_list))+"\n")
	except IOError:
		try:
			sys.stdout.close()
		except IOError:
			pass
		try:
			sys.stderr.close()
		except IOError:
			pass

    

################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	check_input_validity(options, args)
	
	values={}
	all_values={}
	totals={}
	total_assigned={}
	UC={}
	
	files=args
	files.sort()
	
	for file in files:
		values[file]={}
		UC[file]={}
		totals[file]=0.0
		total_assigned[file]=0
		
		species_level=float("Inf")
		fh=open(file, "rU")
		for line in fh:
			line=line.strip()
			words=line.split("\t")
			name=words[5].strip()
			if int(options.value)>1:
				PC=int(words[int(options.value)-1])
			else:
				PC=float(words[int(options.value)-1])
			
			if options.level=="T" and words[3]=="-" and (len(words[5]) - len(words[5].lstrip(' ')))==species_level+2:
				values[file][name]=PC
				if not name in all_values:
					all_values[name]=[]
				all_values[name].append(PC)
				totals[file]+=PC
				total_assigned[file]+=int(words[1])
			elif options.level=="T" and words[3]=="S":
				species_level=(len(words[5]) - len(words[5].lstrip(' ')))
			elif words[3]==options.level:
				values[file][name]=PC
				if not name in all_values:
					all_values[name]=[]
				all_values[name].append(PC)
				totals[file]+=PC
				total_assigned[file]+=int(words[1])
			elif words[3]=="U":
				UC[file]=PC
				totals[file]+=PC
				total_assigned[file]+=int(words[1])
		fh.close()
	
	
	
		
	vl=[]	
	for v in all_values:
		vl.append([sum(all_values[v]), v])
	
	vl.sort()
	vl.reverse()
	
	if options.transpose:
		
		headers=[Match_type[options.level], "Total", "Unclassified"]
		
		for v in vl:
			if max(all_values[v[1]])>=options.cutoff:
				headers.append(v[1])
		print_list_to_stdout(headers)
#		sys.stdout.write('\t'.join(map(str,headers))+"\n")
		
		for file in files:
			outvalues=[file]
			
			if file in total_assigned:
				outvalues.append(total_assigned[file])
			else:
				outvalues.append(0)
			
			if options.scale:
				scale_factor=(1.0/totals[file])*100
			else:
				scale_factor=1
			if file in UC:
				outvalues.append(scale_factor*UC[file])
			else:
				outvalues.append(0)
			for v in vl:
			
				if max(all_values[v[1]])>=options.cutoff:
					if v[1] in values[file]:
						outvalues.append(scale_factor*values[file][v[1]])
					else:
						outvalues.append(0)
			print_list_to_stdout(outvalues)		
#			sys.stdout.write('\t'.join(map(str,outvalues))+"\n")
	
	else:
	
		headers=[Match_type[options.level]]
		
		for file in files:
			headers.append(file)
		print_list_to_stdout(headers)
		#sys.stdout.write('\t'.join(map(str,headers))+"\n")
		
		outvalues=["Total"]
		for file in files:
			if file in total_assigned:
				outvalues.append(total_assigned[file])
			else:
				outvalues.append(0)
		print_list_to_stdout(outvalues)	
#		sys.stdout.write('\t'.join(map(str,outvalues))+"\n")
		
		outvalues=["Unclassified"]
		for file in files:
			if file in UC:
				outvalues.append(UC[file])
			else:
				outvalues.append(0)
		print_list_to_stdout(outvalues)	
#		sys.stdout.write('\t'.join(map(str,outvalues))+"\n")
	
		for v in vl:
			
			if max(all_values[v[1]])>=options.cutoff:
				outvalues=[v[1]]
				for file in files:
					if options.scale:
						scale_factor=(1.0/totals[file])*100
					else:
						scale_factor=1
					if file in UC:
						outvalues.append(scale_factor*UC[file])
					else:
						outvalues.append(0)
					if v[1] in values[file]:
						outvalues.append(scale_factor*values[file][v[1]])
					else:
						outvalues.append(0)
				print_list_to_stdout(outvalues)		
#				sys.stdout.write('\t'.join(map(str,outvalues))+"\n")
	
	
	
			
