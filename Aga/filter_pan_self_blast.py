#!/usr/bin/env python

import os, sys
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

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	group = OptionGroup(parser, "General options")
	group.add_option("-b", "--blast", action="store", dest="blast", help="BLAST output", default="", metavar="FILE")
	group.add_option("-l", "--lengths", action="store", dest="lengths", help="File containing contig length information", default="", metavar="FILE")
	group.add_option("-o", "--output", action="store", dest="output", help="Output fiel name", default="")
	group.add_option("-c", "--contigs", action="store", dest="contigs", help="Contig sequence file", default="", metavar="FILE")
	group.add_option("-s", "--minsize", action="store", dest="minsize", help="minimum size of non match within a contig to retain the contig", default=1000, type="int", metavar="INT")
	group.add_option("-L", "--minmatchlength", action="store", dest="minlength", help="minimum length of blast match", default=100, type="int", metavar="INT")
	group.add_option("-i", "--minid", action="store", dest="minid", help="Minimum identity threshold for a blast match", default=0.9, type="float", metavar="float")
	parser.add_option_group(group)


	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):


#	if options.processors>cpu_count():
#		DoError('You do not have '+str(options.processors)+' processors')
#
	if options.blast=='':
		DoError('No blast file selected')
	elif not os.path.isfile(options.blast):
		DoError('Cannot find file '+options.blast)
	if options.lengths=='':
		DoError('No lengths file selected')
	elif not os.path.isfile(options.lengths):
		DoError('Cannot find file '+options.lengths)
	if options.contigs=='':
		DoError('No contigs file selected')
	elif not os.path.isfile(options.contigs):
		DoError('Cannot find file '+options.contigs)
		
	if options.minid<0 or options.minid>100:
		DoError('Minimum identity threshold must be between 0 and 100')
	
	if options.minsize<0 or options.minid>100000:
		DoError('Minimum size threshold must be between 0 and 100000')
	
	
	return



################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	#Do some checking of the input files
	
	check_input_validity(options, args)


	lines=open(options.contigs, "rU").read().split(">")[1:]
	strain_seqs={}
	for line in lines:
		words=line.split('\n')
		name=words[0].strip().split()[0]
		seq=''.join(words[1:])
		strain_seqs[name]=seq
	
	lengths={}
	for line in open(options.lengths):
		line=line.strip()
		words=line.split()
		lengths[words[0]]=int(words[1])
		
	
	lastsubject=''
	lastquery=''
	keep=False
	first_query=True
	first_subject=True
	trim_start=0
	trim_end=0
	for line in open(options.blast):
		line=line.strip()
		words=line.split()
		
		
		if words[0]==words[1] or float(words[2])<options.minid or float(words[3])<options.minlength:
#			print words
			continue
				
		if words[0]!=lastquery:
			if not first_subject:
				unmatched_run=0
				max_run=0
				for base in match:
					if base=='0':
						unmatched_run+=1
					else:
						if unmatched_run>max_run:
							max_run=unmatched_run
						unmatched_run=0
				if unmatched_run>max_run:
					max_run=unmatched_run
					
				
				if (lengths[lastsubject]>lengths[lastquery] or (lengths[lastsubject]==lengths[lastquery] and lastsubject>lastquery)) and max_run<options.minsize:
					keep=False
					unmatched_run_score=max_run
					failed_match=''.join(map(str,match))
				
				if max_run>longest_run:
					longest_run=max_run
					
#			if keep:
#				print "Keep", lastquery, unmatched_run_score, longest_run, lengths[lastquery]#, match
#			elif not first_query:
#				print "Remove", lastquery, lastsubject, unmatched_run_score, longest_run, lengths[lastquery]#, failed_match
			if lastquery!='':
				
				print lastquery, trim_start, lengths[lastquery]-trim_end
					
			if not first_query and (not keep or lengths[lastquery]<options.minsize):
				del strain_seqs[lastquery]
			lastsubject=''
			first_subject=True
			first_query=False
			keep=True
			searching=True
			unmatched_run_score=0
			longest_run=0
			trim_start=0
			trim_end=lengths[words[0]]
			strim=1
			etrim=lengths[words[0]]
		
		if words[1]!=lastsubject:
			if strim-1>trim_start:
				trim_start=strim-1
			if etrim<trim_end:
				trim_end=etrim
			strim=1
			etrim=lengths[words[0]]
			if not first_subject:
				unmatched_run=0
				max_run=0
				for base in match:
					if base=='0':
						unmatched_run+=1
						
					else:
						if unmatched_run>max_run:
							max_run=unmatched_run
						unmatched_run=0
				if unmatched_run>max_run:
					max_run=unmatched_run
					
				if (lengths[lastsubject]>lengths[lastquery] or (lengths[lastsubject]==lengths[lastquery] and lastsubject>lastquery)) and max_run<options.minsize:
					keep=False
					unmatched_run_score=max_run
					failed_match=''.join(map(str,match))
					
				if max_run>longest_run:
					longest_run=max_run
#				if keep:
#					print "Keep", lastquery, lastsubject,  unmatched_run_score, longest_run, lengths[lastquery]#, match
#				elif not first_query:
#					print "Remove", lastquery, lastsubject, unmatched_run_score, longest_run, lengths[lastquery]#, failed_match
			first_subject=False
			match=[]
			for x in range(lengths[words[0]]):
				match.append('0')
			
		elif keep==False:
			continue
		
		
		
		fromloc=int(words[6])
		toloc=int(words[7])
		
		if fromloc>toloc:
			start=toloc
			end=fromloc
		else:
			start=fromloc
			end=toloc
		
		for x in range(start-1, end):
			match[x]=1
		
		
		sfromloc=int(words[8])
		stoloc=int(words[9])
		
		if sfromloc>stoloc:
			sstart=stoloc
			send=sfromloc
		else:
			sstart=sfromloc
			send=stoloc
		
		if send-sstart>options.minlength:
		
			if (sstart-options.minlength)<=1 and (send+options.minlength)>=lengths[words[1]]:
				s="B"
			elif (sstart-options.minlength)<=1:
				s="S"
			elif (send+options.minlength)>=lengths[words[1]]:
				s="E"
			else:
				s="M"
			if (start-options.minlength)<=strim and (end+options.minlength)>=etrim:
				
				if (lengths[words[1]]>lengths[words[0]] or (lengths[words[1]]==lengths[lastquery] and words[1]>words[0])):
					print "Remove contig", words[0]
					#print "B", s, words[0], words[1], start, end, sstart, send, lengths[words[0]], lengths[words[1]]
			elif (start-options.minlength)<=strim:
				
				if (lengths[words[1]]>lengths[words[0]] or (lengths[words[1]]==lengths[lastquery] and words[1]>words[0])) or s=="M":
					strim=end
					#print "S", s, words[0], words[1], start, end, sstart, send, lengths[words[0]], lengths[words[1]]
			elif (end+options.minlength)>=etrim:
				
				if (lengths[words[1]]>lengths[words[0]] or (lengths[words[1]]==lengths[lastquery] and words[1]>words[0])) or s=="M":
					etrim=start
					#print "E", s, words[0], words[1], start, end, sstart, send, lengths[words[0]], lengths[words[1]]
		
		
		querylen=lengths[words[0]]
		subjectlen=lengths[words[1]]
		
		lastquery=words[0]
		lastsubject=words[1]
	
	if not first_subject:
		unmatched_run=0
		max_run=0
		for base in match:
			if base=='0':
				unmatched_run+=1
			else:
				if unmatched_run>max_run:
					max_run=unmatched_run
				unmatched_run=0
		if unmatched_run>max_run:
			max_run=unmatched_run
			
		
		if (lengths[lastsubject]>lengths[lastquery] or (lengths[lastsubject]==lengths[lastquery] and lastsubject>lastquery)) and max_run<options.minsize:
			keep=False
			unmatched_run_score=max_run
			failed_match=''.join(map(str,match))
		
		if max_run>longest_run:
			longest_run=max_run
		
	if not first_query and (not keep or lengths[lastquery]<options.minsize):
		
#		if keep:
#			print "Keep", lastquery, unmatched_run_score, longest_run, lengths[lastquery]#, match
#		elif not first_query:
#			print "Remove", lastquery, lastsubject, unmatched_run_score, longest_run, lengths[lastquery]#, failed_match
#		sys.exit()

		del strain_seqs[lastquery]
	output=open(options.output, "w")
	for seq in strain_seqs:
		if len(strain_seqs[seq])>options.minsize:
			print >> output, ">"+seq
			print >> output, strain_seqs[seq]
	output.close()
