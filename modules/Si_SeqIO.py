from Bio import AlignIO
from Bio import SeqIO
from Bio.GenBank import Scanner
#from Bio.SeqIO.InsdcIO import *
import string, os, sys
from random import *
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *

class SimonError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)





##########################################################
# Function to read an alignment whichever format it's in #
##########################################################

def tab_parser(handle, quiet=False):
	from Bio.GenBank import _FeatureConsumer
	from Bio.GenBank.utils import FeatureValueCleaner

	
	def Si_parse_tab_features(object, skip=False):
		"""Return list of tuples for the features (if present)

		Each feature is returned as a tuple (key, location, qualifiers)
		where key and location are strings (e.g. "CDS" and
		"complement(join(490883..490885,1..879))") while qualifiers
		is a list of two string tuples (feature qualifier keys and values).
		Assumes you have already read to the start of the features table.
		"""
#		if object.line.rstrip() not in object.FEATURE_START_MARKERS:
#			if object.debug : print "Didn't find any feature table"
#			return []
#	   
#		while object.line.rstrip() in object.FEATURE_START_MARKERS:
#			object.line = object.handle.readline()

		features = []
		line = object.line
		while True:
			if not line:
				break
				raise ValueError("Premature end of line during features table")
			if line[:object.HEADER_WIDTH].rstrip() in object.SEQUENCE_HEADERS:
				if object.debug : print "Found start of sequence"
				break
			line = line.rstrip()
			if line == "//":
				raise ValueError("Premature end of features table, marker '//' found")
			if line in object.FEATURE_END_MARKERS:
				if object.debug : print "Found end of features"
				line = object.handle.readline()
				break
			if line[2:object.FEATURE_QUALIFIER_INDENT].strip() == "":
				print line[2:object.FEATURE_QUALIFIER_INDENT].strip()
				raise ValueError("Expected a feature qualifier in line '%s'" % line)

			if skip:
				line = object.handle.readline()
				while line[:object.FEATURE_QUALIFIER_INDENT] == object.FEATURE_QUALIFIER_SPACER:
					line = object.handle.readline()
			else:
				#Build up a list of the lines making up this feature:
				feature_key = line[2:object.FEATURE_QUALIFIER_INDENT].strip()
				feature_lines = [line[object.FEATURE_QUALIFIER_INDENT:]]
				line = object.handle.readline()
				while line[:object.FEATURE_QUALIFIER_INDENT] == object.FEATURE_QUALIFIER_SPACER or line.rstrip() == "" : # cope with blank lines in the midst of a feature
					
					feature_lines.append(line[object.FEATURE_QUALIFIER_INDENT:].rstrip())
					line = object.handle.readline()
					if len(line)==0:
						break#EOF
				
				feature_lines.append('/seq="N"')
				sys.stdout.flush()
				features.append(object.parse_feature(feature_key, feature_lines))
		object.line = line
		
		return features
	
		
	
	def Si_feed(object, handle, consumer, do_features=True):
		"""Feed a set of data into the consumer.

		This method is intended for use with the "old" code in Bio.GenBank

		Arguments:
		handle - A handle with the information to parse.
		consumer - The consumer that should be informed of events.
		do_features - Boolean, should the features be parsed?
				      Skipping the features can be much faster.

		Return values:
		true  - Passed a record
		false - Did not find a record
		"""		
		#Should work with both EMBL and GenBank files provided the
		#equivalent Bio.GenBank._FeatureConsumer methods are called...
#		object.set_handle(handle)
		
#		if not object.find_start():
#			#Could not find (another) record
#			consumer.data=None
#			print "here"
#			return False
			           
        #We use the above class methods to parse the file into a simplified format.
        #The first line, header lines and any misc lines after the features will be
        #dealt with by GenBank / EMBL specific derived classes.

        #First line and header:
#		object._feed_first_line(consumer, object.line)
#		object._feed_header_lines(consumer, object.parse_header())

		#Features (common to both EMBL and GenBank):
		if do_features:
			object._feed_feature_table(consumer, Si_parse_tab_features(object,skip=False))
		else:
			Si_parse_tab_features(object,skip=True) # ignore the data
			
		#Footer and sequence
#		misc_lines, sequence_string = object.parse_footer()
#		object._feed_misc_lines(consumer, misc_lines)
		sequence_string="N"
		consumer.sequence(sequence_string)
#		Calls to consumer.base_number() do nothing anyway
		consumer.record_end("//")
		
		length=0
		
		for record in consumer.data.features:
			if record.location.nofuzzy_end>length:
				length=record.location.nofuzzy_end
		
		consumer.data.seq="N"*length
		
#		assert object.line == "//"

		#And we are done
		return True
	

	
	myscanner=Scanner.InsdcScanner()
	myscanner.set_handle(handle)

	myscanner.line=myscanner.handle.readline()
	myscanner.FEATURE_QUALIFIER_INDENT=21
	myscanner.FEATURE_QUALIFIER_SPACER = "FT" + " " * (myscanner.FEATURE_QUALIFIER_INDENT-2) 

	myscanner.debug=True

	
	#featuretuples=Si_parse_tab_features(myscanner)
	
	consumer = _FeatureConsumer(use_fuzziness = 1, feature_cleaner = FeatureValueCleaner())
	
	Si_feed(myscanner, handle, consumer)
	
	
	
	return consumer.data





##########################################################
# Function to read an alignment whichever format it's in #
##########################################################

def read_alignment(filename, quiet=False):
	if not quiet:
		print "Reading alignment file..."

	filetype=["phylip", "fasta", "clustal", "nexus", "emboss", "stockholm", "fasta-m10", "ig"]

	if filename.split(".")[-1].lower() in filetype:
		guesstype=filename.split(".")[-1].lower()
	elif filename.split(".")[-1].lower() in ["phy"]:
		guesstype="phylip"
	elif filename.split(".")[-1].lower() in ["fna", "dna", "aa", "aln", "fas"]:
		guesstype="fasta"
	elif filename.split(".")[-1].lower() in ["nxs", "nex", "nexus"]:
		guesstype="nexus"
	else:
		guesstype=""
	
	readok=False
	
	
	if guesstype!="":
		if not quiet:
			print "Guessing file is in "+guesstype+" format"
			print "Trying to open file "+filename+" as "+guesstype
		try:
			alignmentObject = AlignIO.read(open(filename, "rU"), guesstype)
		except StandardError:
			print "Cannot open alignment file as "+guesstype
		else:
			readok=True
			
		filetype.remove(guesstype)

	x=0

	while readok==False and x<len(filetype):
		if not quiet:
			print "Trying to open file "+filename+" as "+filetype[x]

		try:
			alignmentObject = AlignIO.read(open(filename), filetype[x])
		except StandardError:
			print "Cannot open alignment file "+filename+" as "+filetype[x]
		else:
			readok=True

		x=x+1

	if readok==False:
		raise SimonError("Failed to read alignment")
	else:
		if not quiet:
			print "Alignment read successfully"
		return alignmentObject



##########################################################
# Function to read an alignment whichever format it's in #
##########################################################

def read_seq_file(filename, quiet=False):
	if not quiet:
		print "Reading sequence file..."

	filetype=["phylip", "fasta", "clustal", "nexus", "emboss", "stockholm", "fasta-m10", "ig", "ace", "fastq"]

	if filename.split(".")[-1].lower() in filetype:
		guesstype=filename.split(".")[-1].lower()
	elif filename.split(".")[-1].lower() in ["phy", "phylip"]:
		guesstype="phylip"
	elif filename.split(".")[-1].lower() in ["fna", "dna", "aa", "aln", "mfa", "fas"]:
		guesstype="fasta"
	elif filename.split(".")[-1].lower() in ["fastq"]:
		guesstype="fastq"
	elif filename.split(".")[-1].lower() in ["nxs", "nex", "nexus"]:
		guesstype="nexus"
	else:
		guesstype=""
	
	readok=False
	
	if guesstype!="":
		if not quiet:		
			print "Guessing file is "+guesstype
			print "Trying to open file "+filename+" as "+guesstype
		
		try:
			sequenceObject = list(SeqIO.parse(open(filename, "rU"), guesstype))
		except StandardError:
			if not quiet:				
				print "Cannot open alignment file as "+guesstype
		else:
			readok=True
			
		filetype.remove(guesstype)

	x=0

	while readok==False and x<len(filetype):

		if not quiet:		
			print "Trying to open file "+filename+" as "+filetype[x]

		try:
			sequenceObject = list(SeqIO.parse(open(filename), filetype[x]))
		except StandardError:
			if not quiet:		
				print "Cannot open alignment file "+filename+" as "+filetype[x]
		else:
			readok=True

		x=x+1

	if readok==False:
		raise SimonError("1")
	else:
		if not quiet:		
			print "Sequence file read successfully"
		return sequenceObject




############################################
# Function to create an embl style ID line #
############################################

def write_embl_style_id(seqlen, handle):
	#print "ID   X00001; SV 1; linear; unassigned DNA; STD; UNC; "+str(len(sequence.strip()))+" BP.\nXX"
	#print "ID   X00001; SV 1; circular; unassigned DNA; STD; UNC; "+str(seqlen)+" BP.\nXX\n"
	handle.write("ID   X00001; SV 1; circular; unassigned DNA; STD; UNC; "+str(seqlen)+" BP.\nXX\n")
	
	return


############################################
# Function to create an embl style ID line #
############################################

def write_embl_style_header(handle):
	#print "FH   Key             Location/Qualifiers\nFH"
	handle.write("FH   Key             Location/Qualifiers\n")
	handle.write("FH\n")

	return



####################################################
# Function to create an embl style sequence format #
####################################################

def write_embl_style_sequence(sequence, handle):
	
	handle.write("XX\n")
	#sequence=sequence.strip().lower()
	#basehash={"a":0,"c":0,"g":0,"t":0,"Other":0}
	
	a_count = sequence.count('A') + sequence.count('a')
	c_count = sequence.count('C') + sequence.count('c')
	g_count = sequence.count('G') + sequence.count('g')
	t_count = sequence.count('T') + sequence.count('t')
	other = len(sequence) - (a_count + c_count + g_count + t_count)
	
	handle.write("SQ   Sequence "+str(len(sequence))+" BP; "+str(a_count)+" A; "+str(c_count)+" C; "+str(g_count)+" G; "+str(t_count)+" T; "+str(other)+" other;\n")
	
	currentposition=0
	for x in range(0,len(sequence), 60):
		handle.write("     ")
		charactersadded=0
		
		for y in range(0,60,10):
			if x+y+10<len(sequence):
				handle.write(sequence[x+y:x+y+10]+" ")
				charactersadded=charactersadded+11
				currentposition=currentposition+10
			elif x+y<len(sequence):
				#handle.write(" ")
				handle.write(sequence[x+y:]+" ")
				charactersadded=charactersadded+(1+len(sequence))-(x+y)
				currentposition=currentposition+len(sequence)-(x+y)
		handle.write(" "*(75-(charactersadded+len(str(currentposition))))+str(currentposition)+"\n")
		
	
	handle.write("//\n")

	return


###########################################################
# Function to read an annotation whichever format it's in #
###########################################################
	
def open_annotation(filename, sequence="", quiet=False):
	
	sequence=str(sequence)
	
	if not quiet:
		print "Reading annotation file..."

	filetype=["embl", "gb"]
	
	if filename.split(".")[-1].lower() in filetype:
		guesstype=filename.split(".")[-1].lower()
	else:
		guesstype=""
	
	readok=False
	
	if guesstype!="":
		if not quiet:
			print "Guessing file is "+guesstype
			print "Trying to open file "+filename+" as "+guesstype
		try:
			record = SeqIO.read(open(filename, "rU"), guesstype)
		except StandardError:
			if not quiet:
				print "Cannot open annotation file as "+guesstype
		else:
			readok=True
						
		filetype.remove(guesstype)
	
	x=0

	while readok==False and x<len(filetype):
	
		try:
			record = SeqIO.read(open(filename, "rU"), filetype[x])
		except StandardError:
			if not quiet:
				print "Cannot open annotation file as "+filetype[x]
		else:
			readok=True
			
		x=x+1
	
	
	if readok==False:
		emblfile=open(filename, "rU")
		found_id_line=False
		found_header=False
		found_sequence=False
		newembl=[]
		
		for line in emblfile:
			line=line.rstrip()
			if len(line.split())>0:
				if line.split()[0]=="ID":
					found_id_line=True
					
					IDline=line.split(";")
					IDline=map(lambda x: x.strip(),IDline)
					print IDline
					if len(IDline)>=7 and len(IDline[1].split())==2 and IDline[1].split()[0]=="SV" and IDline[2] in ["circular", "linear"] and IDline[4] in ["CON", "PAT", "EST", "GSS", "HTC", "HTG", "MGA", "WGS", "TPA", "STS", "STD"] and IDline[5] in ['PHG', 'ENV', 'FUN', 'HUM', 'INV', 'MAM', 'VRT', 'MUS', 'PLN', 'PRO', 'ROD', 'SYN', 'TGN', 'UNC', 'VRL'] and len(IDline[6].split())==2 and int(IDline[6].split()[0])==len(sequence) and IDline[6].split()[1][:2]=="BP":
						newembl.append(line)
					else:
						print len(IDline)>=7, len(IDline[1].split())==2, IDline[1].split()[0]=="SV", IDline[2] in ["circular", "linear"], IDline[4] in ["CON", "PAT", "EST", "GSS", "HTC", "HTG", "MGA", "WGS", "TPA", "STS", "STD"], IDline[5] in ['PHG', 'ENV', 'FUN', 'HUM', 'INV', 'MAM', 'VRT', 'MUS', 'PLN', 'PRO', 'ROD', 'SYN', 'TGN', 'UNC', 'VRL'], len(IDline[6].split())==2, int(IDline[6].split()[0])==len(sequence), IDline[6].split()[1][:2]=="BP"
						if not quiet:
							print "Found invalid ID line. Will be replaced."
						found_id_line=False
				elif line.split()[0]=="FH":
					found_header=True
					newembl.append(line)
				elif line.split()[0]=="SQ":
					seqlen=0
					found_sequence=True
					newembl.append(line)
				elif line.split()[0]=="FT" and line.split()[1][0]=="/":
					if len(line.split()[1].split("="))==1:
						newembl.append(line+'="True"')
					else:
						newembl.append(line)
				else:
					if found_sequence:
						try:
							seqlen=int(line.split()[-1])
						except StandardError:
							pass
					newembl.append(line)
		
		emblfile.close()
		
		print found_id_line, found_header, found_sequence
		
		if found_id_line==False or found_header==False or found_sequence==False:
			if not quiet:
				print "Last try. Trying to convert file into readable embl format"
			
			chars = string.ascii_letters + string.digits
			tmpname='tmp'+"".join(choice(chars) for x in range(randint(8, 10)))
			stringfile = open(tmpname+".embl", "w")
			
			if found_id_line==False:
				if not quiet:
					print "Adding ID line"
				if found_sequence and seqlen!=0:
					write_embl_style_id(seqlen, stringfile)
				elif len(sequence)!=0:
					write_embl_style_id(len(sequence.strip()), stringfile)
				else:
					raise SimonError("Cannot ascertain sequence length")
			
			if found_header==False:
				if not quiet:
					print "Adding header line"
				write_embl_style_header(stringfile)
			
			print >> stringfile, '\n'.join(newembl)
			
			if found_sequence==False and sequence!="":
				if not quiet:
					print "Appending sequence"
				write_embl_style_sequence(sequence, stringfile)
			elif found_sequence==False and sequence!="":
				if not quiet:
					print "Failed to find sequence in file"
				#stringfile, write_sequence(sequence, stringfile)
			
			stringfile.close()
			
			#record = SeqIO.read(open(tmpname+".embl"), "embl")
			
			#record = SeqIO.read(open(tmpname+".embl"), "embl")
			try:
				record = SeqIO.read(open(tmpname+".embl"), "embl")
				os.system("rm "+tmpname+".embl")
			except StandardError:
				os.system("rm "+tmpname+".embl")
				raise SimonError("1")
			else:
				readok=True
				
				
		else:
			return None
	elif readok==False and sequence=="":
		raise SimonError("2")

	if readok==False:
		raise SimonError("3")
	else:
		if not quiet:
			print "Annotation file successfully read."
		return record




def filter_fastq(forwardfile, reversefile=False, shuffled=False, quality_cutoff=15, length_cutoff=36, GCcutoff=False):
	
	#trims fastq reads to the first base below quality_cutoff. If the remaining read is shorter than length_cutoff, the reads (or pair) will be removed
	
	def check_pairs(readlist1, readlist2):
		if readlist1[0].split("/")[0]!=readlist2[0].split("/")[0]:
			DoError("Found dodgy read in fastq file: \n"+"\n".join(readlist1)+"\n"+"\n".join(readlist2))
	
	def check_read(readlist):
		if len(readlist)!=4 or readlist[0][0]!="@" or readlist[2][0]!="+" or len(readlist[1])!=len(readlist[3]):
			DoError("Found dodgy read in fastq file: "+"\n".join(readlist))
	
	def trim_read(readlist):
		
		for base, basequal in enumerate(readlist[3]):
			if (ord(basequal)-33)<quality_cutoff:
				readlist[3]=readlist[3][:base]
				readlist[1]=readlist[1][:base]
				
				if base<length_cutoff:
					return True
				else:
					return False
		
		return False
	
	def GC_check(readlist,cutoff):
		if not cutoff[0] in ["a", "b"]:
			DoError("Invalid GC cutoff. Must start with a or b")
		elif cutoff[0]=="a":
			greater=True
		else:
			greater=False
		try:
			value=float(cutoff[1:])
		except StandarError:
			DoError("Invalid GC cutoff. Must start with a or b followed by a float")
		
		readlen=len(readlist[1].upper().replace("N",""))
		GClen=len(readlist[1].upper().replace("N","").replace("A","").replace("T",""))
		GCpercent=(float(GClen)/readlen)*100
		
		if greater and GCpercent>value:
			return True
		elif not greater and GCpercent<value:
			return True
		else:
			return False
	
	def print_read(handle,readlist):
		for line in readlist:
			print >> handle, line
			
	
	if shuffled and reversefile:
		DoError("Cannot cope with files that are shuffled and have forward and reverse files")
	
	printtext="Filtering "
	
	if shuffled:
		printtext=printtext+"shuffled "
	
	printtext=printtext+"fastq file"
	
	if reversefile:
		printtext=printtext+"s "
	else:
		printtext=printtext+" "
	
	printtext=printtext+forwardfile
	if reversefile:
		printtext=printtext+" and "+reversefile
	
	print printtext
	sys.stdout.flush()
	
	if not os.path.isfile(forwardfile):
		DoError("Cannot find "+forwardfile)
	else:
		readfilef=open(forwardfile,"rU")
		outputf=open('.'.join(forwardfile.split(".")[:-1])+"_filtered."+forwardfile.split(".")[-1],"w")
		
	if reversefile and not os.path.isfile(reversefile):
		DoError("Cannot find "+reversefile)
	elif reversefile:
		readfiler=open(reversefile,"rU")	
		outputr=open('.'.join(reversefile.split(".")[:-1])+"_filtered."+reversefile.split(".")[-1],"w")
	
	removedcount=0
	keptcount=0
	
	for x in readfilef:
		
		toremove=False
		
		linesf=[x.strip()]
		
		for y in range(0,3):
			linesf.append(readfilef.next().strip())
			
		check_read(linesf)
		
		if shuffled:
			linesr=[]
		
			for y in range(0,4):
				linesr.append(readfilef.next().strip())
			
			check_read(linesr)
			
			check_pairs(linesf,linesr)
	
		
		elif reversefile:
			linesr=[readfiler.next().strip()]
		
			for y in range(0,3):
				linesr.append(readfiler.next().strip())
			
			check_read(linesr)
			
			check_pairs(linesf,linesr)
		
		toremove=trim_read(linesf)
		
		if not toremove and (reversefile or shuffled):
			toremove=trim_read(linesr)
		
		if not toremove and GCcutoff:
			highgc=GC_check(linesf, GCcutoff)
			if not highgc and (reversefile or shuffled):
				highgc=GC_check(linesr, GCcutoff)
			if highgc:
				#print "removed"
				toremove=True
		
		if not toremove:
			print_read(outputf, linesf)
			if reversefile:
				print_read(outputr, linesr)
			elif shuffled:
				print_read(outputf, linesr)
			keptcount+=1
		else:
			removedcount+=1
	
	printtext="Removed "+str(removedcount)+" and kept "+str(keptcount)
	if shuffled or reversefile:
		printtext=printtext+" paired reads."
	else:
		printtext=printtext+" reads."
	print printtext
	sys.stdout.flush()
	outputf.close()
	if reversefile:
		outputr.close()
	
	
	
	return
