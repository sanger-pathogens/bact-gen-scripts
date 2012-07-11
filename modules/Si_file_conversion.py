
########################################################################################################################################
# A function to split really old solexa paired-end files which had both reads in one file concatenated together into one long sequence #
########################################################################################################################################

def split_fastq_with_concatenated_read_pairs(input_handle, forward_handle, reverse_handle):
	
	#for each read in the input file
	for x in input_handle:
		#line x line should be a read header line
		
		if x[0]!="@":
			print "Error, did not find @ symbol at start of what should be a header line. Are there clank lines in your fastq file?"
			sys.exit()
		
		lines=[x.strip()]
		#read the next three lines to get the sequence and quality for the read
		for y in range(0,3):
			lines.append(input_handle.next().strip())
		
		#print the read header to the two output files with appropriate suffixes for direction
		print >> forward_handle, lines[0]+"/1"
		print >> reverse_handle, lines[0]+"/2"
		
		#print the first half of the sequence line to the forward output file
		print >> forward_handle, lines[1][:len(lines[1])/2]
		#print the second half of the sequence line to the forward output file
		print >> reverse_handle, lines[1][len(lines[1])/2:]
		
		#print the line between the sequence and quality to both files
		print >> forward_handle, lines[2]
		print >> reverse_handle, lines[2]
		
		#print the first half of the quality line to the forward output fil
		print >> forward_handle, lines[3][:len(lines[3])/2]
		#print the second half of the quality line to the forward output file
		print >> reverse_handle, lines[3][len(lines[3])/2:]
		
	return