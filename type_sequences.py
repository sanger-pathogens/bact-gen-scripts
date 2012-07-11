#!/usr/bin/env python

#Compares mummer snp taboutput files.  NOTE: must be from the same reference sequence



##################
# Import modules #
##################

import string, re
import os, sys
import tre
from optparse import OptionParser
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
from Si_general import *
from Si_SeqIO import *
import sqlite3



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options] fasta/multifasta input files"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-s", "--sequences", action="store", dest="sequences", help="Filesname of query sequences", default="", metavar="FILE")
	parser.add_option("-d", "--database", action="store", dest="database", help="Database name", default="tmpdb", metavar="FILE")
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.sequences=='':
		DoError('No query sequence file selected!')
	elif not os.path.isfile(options.sequences):
		DoError('Cannot find file '+options.sequences+'!')
#	if options.references=='':
#		DoError('No query sequence file selected!')
#	elif not os.path.isfile(options.references):
#		DoError('Cannot find file '+options.references+'!')
		
	return


################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()

	#Do some checking of the input files
	
	check_input_validity(options, args)
	

	
	#create the sqlite3 database
	
	conn = sqlite3.connect(options.database) #creates a connection
	
	c = conn.cursor() #creates a cursor
	
	#Create table
	try:
		c.execute('''create table sequencetypes(id INTEGER PRIMARY KEY AUTOINCREMENT, sequencetype text, sequence text, old text)''')
		conn.commit()
	except StandardError:
		print "Table already exists"

	
	
	#Read the query sequences file

	try:
		queries=read_seq_file(options.sequences)
	except StandardError:
		DoError("Cannot open "+options.sequences)
	
	print '\t'.join(["Query", "Type_id", "Type_example"])
	
	for record in queries:
		data=(str(record.seq),)
		c.execute('select * from sequencetypes where sequence=?', data)
		
		hits=c.fetchall()
		
		if len(hits)==0:
			data=(record.id, str(record.seq), "FALSE")
			c.execute('insert into sequencetypes (sequencetype, sequence, old) values (?,?,?)', data)
			conn.commit()
			data=(str(record.seq),)
			c.execute('select * from sequencetypes where sequence=?', data)
			hits=c.fetchall()

		for row in hits:
			print '\t'.join(map(str, [record.id, row[0], row[1]]))
	
	print "\n"
	print '\t'.join(["Type_id", "Type_example", "Type_sequence"])
	
	c.execute('select * from sequencetypes')
	hits=c.fetchall()
	for row in hits:
		print '\t'.join(map(str,[row[0], row[1], row[2]]))

	
	