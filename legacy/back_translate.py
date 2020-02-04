#!/usr/bin/env python
import string, re, copy
import os, sys
sys.path.extend(map(os.path.abspath, ['/nfs/pathogen/sh16_scripts/modules/']))
from Si_alignment import *
#from Si_general import *
#from Si_SeqIO import *
#from Si_nexus import *
from optparse import OptionParser


if len(sys.argv)!=3 or "-h" in sys.argv or "--help" in sys.argv:
	print "Back translate an amino acid alignment to nucleotides"
	print "Output is to standard out, so needs to be redirected if you want it in a file"
	print "Usage:"
	print "/nfs/pathogen/sh16_scripts/back_translate <protein alignment file> <nucleotide mfa file>"
	
	sys.exit()

back_translate_protein_alignment(sys.argv[1], sys.argv[2])
