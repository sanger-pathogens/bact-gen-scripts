#!/usr/bin/env python

import os, sys
from optparse import OptionParser
from Bio import SeqIO

def main():
    usage = "usage: %prog [options] <fastq file prefix>"
    parser=OptionParser(usage=usage)
    
    parser.add_option("-f", "--forwardfastq", action="store", dest="ffastq", help="Input forward fastq files", default="")
    parser.add_option("-r", "--reversefastq", action="store", dest="rfastq", help="Input reverse fastq files", default="")
    parser.add_option("-d", "--resistomedb", action="store", dest="resistome", help="Resistance db fasta file", default="")
    
    return parser.parse_args()

def DoError(errorstring):
    print "Error:", errorstring
    sys.exit()

if __name__ == "__main__":

    #Get command line arguments
    
    (options, args) = main()
    
    if options.ffastq=="":
        DoError("Forward fastq file required")
    elif not os.path.isfile(options.ffastq):
        DoError("Cannot find forward fastq file")
    
    if options.resistome=="":
        DoError("No resistance database fasta file selected")
    elif not os.path.isfile(options.resistome):
        DoError("Cannot find resistance database fasta file")

    if options.rfastq=="":
        print "No reverse fastq file selected. Assuming data is single ended"
        DoError("Single ended data not supported yet")
    elif not os.path.isfile(options.rfastq):
        DoError("Cannot find reverse fastq file")

    fastaname=options.ffastq.split(".")[0]+"_resistome.fasta"
    blastname=options.ffastq.split(".")[0]+"_resistome.blast"


    returnval=os.system(' '.join(["/nfs/users/nfs_s/sh16/scripts/resistome/fastool --to-fasta", options.ffastq, options.rfastq, ">", fastaname]))

    print "fastool returnvalue:", returnval
    
    returnval=os.system(' '.join(["makeblastdb -in", fastaname, "-dbtype nucl"]))

    print "makeblastdb returnvalue:", returnval

    returnval=os.system(' '.join(["legacy_blast.pl blastall -p blastn -d ", fastaname, "-o", blastname, "-m 8 -e 0.001 -i", options.resistome, "--path /software/pubseq/bin/ncbi_blast+/"]))

    print "blastn returnvalue:", returnval
    
    reads=set([])
    genes=set([])

    try:
        for line in open(blastname):
            line=line.strip()
            words=line.split()
            reads.add(words[1].strip("/1").strip("/2"))
            genes.add(words[0])

    except StandardError:
        DoError("Failed to read parse blast file")
        
    print len(reads), "reads match", len(genes), "genes"
