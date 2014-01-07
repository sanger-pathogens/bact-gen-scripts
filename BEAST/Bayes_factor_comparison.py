#!/usr/bin/env python


import string
import os, sys
from optparse import OptionParser
import subprocess
import tempfile



##########################################
# Function to Get command line arguments #
##########################################


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-e", "--exclude", action="store", dest="exclude", help="File containing list of MLE files to exclude (e.g. if they have failed to converge)", default="")
	parser.add_option("-s", "--suffix", action="store", dest="suffix", help="Suffix for MLE files (used to identify correct files. [default = %default]", default=".MLE.log")
	parser.add_option("-S", "--separator", action="store", dest="separator", help="Separator used to split repeat identifier suffix from file name prefix. [default = %default]", default="_")
	parser.add_option("-n", "--nosplit", action="store_false", dest="split", help="Do not split read on separator.", default=True)
	parser.add_option("-d", "--directory", action="store", dest="directory", help="Directory containing MLE files. Default is current working directory", default="")
	#Could add more options in here so people can specify similarities etc.
	
	
	return parser.parse_args()


####################
# Get cluster name #
####################

def getclustername():
	mycluster="unknown"
	try:
		lsid_output=subprocess.check_output(["lsid"])
		
		for line in lsid_output.split("\n"):
			words=line.strip().split()
			if len(words)>0:
				if words[1]=="cluster":
					mycluster=words[4]
	
		
	except StandardError:
		return mycluster
	
	return mycluster

	
################
# Main program #
################		

if __name__ == "__main__":

	#Get command line arguments

	(options, args) = main()
	
	
	cluster=getclustername()
#	print "Running on cluster:", cluster
	if cluster in ["farm3"]:
		BEAST_LOC=""
	elif cluster in ["farm2", "pcs4", "pcs5"]:
		BEAST_LOC="/nfs/users/nfs_s/sh16/BEASTv1.7.4/bin/"
	
	if options.exclude!="":
		if not os.path.isfile(options.exclude):
			print >> sys.stderr, "Cannot find file", options.exclude
			sys.exit()
		else:
			print options.exclude
	
	if options.directory=="":
		options.directory=os.getcwd()
	elif not os.path.isdir(options.directory):
		print >> sys.stderr, "Cannot find directory", options.directory
		sys.exit()
		
	mlefiles={}
	for file in os.listdir(options.directory):
		filename=options.directory+"/"+file
		if os.path.isfile(filename) and filename[-1*len(options.suffix):]==options.suffix:
			
			tail = subprocess.check_output(["tail", "-n", "1", filename])
			words=tail.strip().split('\t')
			if len(words)!=6:
				print >> sys.stderr, "Expecting six columns in MLE file data. Found", len(words)
				print >> sys.stderr, "Skipping", file
				continue
			try:
				values=map(float,words)
			except ValueError:
				print >> sys.stderr, "Expecting all values in MLE file to be numeric."
				print >> sys.stderr, "Skipping", file
				continue
			
			if values[4]!=0:
				print >> sys.stderr, "Theta has not reached zero in", file
				print >> sys.stderr, "Skipping", file
				continue
		       	prefix=filename.rstrip(options.suffix)
			if options.split:
				if len(prefix.split(options.separator))<2:
					print >> sys.stderr, "File", filename, "cannot be split by separator:", options.separator
					print >> sys.stderr, "Using", prefix, "..."
					rootname=prefix
				else:
					rootname=options.separator.join(prefix.split(options.separator)[:-1])
			else:
				rootname=prefix

			if not rootname in mlefiles:
				mlefiles[rootname]=[]
			mlefiles[rootname].append(filename)
	
	
	PS=[]
	SS=[]
	
	mlefilesort=mlefiles.keys()
	mlefilesort.sort()
	
	for mlefile in mlefilesort:
		mlefiles[mlefile].sort()
		
		tf = tempfile.NamedTemporaryFile(delete=False) 
		tfname=tf.name
		print >> tf, '<beast>'
		print >> tf, '\t<pathSamplingAnalysis fileName="'+' '.join(mlefiles[mlefile])+'">'
		print >> tf, '\t\t<likelihoodColumn name="pathLikelihood.delta"/>'
		print >> tf, '\t\t<thetaColumn name="pathLikelihood.theta"/>'     
		print >> tf, '\t</pathSamplingAnalysis>'
		print >> tf, '</beast>'
		tf.close()
		
		print '<beast>'
		print '\t<pathSamplingAnalysis fileName="'+' '.join(mlefiles[mlefile])+'">'
		print '\t\t<likelihoodColumn name="pathLikelihood.delta"/>'
		print '\t\t<thetaColumn name="pathLikelihood.theta"/>'     
		print '\t</pathSamplingAnalysis>'
		print '</beast>'
		
		try:
			psbeast = subprocess.check_output([BEAST_LOC+"beast", tfname], stderr=subprocess.STDOUT)
		except subprocess.CalledProcessError, e:
			if cluster in ["farm3"]:
				print "\n\n!!!BEAST failed to run. As you are on farm3, you will need to bsub the job!!!\n"
			else:
				print "\n\n!!!BEAST failed to run!!!\n"
			print "Returncode =", e.returncode
			print e.output
			sys.exit()
				
		for line in psbeast.split("\n"):
			words=line.strip().split("=")
			
			if len(words)==2 and words[0]=="log marginal likelihood (using path sampling) from pathLikelihood.delta ":
				try:
					PS.append([float(words[1]), mlefile.split("/")[-1]])
				except ValueError:
					print >> sys.stderr, "Could not read beast path sampling results output for", mlefile
		
		os.unlink(tfname)
		
		tf = tempfile.NamedTemporaryFile(delete=False) 
		tfname=tf.name
		print >> tf, '<beast>'
		print >> tf, '\t<steppingStoneSamplingAnalysis fileName="'+' '.join(mlefiles[mlefile])+'">'
		print >> tf, '\t\t<likelihoodColumn name="pathLikelihood.delta"/>'
		print >> tf, '\t\t<thetaColumn name="pathLikelihood.theta"/>'     
		print >> tf, '\t</steppingStoneSamplingAnalysis>'
		print >> tf, '</beast>'
		tf.close()
		
		try:
			ssbeast = subprocess.check_output([BEAST_LOC+"beast", tfname], stderr=subprocess.STDOUT)
		except subprocess.CalledProcessError:
			if cluster in ["farm3", "pcs5"]:
				print "\n\n!!!BEAST failed to run. As you are on farm3 or pcs5, you will need to bsub the job!!!\n"
			else:
				print "\n\n!!!BEAST failed to run!!!\n"
			sys.exit()
		
		for line in ssbeast.split("\n"):
			words=line.strip().split("=")
			if len(words)==2 and words[0]=="log marginal likelihood (using stepping stone sampling) from pathLikelihood.delta ":
				try:
					SS.append([float(words[1]), mlefile.split("/")[-1]])
				except ValueError:
					print >> sys.stderr, "Could not read beast stepping stones results output for", mlefile
		
		os.unlink(tfname)
	
	
	
	PS.sort()
	PS.reverse()
	SS.sort()
	SS.reverse()
	foundfirst=False
	firstvalue=0
	
	print "\n\n=====================\nPath Sampling Results\n=====================\n"
	
	print '\t'.join(["Prefix", "log MLE", "log BF", "Strength of Evidence (Kass & Raftery, 1995)"])
	for PSvalue in PS:
		isfirst=False
		if not foundfirst and PSvalue[0]<0:
			foundfirst=True
			firstvalue=PSvalue[0]
			isfirst=True
			
		BF=firstvalue-PSvalue[0]
		
		if BF>0 and BF<1:
			comment="Not worth more than a bare mention"
		elif BF<3:
			comment="Positive"
		elif BF<5:
			comment="Strong"
		elif BF>0:
			comment="Very strong"
		else:
			print >> sys.stderr, "Error, log BF is negative!"
			
			
		if isfirst:
			print '\t'.join(map(str, [PSvalue[1], PSvalue[0], "-", "-"]))
		elif foundfirst:
			print '\t'.join(map(str, [PSvalue[1], PSvalue[0], firstvalue-PSvalue[0], comment]))
		else:
			print '\t'.join(map(str, [PSvalue[1], PSvalue[0], "-", "log MLE is positive, analysis may be using improper priors"]))
	
	foundfirst=False
	firstvalue=0
	
	print "\n\n===============================\nStepping Stone Sampling Results\n===============================\n"
	
	print '\t'.join(["Prefix", "log MLE", "log BF", "Strength of Evidence (Kass & Raftery, 1995)"])
	
	for SSvalue in SS:
		isfirst=False
		if not foundfirst and SSvalue[0]<0:
			foundfirst=True
			firstvalue=SSvalue[0]
			isfirst=True
		
		BF=firstvalue-SSvalue[0]
		
		if BF>0 and BF<1:
			comment="Not worth more than a bare mention"
		elif BF<3:
			comment="Positive"
		elif BF<5:
			comment="Strong"
		elif BF>0:
			comment="Very strong"
		else:
			print >> sys.stderr, "Error, log BF is negative!"
		
		if isfirst:
			print '\t'.join(map(str, [SSvalue[1], SSvalue[0], "-", "-"]))
		elif foundfirst:
			print '\t'.join(map(str, [SSvalue[1], SSvalue[0], firstvalue-SSvalue[0], comment]))
		else:
			print '\t'.join(map(str, [SSvalue[1], SSvalue[0], "-", "log MLE is positive, analysis may be using improper priors"]))
		
	print "\n\n"
	print """If you use the results of this analysis in a paper, please cite:
	
Baele G, Lemey P, Bedford T, Rambaut A, Suchard MA and Alekseyenko AV (2012) 'Improving the accuracy of demographic and molecular clock model comparison while accommodating phylogenetic uncertainty' Molecular Biology and Evolution 29(9):2157-2167

and

Baele G, Li WLS, Drummond AJ, Suchard MA and Lemey P (2013) 'Accurate model selection of relaxed molecular clocks in Bayesian phylogenetics' Molecular Biology and Evolution 30(2):239-243 """
	
	
	
