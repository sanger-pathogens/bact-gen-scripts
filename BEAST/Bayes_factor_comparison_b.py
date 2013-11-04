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


def calculate_path_sampling(infilelist, name):

	PS_value=float("Inf")
	tf = tempfile.NamedTemporaryFile(delete=False) 
	tfname=tf.name
	print >> tf, '<beast>'
	print >> tf, '\t<pathSamplingAnalysis fileName="'+' '.join(infilelist)+'">'
	print >> tf, '\t\t<likelihoodColumn name="pathLikelihood.delta"/>'
	print >> tf, '\t\t<thetaColumn name="pathLikelihood.theta"/>'     
	print >> tf, '\t</pathSamplingAnalysis>'
	print >> tf, '</beast>'
	tf.close()
	
	try:
		psbeast = subprocess.check_output([BEAST_LOC+"beast", tfname], stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError:
		if cluster in ["farm3", "pcs5"]:
			print "\n\n!!!BEAST failed to run. As you are on farm3 or pcs5, you will need to bsub the job!!!\n"
		else:
			print "\n\n!!!BEAST failed to run!!!\n"
		os.unlink(tfname)
		sys.exit()
			
	for line in psbeast.split("\n"):
		words=line.strip().split("=")
		
		if len(words)==2 and words[0]=="log marginal likelihood (using path sampling) from pathLikelihood.delta ":
			try:
				PS_value=float(words[1])
			except ValueError:
				print >> sys.stderr, "Could not read beast path sampling results output for", name
	
	os.unlink(tfname)
	return PS_value

def calculate_stepping_stone_sampling(infilelist, name):
	
	SS_value=float("Inf")
	tf = tempfile.NamedTemporaryFile(delete=False) 
	tfname=tf.name
	print >> tf, '<beast>'
	print >> tf, '\t<steppingStoneSamplingAnalysis fileName="'+' '.join(infilelist)+'">'
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
		os.unlink(tfname)
		sys.exit()
	
	for line in ssbeast.split("\n"):
		words=line.strip().split("=")
		if len(words)==2 and words[0]=="log marginal likelihood (using stepping stone sampling) from pathLikelihood.delta ":
			try:
				SS_value=float(words[1])
			except ValueError:
				print >> sys.stderr, "Could not read beast stepping stones results output for", name
	
	os.unlink(tfname)
	return SS_value

def BF_interpretation(BF):
	if BF<0:
		print >> sys.stderr, "Error, log BF is negative!"
		comment="Err"
	elif BF==0:
		comment="No difference"
	elif BF>0 and BF<1:
		comment="Not worth more than a bare mention"
	elif BF<3:
		comment="Positive"
	elif BF<5:
		comment="Strong"
	elif BF>0:
		comment="Very strong"
		
	
	return comment


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
	if cluster in ["farm3", "pcs5"]:
		BEAST_LOC=""
	elif cluster in ["farm2", "pcs4"]:
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
	
	print "\n\n=========================\nIntra-model Bayes Factors\n========================="
	print "[Key: * better Path Sampling MLE, + better Stepping Stone Sample MLE]"
	sys.stdout.flush()
	for mlefile in mlefiles:
		if len(mlefiles[mlefile])==1:
			continue
		print "\n"+mlefile.split("/")[-1]+":"
		sys.stdout.flush()
		mlefiles[mlefile].sort()
		intraPS=[]
		intraSS=[]
		for x in xrange(len(mlefiles[mlefile])):
			name=mlefiles[mlefile][x].rstrip(options.suffix)
			PS_value=calculate_path_sampling([mlefiles[mlefile][x]], name)
			if PS_value!=float("Inf"):
				intraPS.append([name.split("/")[-1], PS_value])
			SS_value=calculate_stepping_stone_sampling([mlefiles[mlefile][x]], name)
			if SS_value!=float("Inf"):
				intraSS.append([name.split("/")[-1], SS_value])
			
		
#		print "Marginal Likelihood Estimations:"
#		print '\t'.join(["File", "PS MLE", "SSS MLE"])
#		for x in xrange(len(intraPS)):
#			print '\t'.join(map(str,[intraPS[x][0], intraPS[x][1], intraSS[x][1]]))
#		
#		print "Bayes Factors:"
		print '\t'.join(["File 1", "File 2", "PS BF", "PS comment", "SSS BF", "SSS comment"])
		for x in xrange(len(intraPS)):
				for y in xrange(x+1, len(intraPS)):
					if intraSS[x][0]!=intraPS[x][0] or intraSS[y][0]!=intraPS[y][0]:
						print "Error, expecting PS and SS samples to have same names"
						continue
					PSBF=intraPS[x][1]-intraPS[y][1]
					if PSBF==0:
						name1=intraPS[x][0]
						name2=intraPS[y][0]
					elif PSBF<0:
						PSBF=-1.0*PSBF
						name1=intraPS[x][0]
						name2=intraPS[y][0]+"*"
					else:
						name1=intraPS[x][0]+"*"
						name2=intraPS[y][0]
						
					SSBF=intraSS[x][1]-intraSS[y][1]
					
					if SSBF==0:
						name1=name1
						name2=name2
					elif SSBF<0:
						SSBF=-1.0*SSBF
						name1=name1
						name2=name2+"+"
					else:
						name1=name1+"+"
						name2=name2
						
					print '\t'.join(map(str,[name1, name2, PSBF, BF_interpretation(PSBF), SSBF, BF_interpretation(SSBF)]))
					sys.stdout.flush()
		
#		intraSS.sort()
#		
#		print "Stepping Stone Sampling:"
#		for x in xrange(len(intraSS)):
#				for y in xrange(x+1, len(intraSS)):
#					BF=intraSS[x][1]-intraSS[y][1]
#					if BF<0:
#						BF=-1.0*BF
#					print '\t'.join(map(str,[intraSS[x][0], intraSS[y][0], BF, BF_interpretation(BF)]))
#		
#			
#					sys.stdout.flush()
		
	
	
	print "\n\n=========================\nInter-model Bayes Factors\n========================="
	sys.stdout.flush()
	PS=[]
	SS=[]
	
	for mlefile in mlefiles:
		mlefiles[mlefile].sort()
		
		PS_value=calculate_path_sampling(mlefiles[mlefile], mlefile)
		if PS_value!=float("Inf"):
			PS.append([PS_value, mlefile.split("/")[-1]])
		
		SS_value=calculate_stepping_stone_sampling(mlefiles[mlefile], mlefile)
		if SS_value!=float("Inf"):
			SS.append([SS_value, mlefile.split("/")[-1]])
	
	
	
	PS.sort()
	PS.reverse()
	SS.sort()
	SS.reverse()
	foundfirst=False
	firstvalue=0
	
	print "\n\nPath Sampling Results\n=====================\n"
	
	print '\t'.join(["Prefix", "log MLE", "log BF", "Strength of Evidence (Kass & Raftery, 1995)"])
	sys.stdout.flush()
	for PSvalue in PS:
		isfirst=False
		if not foundfirst and PSvalue[0]<0:
			foundfirst=True
			firstvalue=PSvalue[0]
			isfirst=True
			
		BF=firstvalue-PSvalue[0]
		
		
		comment=BF_interpretation(BF)
			
		if isfirst:
			print '\t'.join(map(str, [PSvalue[1], PSvalue[0], "-", "-"]))
		elif foundfirst:
			print '\t'.join(map(str, [PSvalue[1], PSvalue[0], firstvalue-PSvalue[0], comment]))
		else:
			print '\t'.join(map(str, [PSvalue[1], PSvalue[0], "-", "log MLE is positive, analysis may be using improper priors"]))
		sys.stdout.flush()
	foundfirst=False
	firstvalue=0
	
	print "\n\nStepping Stone Sampling Results\n===============================\n"
	
	print '\t'.join(["Prefix", "log MLE", "log BF", "Strength of Evidence (Kass & Raftery, 1995)"])
	sys.stdout.flush()
	for SSvalue in SS:
		isfirst=False
		if not foundfirst and SSvalue[0]<0:
			foundfirst=True
			firstvalue=SSvalue[0]
			isfirst=True
		
		BF=firstvalue-SSvalue[0]
		
		comment=BF_interpretation(BF)
		
		if isfirst:
			print '\t'.join(map(str, [SSvalue[1], SSvalue[0], "-", "-"]))
		elif foundfirst:
			print '\t'.join(map(str, [SSvalue[1], SSvalue[0], firstvalue-SSvalue[0], comment]))
		else:
			print '\t'.join(map(str, [SSvalue[1], SSvalue[0], "-", "log MLE is positive, analysis may be using improper priors"]))
		sys.stdout.flush()
	print "\n\n"
	print """If you use the results of this analysis in a paper, please cite:
	
Baele G, Lemey P, Bedford T, Rambaut A, Suchard MA and Alekseyenko AV (2012) 'Improving the accuracy of demographic and molecular clock model comparison while accommodating phylogenetic uncertainty' Molecular Biology and Evolution 29(9):2157-2167

and

Baele G, Li WLS, Drummond AJ, Suchard MA and Lemey P (2013) 'Accurate model selection of relaxed molecular clocks in Bayesian phylogenetics' Molecular Biology and Evolution 30(2):239-243 """
sys.stdout.flush()	
	
	
