#!/usr/bin/env python


##################
# Import modules #
##################
import os, sys
sys.path.extend(map(os.path.abspath, ['/usr/lib/python2.4/site-packages/']))
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/lib/python2.5/site-packages/']))
from ghmm import *
import ghmmwrapper
from numpy import mean, median, max, std, bincount, argmax
from optparse import OptionParser
from numpy.ma import masked_equal



##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options]"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2012"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	parser.add_option("-a", "--average", action="store", dest="average", help="Average coverage to use for hmm. [Default= %default]", default="mean", type="choice", choices=["mean", "median", "mode"])
	parser.add_option("-d", "--data", action="store", dest="data", help="Test data file (coverage plot, one line per base). If no training data file is selected, the hmm will be trained on this data.", default="")
	parser.add_option("-t", "--training_data", action="store", dest="training", help="Training data file (coverage plot, one line per base). If no training data or genes file is selected, the hmm will be trained on your test data.", default="")
	parser.add_option("-T", "--training_genes", action="store", dest="training_genes", help="Training genes tab file (containing location of genes to use for training, such as MLST genes). If no training genes or data file is selected, the hmm will be trained on your test data.", default="")
	parser.add_option("-g", "--GC_data", action="store", dest="GC", help="GC (moving average GC plot, one line per base). Used to attempt to reduce the bias caused by GC on mapping coverage.", default="")
	parser.add_option("-m", "--min_block", action="store", dest="minblock", help="Minimum size for a block [default=%default].", default=100, type="int")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file name.", default="")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Verbose output. Will show models.", default=False)
	parser.add_option("-s", "--scaling_factor", action="store", dest="scaling_factor", help="Scale mean to be this value to reduce the efect of huge coverage. Use a value of 0 to have no scaling factor. [Default=%default]", default=200, type="int")
	parser.add_option("-f", "--flatten", action="store_true", dest="flatten", help="Flatten extremely large values", default=False)
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.data=="":
		DoError('No test data file specified')
	elif options.data and not os.path.isfile(options.data):
		DoError('Test data file '+options.data+' does not exist')
	if options.training_genes and not os.path.isfile(options.training_genes):
		DoError('Training data file '+options.training_genes+' does not exist')
	if options.GC and not os.path.isfile(options.GC):
		DoError('GC data file '+options.GC+' does not exist')
	if options.training and not os.path.isfile(options.training):
		DoError('Training data file '+options.training+' does not exist')
	if options.output=="":
		DoError('No output file name specified')
	if options.minblock<1:
		DoError('Minimum block size must be greater than 0')
	if options.scaling_factor<0:
		DoError('Scaling factor must be zero or greater')
		
	return


		
########
# Main #
########


if __name__ == "__main__":

	(options, args)=get_user_options()
	#print options, args
	check_input_validity(options, args)
	
	
	#if GC data is provided, try scaling the data values by it
	if options.GC!="":
		GCsum=[0.0]*100
		GCcount=[0.0]*100
		ratios=[0.0]*100
		data=map(int,open(options.data, "rU").read().split("\n")[:-1])
		GCdata=map(int,open(options.GC, "rU").read().split("\n")[:-1])
		
		totalsum=0.0
		totalcount=0
		for x, baseGC in enumerate(GCdata):
			GCsum[baseGC-1]+=data[x]
			GCcount[baseGC-1]+=1
			totalcount+=1
			totalsum+=data[x]
		totalmean=totalsum/totalcount
		if totalmean==0:
			print "No coverage found"
			sys.exit()
		for x, sum in enumerate(GCsum):
			if GCcount[x]>0:
				ratios[x]=(GCsum[x]/GCcount[x])/totalmean
			else:
				ratios[x]=0.0
			
			print x+1, GCsum[x], GCcount[x], ratios[x]
		
		newdata=[]
		
		output=open("newdata500.plot", "w")
		for x, datum in enumerate(data):
			if ratios[GCdata[x]-1]==0:
				newdata.append(0)
			else:
				newdata.append(int(datum/ratios[GCdata[x]-1]))
				print >> output, int(datum/ratios[GCdata[x]-1])
		output.close()
#	print newdata[:100]
	
	#If a training set is provided, calculate the mean using the training set. Otherwise use the test data. Note that the max always needds to be calculated from the test data as well
		
	if options.training=="" and options.training_genes=="":
		print "No training dataset provided\nUsing test dataset to estimate "+options.average+" coverage"
		if options.GC!="":
			data=newdata
		else:
			data=map(int,open(options.data, "rU").read().split("\n")[:-1])
		if options.average=="mean":
			dataaverage=mean(data)
		elif options.average=="median":
			dataaverage=median(data)
		elif options.average=="mode":
			counts=bincount(data)
			if len(counts)==1:
				DoError("Average coverage is zero")
			counts=counts[1:]
			dataaverage=argmax(counts)+1
			print counts
		#datastd=std(data)
		datamax=max(data)
	
	elif options.training!="":
		print "Using training dataset to estimate "+options.average+" coverage"
		data=map(int,open(options.training, "rU").read().split("\n")[:-1])
		if options.average=="mean":
			dataaverage=mean(data)
		elif options.average=="median":
			dataaverage=median(data)
		elif options.average=="mode":
			counts=bincount(data)
			if len(counts)==1:
				DoError("Average coverage is zero")
			counts=counts[1:]
			dataaverage=argmax(counts)
		#datastd=std(data)
		datamax1=max(data)
		data=map(int,open(options.data, "rU").read().split("\n")[:-1])
		datamax2=max(data)
		datamax=max([datamax1,datamax2])
	elif options.training_genes!="":
		if options.GC!="":
			data=newdata
		else:
			data=map(int,open(options.data, "rU").read().split("\n")[:-1])
		training_data=[]
		print "Using training gene locations from test dataset to estimate "+options.average+" coverage"
		for line in open(options.training_genes, "rU"):
			line=line.strip()
			words=line.split()
			if len(words)!=3 or words[0]!="FT":
				continue
			coords=words[2].replace("complement(","").replace(")","").split("..")
			if len(coords)!=2:
				continue
			if int(coords[0])<int(coords[1]):
				start=int(coords[0])
				end=int(coords[1])
			else:
				start=int(coords[1])
				end=int(coords[0])
			training_data=training_data+data[start:end+1]
			
		
		if options.average=="mean":
			dataaverage=mean(training_data)
		elif options.average=="median":
			dataaverage=median(training_data)
		elif options.average=="mode":
			counts=bincount(training_data)
			if len(counts)==1:
				DoError("Average coverage is zero")
			counts=counts[1:]
			dataaverage=argmax(counts)
		datamax1=max(training_data)
		datamin1=min(training_data)
		data=map(int,open(options.data, "rU").read().split("\n")[:-1])
		datamax2=max(data)
		datamax=max([datamax1,datamax2])
		
		print dataaverage, datamax, datamax1, datamin1, len(training_data)
#		sys.exit()
			
	
	if options.minblock>=len(data):
		DoError("Your minimum block length must be smaller than your dataset length")
	
	#If we want to remove non zero data before calculating, do this:
	#nozerodata = masked_equal(data, 0.0, copy=False)
	#nonedata = masked_equal(nozerodata, 1.0, copy=False)
#	datamax=max(nozerodata)
	#dataaverage=mean(nozerodata)
#	dataaverage=median(nozerodata)
#	datastd=std(nozerodata)
#	datamax=max(nozerodata)
	if dataaverage==0:
		DoError("Average coverage is zero")
	

#	dataaverage=1600
	#print options.average,"=", dataaverage, "Stdev =", datastd, "Max =", datamax
	print options.average,"=", dataaverage, ", max =", datamax

	
	if options.scaling_factor==0:
		scalingratio=1
	else:
	
		scalingratio=float(dataaverage)/options.scaling_factor
		data=map(lambda x: float(x)/scalingratio, data)
		
		dataaverage=dataaverage/scalingratio
		datamax=datamax/scalingratio
		print "Scaled", options.average,"=", dataaverage, ", Scaled", " max =", datamax
	
	if options.flatten:
		for x in xrange(len(data)):
			if data[x]>dataaverage*2:	
				data[x]=dataaverage*2
			datamax=dataaverage*2
		print "Flattened max =", datamax
		
	print "Setting up hmm"
	
	F=Float()
				
		
	


	#set up the model parameters for a 5 state continuous hmm
#	transitions = [[0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2],[0.2, 0.2, 0.2, 0.2, 0.2]]
#	ezero=[0.0,1.0]
#	elowintermediate=[dataaverage/2,dataaverage/4]
#	esingle=[dataaverage,dataaverage/2]
#	ehighintermediate=[(3*dataaverage)/2,dataaverage/4]
#	if dataaverage*2>datamax:
#		ehigh=[datamax,dataaverage/4]
#	else:
#		ehigh=[dataaverage*2,1.0]
#		
#	emissions=[ezero, elowintermediate, esingle, ehighintermediate, ehigh]
#	pi=[0.2, 0.2, 0.2, 0.2, 0.2]
	
	transitions = [[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]]
	ezero=[0.0,1.0]
	esingle=[dataaverage,dataaverage]
	ehighintermediate=[(3*dataaverage)/2,dataaverage/4]
	if dataaverage*2>datamax:
		ehigh=[datamax,dataaverage/4]
	else:
		ehigh=[dataaverage*2,1.0]
		
	emissions=[ezero, esingle, ehighintermediate, ehigh]
	pi=[0.25, 0.25, 0.25, 0.25]

	
	#create the hmm from the parameter matrices
	m = HMMFromMatrices(F,GaussianDistribution(F), transitions, emissions, pi)
	
	#fix the emissions so they are not optimised in the Baum Welch
	for state in range(m.cmodel.N):
		s = m.cmodel.getState(state)
		for i in range(s.M):
			emission = s.getEmission(i)
			emission.fixed = 1
	#sys.exit()
	
	#turn the data into an emission sequence
	emissiondata=EmissionSequence(F, data)
	
	if options.verbose:
		print "\nUntrained model:\n"
		#print m
		#s = m.cmodel.getState(state)
		for i in xrange(m.cmodel.N):
			print "\tstate", str(i)+":", m.getEmission(i)
			print "\t\tTransitions:",
			for j in xrange(m.cmodel.N):
				print "->"+str(j), m.getTransition(i,j),
			print
	
	#calculate the initial log likelihood (before optimisation)
	initial_log = m.loglikelihood(emissiondata)*-1	
	print "Untrained -log-likelihood", initial_log

	print "Running Baum-Welch optimisation for transitions and initial frequencies"
	sys.stdout.flush()
	#run the Baum Welch optimisation
	m.baumWelch(emissiondata)
	
	if options.verbose:
		print "\nTrained model:\n"
		#s = m.cmodel.getState(state)
		for i in xrange(m.cmodel.N):
			print "\tstate", str(i)+":", m.getEmission(i)
			print "\t\tTransitions:",
			for j in xrange(m.cmodel.N):
				print "->"+str(j), m.getTransition(i,j),
			print
	
	#calculate the initial log likelihood (after optimisation)
	trained_log = m.loglikelihood(emissiondata)*-1
	print "Trained -log-likelihood", trained_log
	print "Baum-Welch training improved the log-likelihood by ", initial_log - trained_log
	if initial_log==trained_log:
		DoError("Baum-Welch failed")
	sys.stdout.flush()
	
	print "Running Viterbi algorithm"
	sys.stdout.flush()
	#run the viterbi algorithm to calculate the best path through the data
	v = m.viterbi(emissiondata)
	

	print "Analysing Viterbi output"
	inblock="n"
	blocks=[]
	blockstart=0
	currstate=-1
	sys.stdout.flush()
#	for x, state in enumerate(v[0]):
#
#		if state!=currstate and state!=2 and inblock=="n":
#			blockstart=x+1
#			inblock="y"
#		elif state!=currstate and inblock=="y":
#			blocks.append([blockstart,x,currstate])
#			inblock="n"
#			if state!=2:
#				blockstart=x+1
#				inblock="y"
#		currstate=state
#	if inblock=="y":
#		blocks.append([blockstart,x,currstate])
#	
#	
#	#join blocks within xbp of each other, as set in the minblock option. Also rename block 5 to 4
##	for x in xrange(len(blocks)-1, 1, -1):
##		if blocks[x-1][2]==5:
##			blocks[x-1][2]=4
##		if blocks[x][2]==5:
##			blocks[x][2]=4
##		if blocks[x-1][2]==blocks[x][2] and (blocks[x][0]-blocks[x-1][1])<options.minblock:
##			blocks[x-1][1]=blocks[x][1]
##			del blocks[x]
#		
#			
#	#remove any intermediate blocks that are not directly adjacent to blocks 1 or 3 on either side
#	for x in xrange(len(blocks)-1, 1, -1):
#		if blocks[x-1][2]==1:
#			if blocks[x][2]!=blocks[x-2][2] or (blocks[x][0]-blocks[x-1][1])>1:
#				del blocks[x-1]
#		elif blocks[x-1][2]==3:
#			if blocks[x][2]!=blocks[x-2][2] or (blocks[x][0]-blocks[x-1][1])>1:
#				del blocks[x-1]
#	
#	#join adjacent blocks and colour them correctly
#	for x in xrange(len(blocks)-1, 1, -1):
#		if blocks[x-1][2] in [0,1] and blocks[x][2] in [0,1] and (blocks[x][0]-blocks[x-1][1])==1:
#			blocks[x-1][1]=blocks[x][1]
#			blocks[x-1][2]=0
#			del blocks[x]
#		elif blocks[x-1][2] in [3,4] and blocks[x][2] in [3,4] and (blocks[x][0]-blocks[x-1][1])==1:
#			blocks[x-1][1]=blocks[x][1]
#			blocks[x-1][2]=4
#			del blocks[x]
	
	print v[0][35800:36100]
	for x, state in enumerate(v[0]):

		if state!=currstate and state!=1 and inblock=="n":
			blockstart=x+1
			inblock="y"
		elif state!=currstate and inblock=="y":
			blocks.append([blockstart,x,currstate])
			inblock="n"
			if state!=1:
				blockstart=x+1
				inblock="y"
		currstate=state
	
	print blocks[:100]
	
	if inblock=="y":
		blocks.append([blockstart,x,currstate])
	
	
	print "a", blocks[:100]
			
	#remove any intermediate blocks that are not directly adjacent to blocks 1 or 2 on either side
	for x in xrange(len(blocks)-1, 1, -1):
		if blocks[x-1][2]==2:
			if blocks[x][2]!=blocks[x-2][2] or (blocks[x][0]-blocks[x-1][1])>1:
				del blocks[x-1]
	
	print "b", blocks[:100]
	
	#join adjacent blocks and colour them correctly
	for x in xrange(len(blocks)-1, 1, -1):
		if blocks[x-1][2] in [2,3] and blocks[x][2] in [2,3] and (blocks[x][0]-blocks[x-1][1])==1:
			blocks[x-1][1]=blocks[x][1]
			blocks[x-1][2]=3
			del blocks[x]
	
	print "c", blocks[:100]

	
	#sys.exit()
	data=map(int,open(options.data, "rU").read().split("\n")[:-1])
	print "Writing output files"
	tabout=open(options.output,"w")
	for block in blocks:
		#print block
		if block[2] in [0,3] and block[1]-(block[0]-1)>options.minblock:
			blockmean=mean(data[block[0]-1:block[1]])
			blockmedian=median(data[block[0]-1:block[1]])
			blockstd=std(data[block[0]-1:block[1]])
			blockmax=max(data[block[0]-1:block[1]])
			blockmin=min(data[block[0]-1:block[1]])
			blockrelativecoverage=blockmean/(dataaverage*scalingratio)
			print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
			if block[2]==0:
				print >> tabout, "FT                   /colour=1"
				print >> tabout, 'FT                   /note="Predicted region of zero coverage"'
				print >> tabout, 'FT                   /note="Mean coverage = '+str(blockmean)+'"'
				print >> tabout, 'FT                   /note="Coverage standard deviation = '+str(blockstd)+'"'
				print >> tabout, 'FT                   /note="Median coverage = '+str(blockmedian)+'"'
				print >> tabout, 'FT                   /note="Maximum coverage = '+str(blockmax)+'"'
				print >> tabout, 'FT                   /note="Minimum coverage = '+str(blockmin)+'"'
				print >> tabout, 'FT                   /note="Mean coverage relative to training data average = '+str(blockrelativecoverage)+'"'
				
			elif block[2]==3:
				print >> tabout, "FT                   /colour=2"
				print >> tabout, 'FT                   /note="Predicted region of high coverage"'
				print >> tabout, 'FT                   /note="Mean coverage = '+str(blockmean)+'"'
				print >> tabout, 'FT                   /note="Coverage standard deviation = '+str(blockstd)+'"'
				print >> tabout, 'FT                   /note="Median coverage = '+str(blockmedian)+'"'
				print >> tabout, 'FT                   /note="Maximum coverage = '+str(blockmax)+'"'
				print >> tabout, 'FT                   /note="Minimum coverage = '+str(blockmin)+'"'
				print >> tabout, 'FT                   /note="Mean coverage relative to training data average = '+str(blockrelativecoverage)+'"'
				
	
	tabout.close()	
