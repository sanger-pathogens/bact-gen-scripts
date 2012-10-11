#!/usr/bin/env python

from ghmm import *
import ghmmwrapper
import os, sys
from numpy import mean, median, max, std, bincount, argmax
#from numpy.ma import masked_equal

data=map(int,open(sys.argv[1], "rU").read().split("\n")[:-1])

#nozerodata = masked_equal(data, 0.0, copy=False)
datamean=mean(data)
datamedian=median(data)
datastd=std(data)
datamax=max(data)
#datamax=max(nozerodata)
#datamean=mean(nozerodata)
#datamedian=median(nozerodata)
#datastd=std(nozerodata)
#datamax=max(nozerodata)

counts=bincount(data)
counts=counts[1:]
datamode=argmax(counts)


print "Mean =", datamean, "Median =", datamedian, "Max =", datamax, "Stdev =", datastd, "Mode =", datamode


sigma = IntegerRange(0,max(data)+1)
F=Float()
			
	
if datamax>(datamedian*3):
	transitions = [[0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2],[0.2, 0.2, 0.2, 0.2, 0.2]]
	ezero=[0.0,datamedian/4]
	eless=[datamedian/2,datamedian/2]
	esingle=[datamedian,datamedian/2]
	emore=[datamedian*2,datamedian/2]
	ehuge=[datamedian*4,datamax]
	emissions=[ezero, eless, esingle, emore, ehuge]
	pi=[0.2, 0.2, 0.2, 0.2, 0.2]
	
	
#	transitions = [[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]]
#	ezero=[0.0,datamedian/4]
#	esingle=[datamedian,datamedian*2]
#	emore=[datamedian*2,datamedian/2]
#	ehuge=[datamedian*4,datamax]
#	emissions=[ezero, esingle, emore, ehuge]
#	pi=[0.2, 0.2, 0.2, 0.2]
else:
	transitions = [[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]]
	#transitions = [[0.97,0.1,0.1,0.1], [0.1,0.97,0.1,0.1], [0.1,0.1,0.97,0.1], [0.1,0.1,0.1,0.97]]
	ezero=[0.0,datamedian/4]
	eless=[datamedian/2,datamedian/2]
	esingle=[datamedian,datamedian/2]
	emore=[datamedian*2,datamedian/2]
	emissions=[ezero, eless, esingle, emore, ]
	pi=[0.2, 0.2, 0.2, 0.2]
	
transitions = [[0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2], [0.2, 0.2, 0.2, 0.2, 0.2],[0.2, 0.2, 0.2, 0.2, 0.2]]
ezero=[0.0,datamedian/4]
eless=[datamedian/2,datamedian/4]
esingle=[datamedian,datamedian/2]
eintermediate=[(3*datamedian)/2,datamedian/1.5]
emore=[datamedian*2,datamedian]
emissions=[ezero, eless, esingle, eintermediate, emore]
pi=[0.2, 0.2, 0.2, 0.2, 0.2]


transitions = [[1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6], [1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6], [1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6], [1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6],[1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6],[1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6]]
ezero=[0.0,datamean/4]
eless=[datamean/2,datamean/4]
esingle=[datamean,datamean/4]
eintermediate=[(3*datamean)/2,datamean/4]
emore=[datamean*2,datamean/4]
ehuge=[datamax,datamax-((datamean*2))]
emissions=[ezero, eless, esingle, eintermediate, emore, ehuge]
pi=[1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6]




#m = HMMFromMatrices(sigma, DiscreteDistribution(sigma), transitions, emissions, pi)
m = HMMFromMatrices(F,GaussianDistribution(F), transitions, emissions, pi)


#for i in range(m.cmodel.N):
#	state = ghmmwrapper.cstate_array_getRef(m.cmodel.s, i)
#	emission = ghmmwrapper.c_emission_array_getRef(m.cmodel.s.e, 0)
#	emission.fixed = 1
#	emission.setDensity(0)
#	state.e = emissions
#	state.c = ghmmwrapper.list2double_array([1.0])

for state in range(m.cmodel.N):
	s = m.cmodel.getState(state)
	for i in range(s.M):
	    emission = s.getEmission(i)
	    emission.fixed = 1



data=EmissionSequence(F, data)


initial_log = m.loglikelihood(data)
			
print "Untrained log-likelihood", initial_log
print m

#m.baumWelch(EmissionSequence(sigma, data))
m.baumWelch(data)
			
#trained_log = m.loglikelihood(EmissionSequence(sigma, data))
trained_log = m.loglikelihood(data)
print "Baum-Welch training improved the log-likelihood by ", initial_log - trained_log
sys.stdout.flush()
print m

#v = m.viterbi(EmissionSequence(sigma, data))
v = m.viterbi(data)
print v[0][:10], "..", v[0][-10:]

inblock="n"
blocks=[]
blockstart=0
secondstate=-1
firststate=-1
currstate=-1


for x, state in enumerate(v[0]):
	
#	if firststate==-1:
#		firststate=state
#		secondstate=state
#	
#	elif thirdstate!=secondstate:
#		
#		if thirdstate!=firststate:
#			blocks.append([transition,x+1,currstate])
#		
#		#if (secondstate==0 and thirdstate>1) or (secondstate==2 and (thirdstate>3 or thirdstate==0)) or 
#	
#	
#		if secondstate=1 and thirdstate==0 and firststate!=0:
#			firststate=secondstate
#			secondstate=thirdstate
#			transition=x
#		if secondstate=1 and thirdstate==0 and firststate!=0:
#			firststate=secondstate
#			secondstate=thirdstate
#			transition=x
	
	
	
	if state!=currstate and state!=2 and inblock=="n":
		blockstart=x+1
		inblock="y"
	elif state!=currstate and inblock=="y":
		blocks.append([blockstart,x+1,currstate])
		inblock="n"
		if state!=2:
			blockstart=x+1
			inblock="y"
	currstate=state
if inblock=="y":
	blocks.append([blockstart,x+1,currstate])
#print blocks

tabout=open("tabout.tab","w")
for block in blocks:
	print >> tabout, "FT   misc_feature    "+str(block[0])+".."+str(block[1])
	print >> tabout, "FT                   /colour="+str(int(block[2]+1))
	print >> tabout, 'FT                   /note="'+str(block[2])+'"'

tabout.close()	