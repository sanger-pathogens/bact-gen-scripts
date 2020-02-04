#!/usr/bin/env python
import os, sys, string
from pymc import DiscreteUniform, Exponential, deterministic, Poisson, Uniform, MCMC, Model
import numpy as np 

disasters_array = np.array([ 4, 5, 4, 0, 1, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6, 3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5, 2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0, 1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2, 3, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 1, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1]) 


def make_model():
	s = DiscreteUniform('s', lower=0, upper=110, doc='Switchpoint[year]')
	
	e = Exponential('e', beta=1) 
	l = Exponential('l', beta=1) 
	
	@deterministic(plot=False) 
	def r(s=s, e=e, l=l): 
		""" Concatenate Poisson means """ 
		out = np.empty(len(disasters_array)) 
		out[:s] = e 
		out[s:] = l 
		return out 
		
	D = Poisson('D', mu=r, value=disasters_array, observed=True)
#	print D.logp
	return locals()

#print D.parents
#print r.children
#print D.value
#print s.value
#print e.value
#print l.value
#print r.value
#print D.logp
#print s.logp
#print l.logp
#print e.logp

model = Model(make_model())
#print model.D.logp
M = MCMC(model)

#print M.D.logp
M.isample(iter=10000, burn=1000, thin=10)
#print M.trace('s')[:]

from pylab import hist, show
#hist(M.trace('s')[:])
#show()
from pymc.Matplot import plot
plot(M)
show()


