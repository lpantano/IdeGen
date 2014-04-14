import subprocess
import os
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from IPython.core.pylabtools import figsize
import pymc
import numpy

def getposterior(values):
	print(values)
	if numpy.min(values)!=numpy.max(values):
		mu = pymc.Uniform('mu', lower=numpy.min(values), upper=numpy.max(values))
		tau = pymc.Uniform('tau', lower=0.01, upper=10)
		y = pymc.NegativeBinomial('y',mu,tau,value=values,observed=True)
		model=pymc.MCMC([y,mu,tau])
		model.sample(iter=10000, burn=5000,verbose=0)
		return(model)
	return(False)
#print 'mu',mean(model0.trace('mu')[:])
#mu 3.62445881686



def estimatemcmc(gene,inv,values1,values2,data):
	m1 = getposterior(values1)
	m2 = getposterior(values2)
	if (m1 and m2):
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.set_title("Posterior for groups of " + gene)
		ax1.hist(m1.trace('mu')[:], color="blue", bins=30,histtype="stepfilled")
		ax1.hist(m2.trace('mu')[:], color="red", bins=30,histtype="stepfilled")
		fig.savefig(data['outcounts'] + "/" + inv +"_" + gene + ".png")
		m1up = numpy.percentile(m1.trace('mu')[:], 94.2)
		m1down = numpy.percentile(m1.trace('mu')[:], 0.25)
		m2up = numpy.percentile(m2.trace('mu')[:], 94.2)
		m2down = numpy.percentile(m2.trace('mu')[:], 0.25)
		m1mean = numpy.mean(m1.trace('mu')[:])
		m2mean = numpy.mean(m2.trace('mu')[:])
		return([m1mean,m1up,m1down,m2mean,m2up,m2down])
	else:
		return([numpy.min(values1),numpy.min(values1),numpy.min(values1),numpy.min(values2),numpy.min(values2),numpy.min(values2)])