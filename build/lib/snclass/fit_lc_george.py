#!/usr/bin/env python

from __future__ import division

import argparse
import emcee
import os
import numpy as np

import george
from george import kernels

##############################################################

def lnprob2(p, gp, y):
    if np.any((-10000 > p) + (p > 10000)):
        return -np.inf

    lnprior = 0.0

    gp.kernel.pars = np.exp(p)

    return lnprior + gp.lnlikelihood(y, quiet=True)

def fit_LC(data):

    for fil in data['filters']:
        t = data[fil][:,0]
        y = data[fil][:,1]
        yerr = data[fil][:,2]     

        gp = george.GP(30*np.mean(y)*kernels.ExpSquaredKernel(2*max(y)), mean = np.mean(y))
        gp.compute(t, yerr)

        p0 = gp.kernel.vector

        results = gp.optimize(t, y, yerr=yerr)

        data['xarr'][fil] = np.linspace(min(t), max(t), 500)
        data['GP_fit'][fil] = gp.predict(y, data['xarr'][fil])[0]

        if int(data['n_realizations'][0]) > 0:

            nwalkers, ndim = 24, len(gp.kernel)
            sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob2, args=[gp, y])
    
            p1 = [np.log(gp.kernel.pars) + 1e-4 * np.random.randn(ndim) for i in xrange(nwalkers)]

            if bool(data['screen'][0]) == True:  
                print 'Running burn-in for ' + fil + ' band ...'

            p0, _, _ = sampler.run_mcmc(p1, nwalkers)
  
            if bool(data['screen'][0]) == True:
                print 'Running production chain for ' + fil +  ' band...'

            sampler.run_mcmc(p1, nwalkers) 
            data['realizations'][fil] = []

            for k in xrange(int(data['n_realizations'][0])):
                w = np.random.randint(sampler.chain.shape[0])
                n = np.random.randint(nwalkers, sampler.chain.shape[1])
                gp.kernel.pars = np.exp(sampler.chain[w, n])
                data['realizations'][fil].append(gp.sample_conditional(y, data['xarr'][fil]))

    if bool(data['save_samples'][0]) == True:

        if not os.path.exists(data['samples_dir'][0]):
            os.makedirs(data['samples_dir'][0])

        op1 = open(data['samples_dir'][0] + data['file_root'][0] + data['SNID:'][0] + '_samples.dat', 'w')
        for fil in data['filters']:
            for i1 in xrange(len(data['xarr'][fil])):   
                op1.write(fil + str(data['xarr'][fil][i1]) + '    ')
        op1.write('\n')
        for j1 in xrange(len(data['realizations'][fil])):
            for fil in data['filters']:
                for k1 in xrange(len(data['xarr'][fil])):
                    op1.write(str(data['realizations'][fil][j1][k1]) + '    ')
            op1.write('\n') 
        op1.close()

        op2 = open(data['samples_dir'][0] + data['file_root'][0] + data['SNID:'][0] + '_mean.dat', 'w')
        for fil in data['filters']:
            for i2 in xrange(len(data['xarr'][fil])):   
                op2.write(fil + str(data['xarr'][fil][i2]) + '    ')
        op2.write('\n')
        for fil in data['filters']:
            for j2 in xrange(len(data['xarr'][fil])):
                op2.write(str(data['GP_fit'][fil][j2]) + '    ')

        op2.close()        

    return data
       

def main():
  print(__doc__)

if __name__=='__main__':
  main()    

