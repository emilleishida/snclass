#!/usr/bin/env python

from __future__ import division

import argparse
import emcee
import matplotlib.pyplot as plt
import numpy as np

import george
from george import kernels

from snclass.prepare_LC import read_user_input, read_SNANA_lc


def model(params, t):
    amp, loc, sig2 = params
    return amp * np.exp(-0.5 * (t - loc) ** 2 / sig2)


def lnprior_base(p):
    amp, loc, sig2 = p
    if not -10000 < amp < 10000:
        return -np.inf
    if not -50000 < loc < 50000:
        return -np.inf
    if not 0 < sig2 < 300.0:
        return -np.inf
    return 0.0


def lnlike_ind(p, t, y, invar):
    m = model(p[2:], t) + p[0] * t + p[1]
    return -0.5 * np.sum((y - m) ** 2 * invar)


def lnprior_ind(p):
    m, b = p[:2]
    if not -100000 < m < 100000:
        return -np.inf
    if not -10000 < b < 100000:
        return -np.inf
    return lnprior_base(p[2:])


def lnprob_ind(p, t, y, invar):
    lp = lnprior_ind(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike_ind(p, t, y, invar)


def lnlike_gp(p, t, y, yerr):
    a, tau = np.exp(p[:2])
    gp = george.GP(a * kernels.Matern32Kernel(tau))
    gp.compute(t, yerr)
    return gp.lnlikelihood(y - model(p[2:], t))


def lnprior_gp(p):
    lna, lntau = p[:2]
    if not -50000 < lna < 50000:
        return -np.inf
    if not -50000 < lntau < 50000:
        return -np.inf
    return lnprior_base(p[2:])


def lnprob_gp(p, t, y, yerr):
    lp = lnprior_gp(p)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike_gp(p, t, y, yerr)


def fit_gp(initial, data, nwalkers=32):
    ndim = len(initial)
    p0 = [np.array(initial) + 1e-8 * np.random.randn(ndim)
          for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gp, args=data)

    print("Running burn-in")
    p0, lnp, _ = sampler.run_mcmc(p0, 500)
    sampler.reset()

    print("Running second burn-in")
    p = p0[np.argmax(lnp)]
    p0 = [p + 1e-8 * np.random.randn(ndim) for i in xrange(nwalkers)]
    p0, _, _ = sampler.run_mcmc(p0, 500)
    sampler.reset()

    print("Running production")
    p0, _, _ = sampler.run_mcmc(p0, 1000)

    return sampler
def mean_gp(t, y, yerr):
    gp = george.GP( kernels.Matern32Kernel(tau))
    gp.compute(t, yerr)
    x=np.linspace(min(t), max(t), 500)
    mu, cov = gp.predict(y, x)
    std=np.sqrt(np.diag(cov))
    return x, mu, std

def main(args):

    #read_user_input
    user_input = read_user_input(args.input) 

    #read lc data
    lc_data = read_SNANA_lc(user_input)

    #set starting point
    start = [0.0, 0.0, 0.1, 0.1, 0.1]
 
    # Fit assuming GP.
    print("Fitting GP")
    t = lc_data['r'][:,0]
    y = lc_data['r'][:,1]
    yerr = lc_data['r'][:,2]

    data = (t, y, yerr)
    print 'data = ' + str(data)

    truth_gp = [0.0, 0.0] + [-1.0, 0.1, 0.4]
    sampler = fit_gp(truth_gp, data)

    # Plot the samples in data space.
    print("Making plots")
    samples = sampler.flatchain
    x = np.linspace(min(t) - 1.0, max(t) + 1.0, 500)
    xarr, mean, std=mean_gp(t, y, yerr)
    print "Mean gp is:" mean
    plt.figure()
    plt.errorbar(t, y, yerr=yerr, fmt=".k", capsize=0)
    for s in samples[np.random.randint(len(samples), size=24)]:
        gp = george.GP(np.exp(s[0]) * kernels.Matern32Kernel(np.exp(s[1])))
        gp.compute(t, yerr)
        m = gp.sample_conditional(y - model(s[2:], t), x) + model(s[2:], x)
        plt.plot(x, m, color="#4682b4", alpha=0.3)
    plt.plot(xarr, mean, 'r:', linewidth=2)
    plt.ylabel("FLUXCAL")
    plt.xlabel("MJD")
    plt.xlim(min(t) - 1.0, max(t) + 1.0)
    plt.title("George - results with Gaussian process noise model")
    plt.show()
    plt.savefig("gp-results.png", dpi=150)


if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Supernova photometric classification using KPCA.')
    parser.add_argument('-i','--input', help='Input file name',required=True)
    args = parser.parse_args()
   
    main(args)

