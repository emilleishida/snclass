from __future__ import division

import emcee
import os
import numpy as np

from scipy.interpolate import interp1d

import george
from george import kernels

########################################################

class LC(object):

    def __init__(self, data = None):

        self.data = data

    def check_snr(self):
        """
        Check SNR cuts.
        """
    
        cont_fil = 0
        for fil in self.data['filters']:
            cont_epoch = 0
            for item in self.data[fil][:,-1]:
                if item >= float(self.data[fil]['snr_cut'][0]):
                    cont_epoch = cont_epoch + 1

            if cont_epoch >= 3:
                cont_fil = cont_fil + 1
            
        if cont_fil == len(self.data['filters']):
            self.data['snr_cut'] = True
        else:
            self.data['snr_cut'] = False
    

    def fit_george(self):
        """
        Fits light curve using george.
        """

        for fil in self.data['filters']:

            t = self.data[fil][:,0]
            y = self.data[fil][:,1]
            yerr = self.data[fil][:,2]     

            gp = george.GP(30*np.mean(y)*kernels.ExpSquaredKernel(2*max(y)), mean = np.mean(y))
            gp.compute(t, yerr)

            p0 = gp.kernel.vector

            results = gp.optimize(t, y, yerr=yerr)

            self.data['xarr'][fil] = np.linspace(min(t), max(t), 500)
            self.data['GP_fit'][fil] = gp.predict(y, self.data['xarr'][fil])[0]

            if int(self.data['n_realizations'][0]) > 0:

                nwalkers, ndim = 24, len(gp.kernel)
                sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob2, args=[gp, y])
    
                p1 = [np.log(gp.kernel.pars) + 1e-4 * np.random.randn(ndim) for i in xrange(nwalkers)]

                if bool(self.data['screen'][0]) == True:  
                    print 'Running burn-in for ' + fil + ' band ...'

                p0, _, _ = sampler.run_mcmc(p1, nwalkers)
  
                if bool(self.data['screen'][0]) == True:
                    print 'Running production chain for ' + fil +  ' band...'

                sampler.run_mcmc(p1, nwalkers) 
                self.data['realizations'][fil] = []

                for k in xrange(int(self.data['n_realizations'][0])):
                    w = np.random.randint(sampler.chain.shape[0])
                    n = np.random.randint(nwalkers, sampler.chain.shape[1])
                    gp.kernel.pars = np.exp(sampler.chain[w, n])
                    self.data['realizations'][fil].append(gp.sample_conditional(y, data['xarr'][fil]))

    def find_max(self):
        """
        Find maximum flux in all filters and corresponding MJD.
        """

        self.data['peak_flux'] = max([max(self.data['GP_fit'][fil]) for fil in self.data['filters']])
        self.data['realizations_norm'] = {}        

        for fil in self.data['filters']:
            if self.data['peak_flux'] in self.data['GP_fit'][fil]:
                indx = list(self.data['GP_fit'][fil]).index(max_flux)
                self.data['peak_mjd'] = self.data['xarr'][fil][indx]

    def normalize(self):
        """
        Normalize flux to maximum and translate x-axis to peak MJD.
        """

        #transpose MJD and normalize
        self.data['xarr_norm'] = {}
        self.data['GP_norm'] = {}

        for fil in self.data['filters']:
            self.data['xarr_norm'][fil] = self.data['xarr'][fil] - self.data['peak_mjd']
            self.data['GP_norm'][fil] = self.data['GP_fit'][fil]/self.data['peak_flux']

            if int(self.data['n_realizations']) > 0:
                self.data['realizations_norm'][fil] = self.data['realizations'][fil]/self.data['peak_flux']

    def check_epoch(self):
        """
        Check selection cuts on epoch.
        """
   
        #earliest epoch in each filter
        min_mjd = np.array([self.data['xarr_norm'][fil][0] for fil in self.data['filters']])

        #latest epoch in each filter
        max_mjd = np.array([self.data['xarr_norm'][fil][-1] for fil in self.data['filters']])

        if min_mjd.all() <= int(self.data['epoch_cut'][0]) and max_mjd.all() >= int(self.data['epoch_cut'][1]):
            self.data['epoch_cut'] = True
        else:
            self.data['epoch_cut'] = False

    def save_GP_fit(self):
        """
        Save complete GP fit in dat files.
        """

        if not os.path.exists(data['samples_dir'][0]):
            os.makedirs(data['samples_dir'][0])

        for fil in data['filters']:
            op1 = open(self.data['samples_dir'][0] + self.data['file_root'][0] + self.data['SNID:'][0] + '_samples_' + fil + '.dat', 'w')
            for j1 in xrange(len(self.data['realizations'][fil])):
                for k1 in xrange(len(self.data['xarr'][fil])):
                    op1.write(str(self.data['realizations'][fil][j1][k1]) + '    ')
                op1.write('\n') 
            op1.close()

            op2 = open(self.data['samples_dir'][0] + self.data['file_root'][0] + self.data['SNID:'][0] + '_mean_' + fil + '.dat', 'w')
            for elem in len(self.data['xarr'][fil]):   
                op2.write(str(elem) + '\n')
            op2.close()    

            op3 = open(self.data['samples_dir'][0] + self.data['file_root'][0] + self.data['SNID:'][0] + '_xarr_' + fil + '.dat', 'w')
            for elem in  self.data['xarr'][fil]:
                op3.write(str(elem) + '\n')
            op3.close()     

    def lc_for_matrix(self):
        """
        Prepare the light curve format to compose the kernel PCA matrix.  
        """  

        xaxis_matrix = range(int(self.data['epoch_cut'][0]), int(self.data['epoch_cut'][1]))

        self.data['realizations_matrix'] = {} 
        for fil in self.data['filters']:
            f_mean = interp1d(self.data['xarr'][fil], self.data['GP_fit'][fil])
            self.data['matrix_mean'] = f_mean(xaxis_matrix)

            self.data['realizations_matrix'][fil] = []
            if int(self.data['n_realizations'][0]) > 0:
                for member in self.data['realizations'][fil]:
                    f_realizations = interp1d(self.data['xarr'][fil], member)
                    self.data['realizations_matrix'][fil].append(f_realizations(xaxis))  


