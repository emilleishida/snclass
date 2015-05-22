#!/usr/bin/env python

from __future__ import division

import numpy as np
import matplotlib.pylab as plt

from fit_lc_gptools import fit_LC
from scipy import interpolate

##############################################################

class LC(object):
    """
    Light curve object.
    """

    def __init__(self, raw_data, user_choices):
        """"
        Set parameters. 

        input: raw_data -> output from util.read_SNANA_lc
               user_choices -> output from util.read_user_input
        """
  
        self.raw = raw_data				
        self.user_choices = user_choices                


    def check_basic(self):  
        """
        Check selection cuts which must be satisfied before any calculation.

        self.basic_cuts is set to True if object passes basic selection cuts 
        (no calculations at this point, only headers)
        """

        #check if we have observed epochs in all filters
        filters_cut = all(item in self.raw.keys() for item in self.user_choices['filters'])
      
        if filters_cut:

            #apply SNR cuts
            pop = {}
            for fil in self.user_choices['filters']:
                pop[fil] = []
                for line in self.raw[fil]:
                    if float(line[-1]) >= float(self.user_choices['quality_cut'][0]):
                        pop[fil].append(line)
  
            #check if there are at least 3 epochs in each filter
            epoch_cut = all(len(pop[fil]) > 2 for fil in self.user_choices['filters'])

            if epoch_cut:
                self.basic_cuts = True
            else:
                self.basic_cuts = False

        else:
            self.basic_cuts = False

    def fit_GP(self, **kwargs):
        """
        Perform Gaussian Process Fit.

        self.fitted -> dictionary of fitted parameters
                       
        """

        #add extra keys
        self.raw.update(self.user_choices)

        #fit light curve
        self.fitted = fit_LC(self.raw, **kwargs)

    def normalize(self, samples=False):
        "Normalize according to maximum flux in all filters."                      

        #determine maximum flux 
        self.fitted['max_flux'] = max([max(self.fitted['GP_fit'][item]) 
                            for item in self.user_choices['filters']])

        #normalize
        self.fitted['norm_fit'] = {}
        self.fitted['norm_realizations'] = {}
        for fil in self.user_choices['filters']:
            self.fitted['norm_fit'][fil] = [elem/self.fitted['max_flux'] 
                                            for elem in self.fitted['GP_fit'][fil]]
            
            #check if  realizations were calculated
            if samples == True and int(self.user_choices['n_samples'][0]) > 0:     
                self.fitted['norm_realizations'][fil] = [elem/self.fitted['max_flux'] 
                                                         for elem in self.fitted['realizations'][fil]]        

    def mjd_shift(self):
        "Determine day of maximum and shift all epochs."
        
        #determine day of maximum
        self.fitted['peak_mjd_fil'] = [fil for fil in self.user_choices['filters'] 
                                                      if 1.0 in self.fitted['norm_fit'][fil]][0]
        pkmjd_indx = self.fitted['norm_fit'][self.fitted['peak_mjd_fil']].index(1.0)         
        self.fitted['peak_mjd'] = self.fitted['xarr'][self.fitted['peak_mjd_fil']][pkmjd_indx]

        #shift light curve
        self.fitted['xarr_shifted'] = {}
        for fil in self.user_choices['filters']:
            self.fitted['xarr_shifted'][fil] = [elem - self.fitted['peak_mjd'] for elem in self.fitted['xarr'][fil]]

    def check_epoch(self):
        "Check if all filters satisfy epoch coverage requirements."

        #store epoch flags
        epoch_flags = []

        for fil in self.user_choices['filters']:
            if min(self.fitted['xarr_shifted'][fil]) <= int(self.user_choices['epoch_cut'][0]) and max(self.fitted['xarr_shifted'][fil]) >= int(self.user_choices['epoch_cut'][1]):
                epoch_flags.append(True)
            else:
                epoch_flags.append(False)

        self.epoch_cuts = all(test == True for test in epoch_flags)

    def build_steps(self):
        "Build lines for the initial data matrix"

        #create dictionary to store results
        self.flux_for_matrix = {}
 
        for fil in self.user_choices['filters']:

            #create function interpolating previous results
            self.func = interpolate.interp1d(self.fitted['xarr_shifted'][fil], self.fitted['norm_fit'][fil])
            
            #create new horizontal axis
            self.xnew = np.arange(float(self.user_choices['epoch_cut'][0]), 
                             float(self.user_choices['epoch_cut'][1]), 
                             float(self.user_choices['epoch_bin'][0]))

            self.flux_for_matrix[fil] = self.func(self.xnew)
  


    def plot_fitted(self, samples=False, nsamples=0, file_out=None):
        """
        Plotted light curve as it enters the data matrix.

        input:  samples ->   bool, optional
                             Rather or not to plot realizations 
                             from the posterior.
                nsamples ->  int, optional
                             number of samples to draw from the posterior
                             only effective if samples = True
                file_out >   bool, optional
                             File name where to store the final plot.
                             If None shows the plot in the screen.
                             Default is None.

        output: if file_out is str -> plot wrote to file
        """

        #set the number of samples variable according to input
        if samples == False:
            nsamples = 0

        f = plt.figure()
        for i in xrange(len(self.user_choices['filters'])): 
            fil =  self.user_choices['filters'][i] 
            plt.subplot(2, len(self.user_choices['filters'])/2 + 
                           len(self.user_choices['filters'])%2, i + 1)
            plt.title('filter = ' + fil)
            plt.plot(self.fitted['xarr_shifted'][fil], 
                     self.fitted['norm_fit'][fil], color='red')
            
            #plot samples
            if samples == True:
                for s in self.fitted['realizations'][fil]:
                    plt.plot(self.fitted['xarr_shifted'][fil], s/self.fitted['max_flux'], 
                             color="#4682b4", alpha=0.3)
            plt.errorbar(self.raw[fil][:,0] - self.fitted['peak_mjd'], 
                        self.raw[fil][:,1]/self.fitted['max_flux'],
                        yerr=self.raw[fil][:,2]/self.fitted['max_flux'], 
                        color='blue', fmt='o')
            plt.xlabel('days since maximum', fontsize=15)
            plt.ylabel('normalized flux', fontsize=15)
            plt.xlim(float(self.user_choices['epoch_cut'][0]), 
                     float(self.user_choices['epoch_cut'][1]))
        f.tight_layout()
            
        if isinstance(file_out, str):
            plt.savefig(file_out)
            f.close()
        else:             
            plt.show()    


            

        
         
        



        
