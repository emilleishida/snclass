import numpy as np

from fit_lc_george import lnprob2, fit_LC
from scipy import interpolate

##############################################################

class LC(object):
    """
    Light curve object.
    """

    def __init__(self, raw_data, user_choices):
  
        self.raw = raw_data				#output from function util.read_SNANA_lc
        self.user_choices = user_choices                #output from function util.read_user_input

    def check_basic(self):  
        """
        Check selection cuts which must be satisfied before any calculation.
        """

        #check if we have observed epochs in all filters
        filters_cut = all(item in self.raw.keys() for item in self.user_choices['filters'])
      
        if filters_cut:

            #check if there are at least 3 epochs in each filter
            epoch_cut = all(len(self.raw[fil]) > 2 for fil in self.user_choices['filters'])

            if epoch_cut:
                self.basic_cuts = True
            else:
                self.basic_cuts = False

        else:
            self.basic_cuts = False

    def fit_GP(self):
        "Perform Gaussian Process Fit"

        #add extra keys
        self.raw.update(self.user_choices)

        #fit light curve
        self.fitted = fit_LC(self.raw)

    def normalize(self):
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
            
            #check if we are realizations were calculated
            if int(self.user_choices['n_realizations'][0]) > 0:      
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
        self.for_matrix = {}
 
        for fil in self.user_choice['filters'][fil]:

            #create function interpolating previous results
            func = interpolate.interp1d(self.fitted['xarr_shifted'][fil], self.fitted['norm_fit'][fil])
            
            #create new horizontal axis
            xnew = np.arange(float(self.user_choices['epoch_cut'][0]), 
                             float(self.user_choices['epoch_cut'][1]), 
                             float(self.user_choices['epoch_bin'][0]))

            self.for_matrix[fil] = func(xnew)
  



            

        
         
        



        
