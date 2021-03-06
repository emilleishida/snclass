"""
Created by Emille Ishida in May, 2015.

Class to deal with operations on light curves.
"""

from __future__ import division

import numpy as np
import os
import matplotlib.pylab as plt
from scipy import interpolate

from snclass.fit_lc_gptools import fit_lc
from snclass.util import read_fitted, read_snana_lc
from snclass.functions import screen

##############################################################


class LC(object):

    """
    Light curve object.

    Methods:
        - check_basic: Check selection cuts from raw curve.
        - fit_gp: Perform Gaussian Process Fit.
        - load_fit_GP: Load previously calculated GP fit.
        - normalize: Normalize according to maximum flux in all filters.
        - mjd_shift: Determine day of maximum and shift all epochs.
        - check_epoch: Check if all filters satisfy epoch requirements.
        - build_steps: Build lines for the initial data matrix.
        - plot_raw: Plot raw data.
        - plot_fitted: Plotted light curve as it enters the data matrix.

    Attributes:
        - raw, dict: raw data
        - user_choices, dict: user input choices
        - basic_cuts, bool: basic cuts flag
        - fitted, dict: results from GP fit
        - epoch_cuts, bool: epoch cuts flag
        - flux_for_matrix, dict: results for data matrix lines
        - xnew: data matrix cadence
        - samples_for_matrix: normalized random GP realizations
    """

    def __init__(self, raw_data, user_choices):
        """
        Set parameters.

        input: raw_data -> output from util.read_snana_lc
               user_choices -> output from util.read_user_input
        """
        self.raw = raw_data
        self.user_choices = user_choices
        self.basic_cuts = None
        self.fitted = None
        self.epoch_cuts = None
        self.flux_for_matrix = {}
        self.func = None
        self.xnew = None
        self.samples_for_matrix = {}

    def check_basic(self):
        """
        Check selection cuts which must be satisfied before any calculation.

        self.basic_cuts is set to True if object passes basic selection cuts
        (no calculations at this point, only headers)
        """
        # check if we have observed epochs in all filters
        filters_cut = all(item in self.raw.keys()
                          for item in self.user_choices['filters'])

        if filters_cut:

            # apply SNR cuts
            pop = {}
            for fil in self.user_choices['filters']:
                pop[fil] = []
                for line in self.raw[fil]:
                    quality = float(self.user_choices['quality_cut'][0])
                    if float(line[-1]) >= quality and self.user_choices['photon_flag'][0] == 'FLUXCAL':
                        pop[fil].append(line)
                    elif float(line[-1]) >= quality and self.user_choices['photon_flag'][0] == 'mag' and float(line[1]) < 50.0 and float(line[2]) < 50.0:
                        pop[fil].append(line)

            # check if there are at least 3 epochs in each filter
            epoch_cut = all(len(pop[fil]) >= int(self.user_choices['nepoch_cut'][0])
                            for fil in self.user_choices['filters'])

            if epoch_cut:
                self.basic_cuts = True
            else:
                self.basic_cuts = False

        else:
            self.basic_cuts = False

    def fit_GP(self, mean=True, samples=False, screen=False, do_mcmc=True,
               save_mean=True, save_samples=False):
        """
        Perform Gaussian Process Fit.

        input: mean -> bool, optional
                       if True, calculate mean GP fit
                       Default is True

               samples -> bool, optional
                          if True, calculate samples from the final GP
                          Default is False

              screen -> bool, optional
                        if True, print calculation steps into screen
                        Default is False

              do_mcmc -> bool, optional
                         if True, optimize kernel parameters using mcmc
                         Default is True

              save_mean -> bool, optional
                           if True save mean GP fit to file
                           Default is True

              save_samples -> bool, optional
                              if True save GP draws to file
                              Default is False
        """
        # add extra keys
        self.raw.update(self.user_choices)

        # fit light curve
        if self.user_choices['epoch_cut'][0] == '-999' and self.user_choices['photon_flag'][0] == 'FLUXCAL':
            p1 = [int(self.user_choices['epoch_predict'][0]), 
                  int(self.user_choices['epoch_predict'][1])] 
        else:
            p1 = None

        self.fitted = fit_lc(self.raw, mean=mean, samples=samples,
                             screen=screen, do_mcmc=do_mcmc,
                             save_mean=save_mean, save_samples=save_samples,
                             predict=p1)

    def load_fit_GP(self, mean_file):
        """
        Load previously calculated GP fit.

        input: mean_file, str
               file with previously fitted GP result
        """
        # add extra keys
        self.raw.update(self.user_choices)

        # load
        self.fitted = read_fitted(self.raw, mean_file)

    def normalize(self, samples=False, ref_filter=None):
        """
        Normalize GP fit and samples.

        input: samples, bool - optional
               if True, normalize realizations from the GP as well
               Default is False

               ref_filter, str - optional
               Reference filter for normalization
               Default is None
        """
        # determine maximum flux
        if ref_filter == None:
            self.fitted['max_flux'] = max([max(self.fitted['GP_fit'][item])
                                           for item in
                                           self.user_choices['filters']])
        else:
            self.fitted['max_flux'] = max(self.fitted['GP_fit'][ref_filter])

        # normalize
        self.fitted['norm_fit'] = {}
        self.fitted['norm_realizations'] = {}
        for fil in self.user_choices['filters']:
            max_f = self.fitted['max_flux']
            gp_fit = self.fitted['GP_fit'][fil]
            self.fitted['norm_fit'][fil] = np.array([elem / max_f
                                                     for elem in gp_fit])

            # check if  realizations were calculated
            if samples and int(self.user_choices['n_samples'][0]) > 0:
                max_f = self.fitted['max_flux']
                gp_fitted = self.fitted['realizations'][fil]
                self.fitted['norm_realizations'][fil] = np.array([np.array(elem) / max_f
                                                                  for elem in
                                                                  gp_fitted])

    def mjd_shift(self):
        """Determine day of maximum and shift all epochs."""
        # determine day of maximum
        self.fitted['peak_mjd_fil'] = [fil for fil in
                                       self.user_choices['filters'] if 1.0 in
                                       self.fitted['norm_fit'][fil]][0]

        max_fil = self.fitted['peak_mjd_fil']
        pkmjd_indx = list(self.fitted['norm_fit'][max_fil]).index(1.0)
        self.fitted['peak_mjd'] = self.fitted['xarr'][max_fil][pkmjd_indx]

        # shift light curve
        self.fitted['xarr_shifted'] = {}
        for fil in self.user_choices['filters']:
            pkmjd = self.fitted['peak_mjd']
            xlist = self.fitted['xarr'][fil]
            self.fitted['xarr_shifted'][fil] = np.array([elem - pkmjd
                                                         for elem in xlist])

    def check_epoch(self):
        """Check if all filters satisfy epoch coverage requirements."""
        # store epoch flags
        epoch_flags = []

        for fil in self.user_choices['filters']:
            if self.user_choices['epoch_cut'][0] == '-999' and self.user_choices['photon_flag'][0] == 'FLUXCAL':
                epoch_flags.append(True)

            elif (min(self.fitted['xarr_shifted'][fil]) <=
               int(self.user_choices['epoch_cut'][0])) and \
               (max(self.fitted['xarr_shifted'][fil]) >=
               int(self.user_choices['epoch_cut'][1])):
                epoch_flags.append(True)

            else:
                epoch_flags.append(False)

        self.epoch_cuts = all(test for test in epoch_flags)

    def build_steps(self, samples=False):
        """
        Build lines for the initial data matrix.

        input: samples, bool, optional
               if True built steps for normalized realizations
               default is False
        """
        self.mean_for_matrix = []
        for fil in self.user_choices['filters']:


            xaxis = self.fitted['xarr_shifted'][fil]
            yaxis = self.fitted['norm_fit'][fil]
            # create function interpolating previous results
            func = interpolate.interp1d(xaxis, yaxis)

            # create new horizontal axis
            if self.user_choices['epoch_cut'][0] == '-999' and self.user_choices['photon_flag'][0] == 'FLUXCAL':
                xmin = float(self.user_choices['epoch_predict'][0])
                xmax = float(self.user_choices['epoch_predict'][1])
            elif  self.user_choices['epoch_cut'][0] == '-999' and self.user_choices['photon_flag'][0] == 'mag':
                xmin = min(xaxis)
                xmax = max(xaxis)
            else:
                xmin = float(self.user_choices['epoch_cut'][0])
                xmax = float(self.user_choices['epoch_cut'][1])
         
            xstep = float(self.user_choices['epoch_bin'][0])
            xnew = np.arange(xmin, xmax, xstep)

            self.flux_for_matrix[fil] = func(xnew)

            for item in self.flux_for_matrix[fil]:
                self.mean_for_matrix.append(item)

        if samples:
            self.samples_for_matrix = []
            fini = self.user_choices['filters'][0]
            for j in xrange(len(self.fitted['norm_realizations'][fini])):   
                line = [] 
                for fil in self.user_choices['filters']:
                    xaxis2 = self.fitted['xarr_shifted'][fil]
                    item = self.fitted['norm_realizations'][fil][j]
                    # create function interpolating previous results
                    func_samp = interpolate.interp1d(xaxis2, item)
                    # calculate sample grid in epochs
                    new_grid = func_samp(xnew)

                    # store
                    for element in new_grid:
                        line.append(element)
                self.samples_for_matrix.append(line)

    def plot_raw(self, file_out=None):
        """
        Plot raw data.

        input: file_out ->   bool, optional
                             File name where to store the final plot.
                             If None shows the plot in the screen.
                             Default is None.

        output: if file_out is str -> plot wrote to file
        """
        # set the number of samples variable according to input

        if self.user_choices['photon_flag'][0] == 'FLUXCAL':
            xmin = float(self.user_choices['epoch_predict'][0])
            xmax = float(self.user_choices['epoch_predict'][1])
            sign = 1.0

        elif self.user_choices['photon_flag'][0] == 'MAG':
            sign = -1.0

        my_fig = plt.figure(figsize=(10,10))
        plt.suptitle('SNID ' + self.raw['SNID:'][0])
        for i in xrange(len(self.user_choices['filters'])):

            fil = self.user_choices['filters'][i]

            plt.subplot(2, len(self.user_choices['filters']) / 2 +
                        len(self.user_choices['filters']) % 2, i + 1)
            my_axis = plt.gca()
            plt.title('filter = ' + fil)

            if sign == 1.0 and len(self.raw[fil][:,0]) > 0:
                plt.errorbar(self.raw[fil][:, 0],
                             sign * self.raw[fil][:, 1],
                             yerr=self.raw[fil][:, 2],
                             color='blue', fmt='o')
            elif len(self.raw[fil][:,0]) > 0:
                plt.errorbar(self.raw[fil][:, 0],
                             sign * self.raw[fil][:, 1],
                             yerr=self.raw[fil][:, 2], color='blue', fmt='o')

            plt.xlabel('days since maximum', fontsize=15)

            if sign == -1.0 and len(self.raw[fil][:,0]) > 0:
                ax = plt.gca()
                plt.ylim(min(sign * self.raw[fil][:, 1]) - 1.5*max(self.raw[fil][:, 2]),max(sign *self.raw[fil][:, 1]) + 1.5*max(self.raw[fil][:, 2]))  
                ax.invert_yaxis()

        my_fig.tight_layout()

        if isinstance(file_out, str):
            plt.savefig(self.user_choices['path_output_plot'][0] + file_out)
            plt.close()
        else:
            plt.show()
                

    def plot_fitted(self, file_out=None):
        """
        Plotted light curve as it enters the data matrix.

        input: file_out >   bool, optional
                             File name where to store the final plot.
                             If None shows the plot in the screen.
                             Default is None.

        output: if file_out is str -> plot wrote to file
        """
        # set the number of samples variable according to input
        samples = bool(int(self.user_choices['n_samples'][0]))

        if self.user_choices['photon_flag'][0] == 'FLUXCAL':
            xmin = float(self.user_choices['epoch_predict'][0])
            xmax = float(self.user_choices['epoch_predict'][1])
            sign = 1.0

        elif self.user_choices['photon_flag'][0] == 'mag':
            sign = -1.0

        my_fig = plt.figure(figsize=(10,10))
        plt.suptitle('SNID ' + self.raw['SNID:'][0])
        for i in xrange(len(self.user_choices['filters'])):

            fil = self.user_choices['filters'][i]

            plt.subplot(2, len(self.user_choices['filters']) / 2 +
                        len(self.user_choices['filters']) % 2, i + 1)
            my_axis = plt.gca()
            plt.title('filter = ' + fil)

            if sign == 1.0:
                func = interpolate.interp1d(self.fitted['xarr_shifted'][fil],
                                        self.fitted['norm_fit'][fil])
                plt.plot(self.fitted['xarr_shifted'][fil],
                         self.fitted['norm_fit'][fil], color='red')
                plt.ylabel('normalized flux', fontsize=15)
                plt.errorbar(self.raw[fil][:, 0] - self.fitted['peak_mjd'],
                             sign * self.raw[fil][:, 1] / self.fitted['max_flux'],
                             yerr=self.raw[fil][:, 2] / self.fitted['max_flux'],
                             color='blue', fmt='o')
            else:
                func = interpolate.interp1d(self.fitted['xarr_shifted'][fil],
                                        self.fitted['GP_fit'][fil])
                plt.plot(self.fitted['xarr_shifted'][fil],
                         sign * np.array(self.fitted['GP_fit'][fil]), color='red')
                plt.ylabel('magnitude', fontsize=15)
                plt.errorbar(self.raw[fil][:, 0] - self.fitted['peak_mjd'],
                             sign * self.raw[fil][:, 1],
                             yerr=self.raw[fil][:, 2], color='blue', fmt='o')

            # plot samples
            if samples and sign == 1.0:
                for curve in self.fitted['realizations'][fil]:
                    plt.plot(self.fitted['xarr_shifted'][fil],
                             sign * np.array(curve) / self.fitted['max_flux'],
                             color='gray', alpha=0.3)
            elif samples and sign == -1.0:
                for curve in self.fitted['realizations'][fil]:
                    plt.plot(self.fitted['xarr_shifted'][fil],
                             sign * np.array(curve), color='gray', alpha=0.3)
                
            plt.xlabel('days since maximum', fontsize=15)
            plt.xlim(min(self.fitted['xarr_shifted'][fil]) - 1.0,
                     max(self.fitted['xarr_shifted'][fil]) + 1.0)
            if self.user_choices['photon_flag'][0] == 'FLUXCAL' and self.user_choices['epoch_cut'][0] != '-999':
                plt.vlines(xmin, my_axis.get_ylim()[0], func(xmin), color='black',
                           linestyles='dashed')
                plt.vlines(xmax, my_axis.get_ylim()[0], func(xmax), color='black',
                           linestyles='dashed')

            if sign == -1.0:
                ax = plt.gca()
                plt.ylim(min(sign * self.raw[fil][:, 1]) - 1.5*max(self.raw[fil][:, 2]),max(sign *self.raw[fil][:, 1]) + 1.5*max(self.raw[fil][:, 2]))  
                ax.invert_yaxis()

        my_fig.tight_layout()

        if isinstance(file_out, str):
            plt.savefig(self.user_choices['path_output_plot'][0] + file_out)
            plt.close()
        else:
            plt.show()


def fit_objs(user_choices, plot=False, calc_mean=True, calc_samp=False,
             save_samp=False):
    """
    Perform a GP fit in a set of objects.

    input: user_choices
           output from read_user_input

           plot - bool, optional
           rather or not to generate an output png file
           default is False

           calc_mean - bool, optional
           rather or not to calculate the main GP fit
           default is True

           calc_samp - bool, optional
           rather or not to calulate realizations of the final fit
           default is False

           save_samp - bool, optional
           rather or not to save realizations of the final fit
           default is False
    """
    if not os.path.exists(user_choices['samples_dir'][0]):
        os.makedirs(user_choices['samples_dir'][0])

    # read list of SN in sample
    f_open = open(user_choices['snlist'][0], 'r')
    lin = f_open.readlines()
    f_open.close()

    snlist = [elem.split()[0] for elem in lin]

    fit_method = bool(int(user_choices['do_mcmc'][0]))

    for supernova in snlist:

        # update object
        user_choices['path_to_lc'] = [supernova]

        # read light curve raw data
        raw = read_snana_lc(user_choices)

        if not os.path.isfile(user_choices['samples_dir'][0] + \
                              user_choices['file_root'][0] + \
                              raw['SNID:'][0] + '_' + user_choices['photon_flag'][0] +  '_mean.dat'):

            # initiate light curve object
            my_lc = LC(raw, user_choices)

            screen('Fitting SN' + raw['SNID:'][0], user_choices)

            # perform basic check
            my_lc.check_basic()

            # check if satisfy minimum cut
            if my_lc.basic_cuts:
                screen('... Passed basic cuts', user_choices)

                # fit
                my_lc.fit_GP(mean=calc_mean, samples=calc_samp, 
                             do_mcmc=fit_method, save_samples=save_samp, 
                             screen=bool(int(user_choices['screen'][0])))

                if plot:
                    my_lc.normalize()
                    my_lc.mjd_shift()
                    my_lc.plot_fitted(file_out=user_choices['path_output_plot'][0] + 'gp-SN' + raw['SNID:'][0] + '_' + user_choices['photon_flag'][0] + '.png')

                print '\n'

            else:
                screen('Failed to pass basic cuts!\n', user_choices)

        else:
            screen('Found fitted SN' + raw['SNID:'][0], user_choices)


def main():
    """Print documentation."""
    print __doc__

if __name__ == '__main__':
    main()
