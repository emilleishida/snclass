"""
Created by Emille Ishida in May, 2015.

Function for performing Gaussian Process fit using gptools.

- imp_gptools:
         Perform Gaussian Process with gptools through MCMC.

- fit_LC:
         Gaussian Process fit using gptools.
"""

import numpy as np
import gptools
import os


def imp_gptools(data, fil, mcmc=False):
    """
    Perform Gaussian Process with gptools through MCMC.

    input: data, dict
           dictionary of raw data
           output from read_snana_lc
           keys: filters

           fil, str
           filter
        
           mcmc, bool, optional
           if True, optimize kernel parameters using mcmc
           Default is False

    output: data, dict
            updated dictionary with GP results
    """
    # format data
    mjd = data[fil][:, 0]
    flux = data[fil][:, 1]
    fluxerr = data[fil][:, 2]

    # setup GP
    k_obj = gptools.SquaredExponentialKernel(param_bounds=[(0, max(flux)),
                                             (0, np.std(mjd))])
    gp_obj = gptools.GaussianProcess(k_obj)
    gp_obj.add_data(mjd, flux, err_y=fluxerr)

    data['xarr'][fil] = np.arange(min(mjd), max(mjd), 0.2)

    if mcmc:
        out = gp_obj.predict(data['xarr'][fil], use_MCMC=True,
                             num_proc=int(data['n_proc'][0]), nsamp=200,
                             plot_posterior=False,
                             plot_chains=False, burn=100, thin=10)

        data['GP_obj'][fil] = gp_obj 

    else:
        out = gp.predict(x_star, use_MCMC=False)
        data['GP_obj'] = gp.k.params

    data['GP_fit'][fil] = out[0]
    data['GP_std'][fil] = out[1]    

    return data


def save_result(data, mean=True, samples=False):
    """
    Save results of GP fit to file.

    input: data, dict
           dictionary of raw data
           output from read_snana_lc
           keys: filters

           mean, bool - optional
           if True, save mean GP fit
           Default is True

           samples, bool - optional
           if True, save draws from GP fit
           Default is False
    """
    # check if storage directory exsts
    if not os.path.exists(data['samples_dir'][0]):
        os.makedirs(data['samples_dir'][0])

    if bool(int(data['save_samples'][0])) and samples:
        op1 = open(data['samples_dir'][0] + data['file_root'][0] +
                   data['SNID:'][0] + '_samples.dat', 'w')
        op1.write('filter    MJD    ')
        for j in xrange(int(data['n_samples'][0])):
            op1.write('samp' + str(j + 1))
        op1.write('\n')
        for fil in data['filters']:
            for i in xrange(len(data['xarr'][fil])):
                op1.write(fil + '    ' +
                          str(data['xarr'][fil][i]) + '    ')
                for j in xrange(int(data['n_samples'][0])):
                    op1.write(str(data['realizations'][fil][j][i]) +
                              '    ')
                op1.write('\n')
        op1.close()

    if mean:
        op2 = open(data['samples_dir'][0] + data['file_root'][0] +
                   data['SNID:'][0] + '_mean.dat', 'w')
        op2.write('filter    MJD    GP_fit     GP_std\n')
        for fil in data['filters']:
            for k in xrange(len(data['xarr'][fil])):
                op2.write(fil + '    ' + str(data['xarr'][fil][k]) +
                          '    ' + str(data['GP_fit'][fil][k]) +
                          '    ' + str(data['GP_std'][fil][k]) + '\n')
        op2.close()


def fit_lc(data, mean=True, samples=False, screen=False, do_mcmc=False):
    """
    Gaussian Process fit using gptools.

    input:  data -> dictionary of raw data
                    output from read_snana_lc
                    keys: filters

            mean -> bool, optional
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
                    Default is False

    output: data -> update dictionary with new keyword:
                    realizations
    """
    key_list = ['realizations', 'xarr', 'GP_std', 'GP_obj']

    for name in key_list:
        if name not in data.keys():
            data[name] = {}

    for fil in data['filters']:
        if screen:
            print '... filter: ' + fil

        if mean:
            data = imp_gptools(data, fil, mcmc=do_mcmc)

        if samples and int(data['n_samples'][0]) > 0:
            if screen:
                print '... ... calculate samples'

            new_obj = data['GP_obj'][fil]
            draws = new_obj.draw_sample(data['xarr'][fil],
                                        num_samp=int(data['n_samples'][0]))

            data['realizations'][fil] = draws.T

    save_result(data, mean=mean, samples=samples)

    if screen:
        print '\n'

    return data


def main():
    """Print documentation."""
    print __doc__

if __name__ == '__main__':
    main()
