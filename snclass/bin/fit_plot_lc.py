# Copyright 2015 Emille Ishida
# This program is distributed under the terms of the GNU General Purpose License (GPL).
# Refer to http://www.gnu.org/licenses/gpl.txt
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
Fit and plot a light curve using Gaussian Process. 

usage: 

    In order to fit a GP and plot the result, do

    $ fit_plot_lc.py -i <user.input> -c 1

    in case you are only interested in plotting a previously calculated
    result, 

    $ fit_plot_lc.py -i <user.input> -c 0
"""

#!/usr/bin/env python

from __future__ import division

import argparse
import matplotlib.pyplot as plt
import numpy as np

from snclass.util import read_user_input, read_snana_lc
from snclass.fit_lc_gptools import fit_lc
from snclass.functions import screen


def main(args):
    """Read user input, fit and plot a GP and the raw data."""
    # read_user_input
    user_input = read_user_input(args.input) 

    # read lc data
    lc_data = read_snana_lc(user_input)

    # add extra keys
    lc_data.update(user_input)

    # set screen output
    out = bool(int(user_input['screen'][0]))

    if user_input['measurement'][0] == 'flux':     
        ylabel = 'flux'
        sign = 1.0
    else:
        ylabel = 'magnitude'
        sign = -1.0

    if bool(int(args.calculate)):
       
        screen('Fitting SN' + lc_data['SNID:'][0], user_input)

        if user_input['measurement'][0] == 'flux':
            p1 = [int(user_input['epoch_predict'][0]), 
                  int(user_input['epoch_predict'][1])]
            sign2 = 1.0

        else:
            p1 = None
            sign2 = -1.0

        # fit lc
        lc_data = fit_lc(lc_data, samples=bool(int(lc_data['n_samples'][0])),
                         save_samples=bool(int(user_input['save_samples'][0])),
                         screen=out, 
                         do_mcmc=bool(int(user_input['do_mcmc'][0])),
                         predict=p1)
    else:
        sign2 = 1.0 

        if bool(int(lc_data['n_samples'][0])):
            op1 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + \
                       lc_data['SNID:'][0] + '_' + user_input['measurement'][0] + '_samples.dat', 'r')
            lin1 = op1.readlines()
            op1.close()

            d1 = [elem.split() for elem in lin1]


            for fil in lc_data['filters']:
                lc_data['xarr'][fil] = []

                if bool(int(lc_data['n_samples'][0])):
                    lc_data['realizations'][fil] = [[float(d1[kk][jj]) 
                                                     for kk in xrange(len(d1)) 
                                                     if d1[kk][0]==fil] 
                                                     for jj in xrange(2, 
                                                     int(lc_data['n_samples'][0]) + 2)]
 
            for i1 in xrange(len(d1)):
                if d1[i1][0] == fil:
                    lc_data['xarr'][fil].append(float(d1[i1][1]))
                
        op2 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + \
                   lc_data['SNID:'][0] + '_' + user_input['measurement'][0] + '_mean.dat', 'r')
        lin2 = op2.readlines()
        op2.close()

        d2 = [elem.split() for elem in lin2]

        lc_data['GP_std'] = {}
        for fil in lc_data['filters']:
            lc_data['xarr'][fil] = []
            lc_data['GP_fit'][fil] = np.array([float(d2[j][2]) 
                                               for j in xrange(1,len(d2)) if d2[j][0] == fil])
            lc_data['GP_std'][fil] = np.array([float(d2[j][3]) 
                                               for j in xrange(1,len(d2)) if d2[j][0] == fil])
            lc_data['xarr'][fil] = np.array([float(d2[j][1])
                                             for j in xrange(1,len(d2)) if d2[j][0] == fil])

    #initiate figure
    f = plt.figure()
    print 'sign = ' + str(sign)
    for fil in user_input['filters']:

        # Plot the samples in data space.
        plt.subplot(2, len(lc_data['filters'])/2 + 
                        len(lc_data['filters'])%2, 
                        lc_data['filters'].index(fil) + 1)
        if bool(int(lc_data['n_samples'][0])):
            for s in lc_data['realizations'][fil]:
                plt.plot(lc_data['xarr'][fil], sign2 * np.array(s), color="gray", alpha=0.3)
        plt.errorbar(lc_data[fil][:,0], sign * lc_data[fil][:,1], 
                     yerr=lc_data[fil][:,2], fmt="o", color='blue', label=fil)
        plt.plot(lc_data['xarr'][fil], sign2 * lc_data['GP_fit'][fil], 
                 color='red', linewidth=2)
        plt.ylabel(ylabel)
        plt.xlabel("MJD")
        plt.legend()
        plt.xlim(min(lc_data['xarr'][fil]) - 1.0, max(lc_data['xarr'][fil]) + 1.0)
        if user_input['measurement'][0] == 'mag':
            plt.ylim(min(sign * lc_data[fil][:,1]) - 1.5*max(lc_data[fil][:,2]),max(sign * lc_data[fil][:,1]) + 1.5*max(lc_data[fil][:,2]))  
            ax = plt.gca()
            ax.invert_yaxis()

    f.tight_layout()
    plt.savefig("gp-SN" + lc_data['SNID:'][0] + "_" + user_input['measurement'][0] + ".png", dpi=350)
    plt.close()


if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Supernova photometric ' + \
                                     'classification using KPCA.')
    parser.add_argument('-i','--input', help='Input file name', 
                        required = True)
    parser.add_argument('-c', '--calculate', help='Read or calculate GP fit',
                        required=True) 
    args = parser.parse_args()
   
    main(args)

