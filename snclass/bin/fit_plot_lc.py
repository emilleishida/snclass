"""
Created by Emille Ishida in May, 2015.

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

    if bool(int(args.calculate)):
       
        screen('Fitting SN' + lc_data['SNID:'][0], user_input)

        # fit lc
        lc_data = fit_lc(lc_data, samples=bool(int(lc_data['n_samples'][0])),
                         screen=out)
    else:

        if bool(int(lc_data['n_samples'][0])):
            op1 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + \
                       lc_data['SNID:'][0] + '_samples.dat', 'r')
            lin1 = op1.readlines()
            op1.close()

            d1 = [elem.split() for elem in lin1]

        for fil in lc_data['filters']:
            lc_data['xarr'][fil] = []
            lc_data['realizations'][fil] = [[float(d1[kk][jj]) 
                                             for kk in xrange(len(d1)) 
                                             if d1[kk][0]==fil] 
                                             for jj in xrange(2, 
                                             int(lc_data['n_samples'][0]) + 2)]
 
            for i1 in xrange(len(d1)):
                if d1[i1][0] == fil:
                    lc_data['xarr'][fil].append(float(d1[i1][1]))
                
        op2 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + \
                   lc_data['SNID:'][0] + '_mean.dat', 'r')
        lin2 = op2.readlines()
        op2.close()

        d2 = [elem.split() for elem in lin2]

        lc_data['GP_std'] = {}
        for fil in lc_data['filters']:
            lc_data['GP_fit'][fil] = [float(d2[j][2]) 
                                      for j in xrange(1,len(d2)) if d2[j][0] == fil]
            lc_data['GP_std'][fil] = [float(d2[j][3]) 
                                      for j in xrange(1,len(d2)) if d2[j][0] == fil]

    #initiate figure
    f = plt.figure()
    
    for fil in user_input['filters']:

        # Plot the samples in data space.
        plt.subplot(2, len(lc_data['filters'])/2 + 
                        len(lc_data['filters'])%2, 
                        lc_data['filters'].index(fil) + 1)
        for s in lc_data['realizations'][fil]:
            plt.plot(lc_data['xarr'][fil], s, color="gray", alpha=0.3)
        plt.errorbar(lc_data[fil][:,0], lc_data[fil][:,1], 
                     yerr=lc_data[fil][:,2], fmt="o", color='blue', label=fil)
        plt.plot(lc_data['xarr'][fil], lc_data['GP_fit'][fil], 
                 color='red', linewidth=2)
        plt.ylabel("FLUXCAL")
        plt.xlabel("MJD")
        plt.legend()
        plt.xlim(min(lc_data[fil][:,0]) - 1.0, max(lc_data[fil][:,0]) + 1.0)

    f.tight_layout()
    plt.savefig("gp-SN" + lc_data['SNID:'][0] + ".png", dpi=350)
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

