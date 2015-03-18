#!/usr/bin/env python

from __future__ import division

import argparse
import emcee
import matplotlib.pyplot as plt
import numpy as np
import sys

import george
from george import kernels

from snclass.util import read_user_input, read_SNANA_lc
from snclass.fit_lc_george import lnprob2, fit_LC


def main(args):

    #read_user_input
    user_input = read_user_input(args.input) 

    #read lc data
    lc_data = read_SNANA_lc(user_input)

    #add extra keys
    lc_data.update(user_input)

    if bool(int(args.calculate)):
        #fit lc
        lc_data = fit_LC(lc_data)
    else:
        op1 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + lc_data['SNID:'][0] + '_samples.dat', 'r')
        lin1 = op1.readlines()
        op1.close()

        d1 = [elem.split() for elem in lin1]

        lc_data['xarr_indx'] = {}
        for fil in lc_data['filters']:
            lc_data['xarr'][fil] = []
            lc_data['xarr_indx'][fil] = []
            lc_data['realizations'][fil] = []
            for item in d1[0]:
                if fil == item[0]:
                    lc_data['xarr'][fil].append(float(item[1:]))
                    lc_data['xarr_indx'][fil].append(d1[0].index(item))

        for line in d1[1:]:
            for fil in lc_data['filters']:
                lc_data['realizations'][fil].append([float(line[i]) for i in lc_data['xarr_indx'][fil]])        

        op2 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + lc_data['SNID:'][0] + '_samples.dat', 'r')
        lin2 = op2.readlines()
        op2.close()

        d2 = [elem.split() for elem in lin2]

        for fil in lc_data['filters']:
            lc_data['GP_fit'][fil] = [float(d2[1][j]) for j in lc_data['xarr_indx'][fil]]

    #initiate figure
    plt.figure()
    
    for fil in user_input['filters']:

        # Plot the samples in data space.
        ax = plt.subplot(len(user_input['filters']), 1, user_input['filters'].index(fil) + 1)
        for s in lc_data['realizations'][fil]:
            plt.plot(lc_data['xarr'][fil], s, color="#4682b4", alpha=0.3)
        plt.errorbar(lc_data[fil][:,0], lc_data[fil][:,1], yerr=lc_data[fil][:,2], fmt=".k", capsize=0, label=fil)
        plt.plot(lc_data['xarr'][fil], lc_data['GP_fit'][fil], 'r:', linewidth=2)
        plt.ylabel("FLUXCAL")
        plt.xlabel("MJD")
        plt.legend()
        plt.xlim(min(lc_data[fil][:,0]) - 1.0, max(lc_data[fil][:,0]) + 1.0)

    plt.suptitle("George - results with Gaussian process noise model")
    plt.savefig("gp-results.png", dpi=350)
    plt.close()


if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description = 'Supernova photometric classification using KPCA.')
    parser.add_argument('-i','--input', help = 'Input file name', required = True)
    parser.add_argument('-c', '--calculate', help = 'Read or calculate GP fit', required = True) 
    args = parser.parse_args()
   
    main(args)

