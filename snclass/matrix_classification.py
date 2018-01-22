"""
Created by Emille Ishida in June 2015.

Functions related to matrix classification and KernelPCA projectons.

- store_test_matrix:
        Store results for different filters in one matrix.
- test_samples:
        Test if samples satisfy selection criterias
- classify_test:
        Classify one photometric supernova using a trained KernelPCA matrix.
- plot_matrix:
        Plot training KernelPCA matrix and eventual test object.
"""

import numpy as np
from snclass.matrix import DataMatrix
from snclass.treat_lc import LC
from snclass.util import read_user_input, read_snana_lc, translate_snid
from snclass.functions import nneighbor

import shutil
import os

##############################################################################


def store_test_matrix(new_lc):
    """
    Store results for different filters in one matrix.

    input: new_lc, snclass.treat_lc.LC obj
           treated light curve with low dimensional representation

    output: samples data matrix
    """
    # initiate matrix
    test_matrix = []
    ini_fil = new_lc.user_choices['filters'][0]

    # concatenate filters
    for j in xrange(len(new_lc.samples_for_matrix)):
        matrix_lines = []
        for item in new_lc.samples_for_matrix[j]:
            matrix_lines.append(item)
        test_matrix.append(matrix_lines)

    return np.array(test_matrix)

def store_test_mean(new_lc):
    """
    Store results for different filters in one matrix.

    input: new_lc, snclass.treat_lc.LC obj
           treated light curve with low dimensional representation

    output: samples data matrix
    """

    # concatenate filters
    matrix_lines = []
    for j in new_lc.raw['filters']:
        for item in new_lc.flux_for_matrix[j]:
            matrix_lines.append(item)

    return np.array([matrix_lines])


def test_samples(new_lc, calc_samples=True):
    """
    Test if samples satisfy selection criterias.

    input: new_lc, snclass.treat_lc.LC object
           information from test object

           calc_samples, bool, optional
           If True, fit GP object and generate sample file as output
           otherwise reads samples from file
           Default is True
    """
    if calc_samples:
        new_lc.fit_GP(mean=True, samples=True, screen=True,
                      save_mean=True, save_samples=True)
        new_lc.normalize(samples=True)
        new_lc.mjd_shift()
        new_lc.build_steps(samples=True)

        #store test matrix in array
        test_matrix = store_test_matrix(new_lc)

    else:
        new_lc.build_steps(samples=calc_samples)

        #store test matrix in array
        test_matrix = store_test_mean(new_lc)

    return test_matrix


def classify_test(test_name, matrix, user_input, test_dir='test_samples/',
                  csamples=True):
    """
    Classify one photometric supernova using a trained KernelPCA matrix.

    input: test_name, str
           name of mean GP fit file

           matrix, snclass.matrix.DataMatrix object
           trained KernelPCA matrix

           user_input, dict
           output from snclass.util.read_user_input

           test_dir, str, optional
           name of directory to store samples from test light curve
           Default is 'test_samples/'

           csamples, bool, optional
           If True, fit GP object and generate sample file as output
           otherwise reads samples from file
           Default is True

    return: new_lc, snclass.treat_lc.LC object
            updated with test projections and probability of being Ia
    """
    # update path to raw light curve
    user_input['path_to_lc'] = [translate_snid(test_name, 'FLUXCAL')[0]]

    # store number of samples for latter tests
    nsamples = user_input['n_samples'][0]

    # reset the number of samples for preliminary tests
    user_input['n_samples'] = ['0']

    # read raw data
    raw = read_snana_lc(user_input)

    # load GP fit and test epoch cuts
    new_lc = LC(raw, user_input)
    new_lc.load_fit_GP(user_input['samples_dir'][0] + test_name)
    new_lc.normalize()
    new_lc.mjd_shift()
    new_lc.check_epoch()

    if new_lc.epoch_cuts:
        # update test sample directory
        user_input['samples_dir'] = [test_dir]

        # update user choices
        new_lc.user_choices = user_input

        # update number of samples
        new_lc.user_choices['n_samples'] = [nsamples]

        # fit GP or normalize/shift fitted mean
        test_matrix = test_samples(new_lc, calc_samples=bool(csamples))   

        # project test
        new_lc.test_proj = matrix.transf_test.transform(test_matrix)
      
        # classify
        new_lc.new_label = nneighbor(new_lc.test_proj, matrix.low_dim_matrix,
                              matrix.sntype, matrix.user_choices)

        if csamples:
            new_lc.prob_Ia = sum([1 for item in new_label 
                                  if item == '0'])/float(nsamples)

        return new_lc

    else:
        return None

