"""
Created by Emille Ishida in May, 2015.

Stand alone functions for supernova classification.

- screen:
        Print messages to screen according to user choice.

- kpca:
        Perform dimensionality reduction using kernel PCA.

- nneighbor:
        Classify a given data matrix according to its n nearst neighbours.

- set_types:
        Transform the original vector of types.

- calc_scores:
        Calculate classification results for 1 data matrix.

- core_cross_val:
        Perform 1/3 validation.
"""

from __future__ import division

import numpy as np
from scipy.stats import uniform

from sklearn.decomposition import KernelPCA
from sklearn import neighbors

from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence

import snclass

#########################################


def screen(message, choices):
    """
    Print message on screen according to users choice.

    input:   message, str
             message to be printed

             choices, dict
             dictionary of users choices
    """
    if bool(int(choices['screen'][0])):
        print message


def kpca(data_matrix, pars, transform=False):
    """
    Perform dimensionality reduction using kernel PCA.

    input: data_matrix, array

           other arguments are passed directly to
           sklearn.decomposition.KernelPCA rotine.

           pars, dict
           dictionary of parameters.
           keywords: 'kernel', 'gamma'

           transform, bool
           if True return the kpca object for further projections.
           Default is False.

    output: X_kpca, array
            lines are objects.
            collumns are projections over different kPCs.
    """
    obj_kpca = KernelPCA(kernel=pars['kernel'], gamma=pars['gamma'],
                         n_components=pars['ncomp'])
    x_kpca = obj_kpca.fit_transform(data_matrix)

    if transform:
        return obj_kpca
    else:
        return x_kpca


def nneighbor(test, data_matrix, types, pars):
    """
    Classify a given data matrix according to its first nearst neighbour.

    input: test, vector (array)
           coordinates of object(s) to be classified

           data_matrix, array
           full data matrix for training

           types, vector of str
           types for each element on data_matrix

           pars - dict
           Dictionary of parameters
           keywords: 'n', 'weights'

    output: type, list of str
            classification of point
    """
    # initia NN object
    clf = neighbors.KNeighborsClassifier(pars['n'], weights=pars['weights'])

    # fit model
    clf.fit(data_matrix, types)

    # predict type
    new_label = clf.predict(test)

    return new_label


def set_types(types):
    """
    Transform the original vector of types.

    All non-Ia objects are typed '1' and all Ia objects are typed as '0'.
    This confusing nomenclature is to abide with SNANA notation.

    input: type - vector of str
           type of each object in the data matrix

    output: new_type - vector of str
            translation to appropriate types
    """
    new_type = []
    for item in types:
        if item == '0':
            new_type.append(item)
        else:
            new_type.append('1')

    return np.array(new_type)


def calc_scores(matrix2, ncomp, dist):
    """
    Calculate classification results for 1 data matrix.

    input: matrix2, DataMatrix object
           output from DataMatrix.build()

           ncomp, int
           number of PCs to calculate

           dist, scipy.stats.uniform distribution
           prior over gamma parameter
    """
    np.random.seed()

    # reduce dimensionality
    matrix2.user_choices['gamma'] = dist.rvs()
    matrix2.user_choices['ncomp'] = ncomp

    matrix2.reduce_dimension()

    # project test
    test_proj = matrix2.transf_test.transform(matrix2.data_test)

    # classify
    new_label = nneighbor(test_proj, matrix2.low_dim_matrix,
                          matrix2.sntype, matrix2.user_choices)

    # calculate score
    score = sum(new_label == matrix2.test_type)

    return int(ncomp), matrix2.user_choices['gamma'], score


def core_cross_val(data, types, user_choices):
    """
    Perform 1/3 validation.

    input: data, array
           data matrix

           types, vector
           vector of types

           user_choices, dict
           output from read_user_input()

    output: vector of floats
            parameters with higher classification success
            [n_components, gamma, n_successes]
    """
    # split sample in 3
    indx_list1 = np.random.randint(0, len(data), size=int(2 * len(data) / 3))
    indx_list2 = [elem for elem in xrange(len(data))
                  if elem not in indx_list1]

    # set train data matrix and types
    matrix2 = snclass.matrix.DataMatrix()
    matrix2.user_choices = user_choices
    matrix2.datam = np.array([data[indx] for indx in indx_list1])
    matrix2.sntype = np.array([types[indx] for indx in indx_list1])

    # set test data matrix and types
    matrix2.data_test = np.array([data[indx] for indx in indx_list2])
    matrix2.test_type = np.array([types[indx] for indx in indx_list2])

    ploc = matrix2.user_choices['gamma_lim'][0]
    pscale = matrix2.user_choices['gamma_lim'][1] - ploc
    dist = uniform(loc=ploc, scale=pscale)

    results = []
    for ncomp in xrange(2, 11):

        screen('... ncomp = ' + str(ncomp), user_choices)

        k = 0
        while k < user_choices['gamma_nparticles']:
            try:
                results.append(calc_scores(matrix2, ncomp, dist))

                # update counter
                k = k + 1

            except ArpackNoConvergence:
                screen('Arparck fail to converge!', user_choices)

    results = np.array(results)
    indx_max = list(results[:, -1]).index(max(results[:, -1]))

    return results[indx_max]

