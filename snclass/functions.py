#!/usr/bin/env python

from __future__ import division

import numpy as np
import matplotlib.pylab as plt
from itertools import combinations
from scipy.stats import uniform

from sklearn.decomposition import KernelPCA
from sklearn import neighbors

from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence

import snclass
from snclass.util import screen

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

    kpca = KernelPCA(kernel=pars['kernel'], gamma=pars['gamma'], n_components=pars['ncomp'])
    X_kpca = kpca.fit_transform(data_matrix)

    if transform:
        return kpca
    else:
        return X_kpca

def nn(test, data_matrix, types, pars):
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

    #initia NN object
    clf = neighbors.KNeighborsClassifier(pars['n'], weights=pars['weights'])

    #fit model
    clf.fit(data_matrix, types)

    #predict type
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
    np.random.seed()

    #split sample in 3
    indx_list1 = [np.random.randint(0, len(data)) for j in xrange(int(2*len(data)/3))]
    indx_list2 = [elem for elem in xrange(len(data)) if (elem not in indx_list1)]

    #set train data matrix and types
    d2 = snclass.DataMatrix()
    d2.user_choices = user_choices
    d2.datam = np.array([data[indx] for indx in indx_list1]) 
    d2.sntype = np.array([types[indx] for indx in indx_list1])

    #set test data matrix and types
    d2.data_test = np.array([data[indx] for indx in indx_list2]) 
    test_type = np.array([types[indx] for indx in indx_list2])

    ploc = d2.user_choices['gamma_lim'][0]
    pscale = d2.user_choices['gamma_lim'][1] - ploc
    dist = uniform(loc=ploc, scale=pscale)

    results = []
    for ncomp in xrange(2, 11):

        screen('... ncomp = ' + str(ncomp), user_choices)
        for k in xrange(user_choices['gamma_nparticles']):

            try:
                #reduce dimensionality
                d2.user_choices['gamma'] = dist.rvs()
                d2.user_choices['ncomp'] = ncomp

                screen('... ... gamma = ' + str(d2.user_choices['gamma']), user_choices)
    
                d2.reduce_dimension()

                #project test 
                test_proj = d2.transf_test.transform(d2.data_test)

                #classify
                new_label = nn(test_proj, d2.low_dim_matrix, d2.sntype, d2.user_choices)

                #calculate score
                score = sum(new_label == test_type)

                #store
                results.append([ncomp, d2.user_choices['gamma'], score])

            except ArpackNoConvergence:
                pass

    results = np.array(results)
    indx_max = list(results[:,-1]).index(max(results[:,-1]))
    
    return results[indx_max]


"""
import numpy as np
import snclass
from snclass.util import read_user_input

path = '/home/emille/Dropbox2/Dropbox/meu_KPCA/artigo_type_Ia_n2/calculations/post-SNPCC/output/matrix/griz/extrapolate_no/GP/SNR0_3_24_s1/data_matrix_spec.dat'

op1 = open(path, 'r')
lin1 = op1.readlines()
op1.close()

data1 = [elem.split() for elem in lin1]

datam=np.array([[float(elem) for elem in line[4:]] for line in data1[1:]])
type = np.array([line[1] for line in data1[1:]])

user_input=snclass.util.read_user_input('fit_lc_input.dat')
d=snclass.DataMatrix()
d.datam = datam
d.sntype = type
d.user_choices = user_input
d.cross_val()

y = np.array(d.sntype[1:])

snIa = y == '0'
snother = y != '0'


plt.figure()
plt.scatter(matrix[snIa, 0], matrix[snIa, 1], color='red', label='Ia')
plt.scatter(matrix[snother, 0], matrix[snother, 1], color='blue', label='nonIa')
plt.scatter([item[0][0]], [item[0][1]], color='green', label=new_label)
plt.legend()
plt.show()
"""
