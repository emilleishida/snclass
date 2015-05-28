"""
Created by Emille Ishida in May, 2015.

Miscelaneous functions for supernova classification.

- screen:
        Print messages to screen according to user choice

- check_reduction:
        Check dimensionality reduction function input choices.

- check_classifier:
        Check classifier function input choices.

- check_types:
        Check set types function input choices.

- check_crossval:
        Check cross-validation function input choices.

- read_user_input:
        read user choices from input file

- read_SNANA_lc:
        read raw SNANA light curve .DAT file

- choose_sn:
        Builds a list of SN satifying basic selction criteria

- read_fitted:
        reads previously calculated GP results
"""

import numpy as np
import os

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


def check_reduction(params):
    """
    Check dimensionality reduction function input choices.

    input: params, dict
           dictionary of input parameters

    output: params, dict
            updated dictionary of input parameters
    """
    if 'dim_reduction_func' in params.keys():
        if  params['dim_reduction_func'][0] == 'kpca':
            from snclass.functions import kpca
            params['dim_reduction_func'] = kpca

            for i in xrange(len(params['kpca_pars'])):
                par = params['kpca_pars'][i]
                try:
                    params[par] = int(params['kpca_val'][i])
                except ValueError:
                    try:
                        params[par] = float(params['kpca_val'][i])
                    except ValueError:
                        params[par] = params['kpca_val'][i]

        elif params['dim_reduction_func'][0] == 'None':
            params['dim_reduction_func'] = None

    return params


def check_classifier(params):
    """
    Check classifier function input choices.

    input: params, dict
           dictionary of input parameters

    output: params, dict
            updated dictionary of input parameters
    """
    if 'classifier_func' in params.keys():
        if params['classifier_func'][0] == 'nn':
            from snclass.functions import nn
            params['classifier_func'] = nn
            for i in xrange(len(params['classifier_pars'])):
                pvar = params['classifier_pars'][i]
                try:
                    params[pvar] = int(params['classifier_val'][i])
                except ValueError:
                    try:
                        params[pvar] = float(params['classifier_val'][i])
                    except ValueError:
                        params[pvar] = params['classifier_val'][i]

        elif params['classifier_func'][0] == 'None':
            params['classifier_func'][0] = None

    return params


def check_types(params):
    """
    Check set types function input choices.

    input: params, dict
           dictionary of input parameters

    output: params, dict
            updated dictionary of input parameters
    """
    if 'transform_types_func' in params.keys():
        if params['transform_types_func'][0] == 'set_types':
            from snclass.functions import set_types
            params['transform_types_func'] = set_types
        elif params['transform_types_func'][0] == 'None':
            params['transform_types_func'][0] = None

    return params


def check_crossval(params):
    """
    Check cross-validation function input choices.

    input: params, dict
           dictionary of input parameters

    output: params, dict
            updated dictionary of input parameters
    """
    if 'cross_validation_func' in params.keys():
        if params['cross_validation_func'][0] == 'cross_val':
            from snclass.functions import core_cross_val
            params['cross_validation_func'] = core_cross_val
            params['gamma_nparticles'] = int(params['gamma_nparticles'][0])

            name = 'n_cross_val_particles'
            params[name] = int(params[name][0])
            for i in xrange(len(params['cross_val_par'])):
                pvar = params['cross_val_par'][i] + '_lim'
                try:
                    params[pvar] = [int(params[pvar][j]) for j in range(2)]
                except ValueError:
                    params[pvar] = [float(params[pvar][j]) for j in range(2)]

        elif params['cross_validation_func'][0] == 'None':
            params['cross_validation_func'] = None

    return params


def read_user_input(filename):
    """
    Read user input from file and construct initial dictionary parameter.

    input:    filename (string) -> user input file parameter

    output:   dictionary with formated user choices
    """
    op1 = open(filename, 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [elem.split() for elem in lin1]

    # store options in params dictionary
    params = dict([(line[0], line[2:line.index('#')])
                   for line in data1 if len(line) > 1])

    params = check_reduction(params)
    params = check_classifier(params)
    params = check_crossval(params)

    params['GP_fit'] = {}
    params['realizations'] = {}
    params['xarr'] = {}
    params['GP_obj'] = {}
    params['GP_std'] = {}

    # check if ``observer'' data already exists
    if not os.path.isdir(params['path_to_obs'][0]):
        raise TypeError('Variable "path_to_obs" is not a valid directory!')

    return params


def read_SNANA_lc(params):
    """
    Read light curve in SNANA format and returns a dictionary.

    input:     params -> dictionary of input parameters

    output:    mdata -> data from light curve
    """

    # read light curve data
    op1 = open(params['path_to_obs'][0] + params['path_to_lc'][0], 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [elem.split() for elem in lin1]
    raw_data = dict([[line[0], line[1:]] for line in data1 if len(line) > 1])

    par = params['param_list'][0]

    # determine MJD index
    mjd_indx = raw_data[par].index(params['mjd_flag'][0]) + 1

    # determine filter index
    filter_indx = raw_data[par].index(params['filter_flag'][0]) + 1

    # determine photon count index
    photon_indx = raw_data[par].index(params['photon_flag'][0]) + 1

    # determine photon count error index
    photonerr_indx = raw_data[par].index(params['photonerr_flag'][0]) + 1

    #determine quality criteria index
    quality_indx = raw_data[par].index(params['quality_flag'][0]) + 1

    #build measurement list for each filter
    mdata = dict([[item, np.array([[float(line[ mjd_indx]),
                                    float(line[photon_indx]),
                                    float(line[photonerr_indx]),
                                    float(line[quality_indx])]
            for line in data1
            if len(line) > 1
            and line[0] == params['epoch_flag'][0]
            and line[filter_indx] == item
            and float(line[photon_indx]) >= 0.0
            and float(line[quality_indx]) >= \
                     float(params['quality_cut'][0])
                 ])]
            for item in params['filters']])

    #add usefull header information to output dictionary
    for item in params['header']:
        if item not in params['param_list']:
            mdata[item] = raw_data[item]

    return mdata


def choose_sn(params, output_file='snlist.dat'):
    """
    Read through all files within a given directory and
    select those satisfying criteria in user input file.

    input:  params (dict)

    output: txt file with the name of objects surviving selections cuts.
    """

    #take all file names in data directory
    filename = os.listdir(params['path_to_obs'][0])

    #store files for light curves surviving selection cuts
    final_list = []

    for name in filename:

        if params['file_root'][0] in name:

            #read light curve data
            op1 = open(params['path_to_obs'][0] + name, 'r')
            lin1 = op1.readlines()
            op1.close()

            data1 = [elem.split() for elem in lin1]

            #take header parameters
            header = {}
            for line in data1:
                if len(line) > 1 and line[0] in params['header']:
                    header[line[0]] = line[1:]

            #check for request SN type
            if  params['type_cut'][0]  != 'None':
                if header[params['type_flag'][0]][0] in params['type_cut']:
                    type_surv = True

                else:
                    type_surv = False
            else:
                type_surv = True

            #check for requested SN sample
            if params['sample_cut'][0] != 'None':

                if str(header[params['sample_flag'][0]][0]) in \
                              params['sample_cut']:
                    sample_surv = True
                else:
                    sample_surv = False

            else:
                sample_surv = True

            #store only if all requirements are satisfied
            if type_surv == True and sample_surv == True:
                final_list.append(name)

    op2 = open(output_file, 'w')
    for item in final_list:
        op2.write(item  + '\n')
    op2.close()

    screen('Found ' + str(len(final_list)) + \
           ' SN satisfying sample and type cuts.', params)
    screen('Surviving objects are listed in file ' + output_file, params)


def read_fitted(lc_data):
    """
    Read GP results and populate dictionary parameters.

    input:  user_input, dic
            output from read_SNANA_lc()

    output: updated dictionary of parameters.
    """
    loaded = {}

    if bool(int(lc_data['n_samples'][0])):
        op1 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + \
                   lc_data['SNID:'][0] + '_samples.dat', 'r')
        lin1 = op1.readlines()
        op1.close()

        d1 = [elem.split() for elem in lin1]

        loaded['realizations'] = {}
        loaded['xarr'] = {}
        for fil in lc_data['filters']:
            loaded['realizations'][fil] = [[float(d1[kk][jj]) 
                                            for kk in xrange(len(d1)) 
                                            if d1[kk][0]==fil]
                        for jj in xrange(2, int(lc_data['n_samples'][0]) + 2)]

            loaded['xarr'][fil] = []
            for i1 in xrange(len(d1)):
                if d1[i1][0] == fil:
                    loaded['xarr'][fil].append(float(d1[i1][1]))

    op2 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + \
               lc_data['SNID:'][0] + '_mean.dat', 'r')
    lin2 = op2.readlines()
    op2.close()

    d2 = [elem.split() for elem in lin2]

    loaded['GP_std'] = {}
    loaded['GP_fit'] = {}
    loaded['xarr'] = {}
    for fil in lc_data['filters']:
        loaded['xarr'][fil] = [float(d2[j][1]) for j in xrange(1,len(d2))
                                               if d2[j][0] == fil]

        loaded['GP_fit'][fil] = [float(d2[j][2]) for j in xrange(1,len(d2))
                                                 if d2[j][0] == fil]

        loaded['GP_std'][fil] = [float(d2[j][3]) for j in xrange(1,len(d2))
                                                 if d2[j][0] == fil]

    return loaded


def main():
    print(__doc__)

if __name__=='__main__':
    main()    
