"""
Created by Emille Ishida in May, 2015.

Miscelaneous functions for supernova classification.


- translate_snid:
    Translate a general sn id between SNANA and raw SNID format.

- check_reduction:
        Check dimensionality reduction function input choices.

- check_classifier:
        Check classifier function input choices.

- check_types:
        Check set types function input choices.

- check_crossval:
        Check cross-validation function input choices.

- read_user_input:
        Read user choices from input file

- build_indx:
        Build dictionary of index for header variables.

- read_snana_lc:
        Read raw snana light curve .DAT file

- compare_type:
        Compare type in user requests with raw SN data.

- choose_sn:
        Builds a list of SN satifying basic selction criteria

- read_fitted:
        Read previously calculated GP results.
"""

import numpy as np
import os

from snclass.functions import screen


def translate_snid(snid, meas):
    """
    Translate a general sn id between SNANA and raw SNID format.

    input: snid, str
           name of file or other string containing sind
 
           meas, str
           type of measurment (e.g. flux, mag, photon)

    output: raw snana file name and/or shorten snid
    """
    if 'mean' in snid:
        name1 = snid[:-len('_' + meas + '_mean.dat')]
    elif 'samples' in snid:
        name1 = snid[:-len('_' + meas + '_samples.dat')]
    elif '.DAT' in snid:
        name1 = snid[:-len('.DAT')]
    else:
        name1 = snid

    if 'X' in name1:
        name2 = name1[name1.index('X') + 1:]
    elif 'DES_SN' in snid:
        name2 = name1[len('DES_SN'):]
    elif 'SDSS_SN' in snid:
        name2 = name1[len('SDSS_SN'):]
    else:
        name2 = name1

    if snid[:len('DES_SN')] == 'DES_SN':
        if snid[-len('.DAT'):] == '.DAT' and len(name2) == 6:
            if '0' in name2:
                name3 = name2.lstrip('0')
            else:
                name3 = name2

            return name3

        else:
            if len(name2) == 6:
                name3 = name2
            elif len(name2) == 5:
                name3 = '0' + name2
            elif len(name2) == 4:
                name3 = '00' + name2
            elif len(name2) == 3:
                name3 = '000' + name2
            elif len(name2) == 2:
                name3 = '0000' + name2
            else:
                print snid
                raise NameError('File name not following convention!')

            return 'DES_SN' + name3 + '.DAT', name2

    elif snid[:len('SDSS_SN')] == 'SDSS_SN':
        if snid[-len('.DAT'):] == '.DAT' and len(name2) == 8:
            if '0' in name2:
                name3 = name2.lstrip('0')
            else:
                name3 = name2

            return name3

        else:
            if len(name2) == 8:
                name3 = name2
            elif len(name2) == 7:
                name3 = '0' + name2
            elif len(name2) == 6:
                name3 = '00' + name2
            elif len(name2) == 5:
                name3 = '000' + name2
            elif len(name2) == 4:
                name3 = '0000' + name2
            elif len(name2) == 3:
                name3 = '00000' + name2
            elif len(name2) == 2:
                name3 = '000000' + name2

            return 'SDSS_SN' + name3 + '.DAT', name2


def check_reduction(params):
    """
    Check dimensionality reduction function input choices.

    input: params, dict
           dictionary of input parameters

    output: params, dict
            updated dictionary of input parameters
    """
    if 'dim_reduction_func' in params.keys():
        if params['dim_reduction_func'][0] == 'kpca':
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
        if params['classifier_func'][0] == 'nneighbor':
            from snclass.functions import nneighbor
            params['classifier_func'] = nneighbor
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
    params = check_types(params)

    params['GP_fit'] = {}
    params['realizations'] = {}
    params['xarr'] = {}
    params['GP_obj'] = {}
    params['GP_std'] = {}

    # check if ``observer'' data already exists
    if not os.path.isdir(params['path_to_obs'][0]):
        raise TypeError('Variable "path_to_obs" is not a valid directory!')

    return params


def build_indx(params, raw):
    """
    Build dictionary of index for header variables.

    input: params, dict
           dictionary of user choices

           raw, dict
           dictionary from raw light curve

    output: indx, dict
            dictionary of index for header variables
    """
    par = params['param_list'][0]

    indx = {}

    # determine MJD index
    indx['mjd'] = raw[par].index(params['mjd_flag'][0]) + 1

    # determine filter index
    indx['filter'] = raw[par].index(params['filter_flag'][0]) + 1

    # determine photon count index
    indx['photon'] = raw[par].index(params['photon_flag'][0]) + 1

    # determine photon count error index
    indx['photonerr'] = raw[par].index(params['photonerr_flag'][0]) + 1

    # determine quality criteria index
    indx['quality'] = raw[par].index(params['quality_flag'][0]) + 1

    return indx


def read_snana_lc(params):
    """
    Read light curve in snana format and returns a dictionary.

    input:     params -> dictionary of input parameters

    output:    mdata -> data from light curve
    """
    # read light curve data
    op1 = open(params['path_to_obs'][0] + params['path_to_lc'][0], 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [elem.split() for elem in lin1]
    raw_data = dict([[line[0], line[1:]] for line in data1 if len(line) > 1])

    # find correct indexes for header parameters
    pindx = build_indx(params, raw_data)

    # build measurement list for each filter
    fils = params['filters']

    if params['photon_flag'][0] == 'mag':
        sign = -1
    else:
        sign = 1
    
    mdata = {}
    for fil in fils:
        l1 = []
        for line in data1:
            if  len(line) > 1 and line[0] == params['epoch_flag'][0] and line[pindx['filter']] == fil and float(line[pindx['quality']]) >= float(params['quality_cut'][0]):
                if params['photon_flag'][0] == 'mag' and float(line[pindx['photonerr_flag']]) < 50.0 and float(line[pindx['photonerr']]) < 50.0:
                    l1.append([float(line[pindx['mjd']]),
                                    -float(line[pindx['photon']]),
                                    float(line[pindx['photonerr']]),
                                    float(line[pindx['quality']])])
                elif params['photon_flag'][0] == 'FLUXCAL' and float(line[pindx['photon']]) > 0:
                    l1.append([float(line[pindx['mjd']]),
                                    float(line[pindx['photon']]),
                                    float(line[pindx['photonerr']]),
                                    float(line[pindx['quality']])])

        mdata[fil] = np.array(l1)

    # add usefull header information to output dictionary
    for item in params['header']:
        if item not in params['param_list']:
            mdata[item] = raw_data[item]

    return mdata


def compare_type(params, header):
    """
    Compare type in user requests with raw SN data.

    input:  params, dict
            dictionary of user choices

            header, dict
            dictionary of header variables

    output: type_surv, bool
            if True, object satisfies type selection cuts
    """
    if params['type_cut'][0] != 'None':
        if header[params['type_flag'][0]][0] in params['type_cut']:
            type_surv = True
        else:
            type_surv = False
            print 'params[type_flag][0] = ' + params['type_flag'][0]
            print 'params[type_cut] = ' + str(params['type_cut'])
            print 'type_surv = ' + str(type_surv)
    else:
        type_surv = True

    return type_surv


def compare_sample(params, header):
    """
    Compare sample in user requests with raw SN data.

    input:  params, dict
            dictionary of user choices

            header, dict
            dictionary of header variables

    output: sample_surv, bool
            if True, object satisfies sample selection cuts
    """
    if params['sample_cut'][0] != 'None':
        par = params['sample_flag'][0]
        if str(header[par][0]) in params['sample_cut']:
            sample_surv = True
        else:
            sample_surv = False
    else:
        sample_surv = True

    return sample_surv


def choose_sn(params, output_file='snlist.dat'):
    """
    Select objects satisfying criterias in user input file.

    input:  params (dict)

    output: txt file with the name of objects surviving selections cuts.
    """
    # take all file names in data directory
    filename = os.listdir(params['path_to_obs'][0])

    # store files for light curves surviving selection cuts
    final_list = []

    for name in filename:

        if params['file_root'][0] in name:

            # read light curve data
            op1 = open(params['path_to_obs'][0] + name, 'r')
            lin1 = op1.readlines()
            op1.close()

            data1 = [elem.split() for elem in lin1]

            # take header parameters
            header = {}
            for line in data1:
                if len(line) > 1 and line[0] in params['header']:
                    header[line[0]] = line[1:]

            # check type
            type_surv = compare_type(params, header)

            # check sample
            sample_surv = compare_sample(params, header)

            # store only if all requirements are satisfied
            if type_surv and sample_surv:
                final_list.append(name)

    op2 = open(output_file, 'w')
    for item in final_list:
        op2.write(item + '\n')
    op2.close()

    screen('Found ' + str(len(final_list)) +
           ' SN satisfying sample and type cuts.', params)
    screen('Surviving objects are listed in file ' + output_file, params)


def read_fitted(lc_data, mean_file):
    """
    Read GP results and populate dictionary parameters.

    input:  lc_data, LC object
            updated with user choices

            mean_file, str
            name of file storing previous calculated GP results.

    output: updated dictionary of parameters.
    """
    loaded = {}

    if bool(int(lc_data['n_samples'][0])):
        op1 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + \
                   lc_data['SNID:'][0] + '_samples.dat', 'r')
        lin1 = op1.readlines()
        op1.close()

        data1 = [elem.split() for elem in lin1]

        loaded['realizations'] = {}
        loaded['xarr'] = {}
        for fil in lc_data['filters']:
            par = lc_data['n_samples'][0]
            loaded['realizations'][fil] = [[float(data1[kk][jj])
                                            for kk in xrange(len(data1))
                                            if data1[kk][0] == fil]
                                            for jj in xrange(2, int(par) + 2)]

            loaded['xarr'][fil] = []
            for i in xrange(len(data1)):
                if data1[i][0] == fil:
                    loaded['xarr'][fil].append(float(data1[i][1]))

    op2 = open(mean_file, 'r')
    lin2 = op2.readlines()
    op2.close()

    data2 = [elem.split() for elem in lin2]

    loaded['GP_std'] = {}
    loaded['GP_fit'] = {}
    loaded['xarr'] = {}
    for fil in lc_data['filters']:
        loaded['xarr'][fil] = [float(data2[j][1])
                               for j in xrange(1, len(data2))
                               if data2[j][0] == fil]

        loaded['GP_fit'][fil] = [float(data2[j][2])
                                 for j in xrange(1, len(data2))
                                 if data2[j][0] == fil]

        loaded['GP_std'][fil] = [float(data2[j][3])
                                 for j in xrange(1, len(data2))
                                 if data2[j][0] == fil]

    return loaded


def main():
    """Print docstring."""
    print __doc__

if __name__ == '__main__':
    main()
