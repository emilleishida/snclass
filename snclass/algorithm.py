"""
Created by Emille Ishida in July 16th, 2015.

Collection of functions to implement supernova photometric classification.

- set_lc_list:
    Build a list of all objects satisfying selection cuts and plot them.
- build_sample:
    Build a directory holding all raw data passing selection cuts.
- sample_pop:
    Count number of each type in sample.
- photo_frac:
    Determine the fraction of each class in photo sample.
- get_names:
    Separate object identification according to class.
- set_parameters:
    Set extra sample parameters and copy raw files to new directory.
- select_GP:
    Select original objs to build a synthetic spectroscopic sample.
- build_spec_matrix:
    Build spectroscopic data matrix.
- plot_proj:
    Plot kPCA projections of training sample and test objects.
- read_matrix
    Read spectroscopic data matrix from file.
- read_hyperpar:
    Read hyperparameters result from cross-validation.
- set_kpca_obj:
    Set kpca object based on cross-validation results.
- classify_1obj:
    Perform classification of 1 supernova.
- classify:
    Classify all objects in photometric sample.
"""

def read_file(fname):
    """
    Read text file and return a list of all elements.

    input: fname, str
           file name

    output: data, list of list of str
            content of the file
    """
    op = open(fname, 'r')
    lin = op.readlines()
    op.close()

    data = [elem.split() for elem in lin]

    # remove empty lines
    while data.count([]) > 0:
        data.remove([])

    return data

def set_lclist(params):
    """
    Build a list of all objects satisfying selection cuts and plot them.

    input: params, dict
           keywords: plot_dir
                     path to store plots. If None do not build plots.
                     if None plots are not generated

                     fitted_data_dir
                     path to fitted data

                     list_dir
                     path to list directory 

                     sample
                     'spec' or 'photo'

                     user_choices, dict
                     output from snclass.read_user_input                
    """
    import numpy as np
    import pylab as plt
    import os

    from snclass.treat_lc import LC
    from snclass.util import translate_snid, read_snana_lc
    from snclass.functions import screen
    import sys

    # create plot directory
    if params['plot_dir'] is not None and \
    not os.path.isdir(params['plot_dir']):
        os.makedirs(params['plot_dir'])

    flist = os.listdir(params['fitted_data_dir'])

    photo_list = []
    cont = 0

    rfil = params['user_choices']['ref_filter'][0]

    for obj in flist:

        if 'mean' in obj and '~' not in obj and 'Y' not in obj:

            screen(obj, params['user_choices'])

            rname = translate_snid(obj)[0]
            params['user_choices']['path_to_lc'] = [rname]
            params['user_choices']['n_samples'] = ['0']

            raw = read_snana_lc(params['user_choices'])
            new_lc = LC(raw, params['user_choices'])

            if (params['user_choices']['file_root'][0] + raw['SNID:'][0] + \
               '_samples.dat' in flist):
                new_lc.user_choices['n_samples'] = ['100']
                new_lc.user_choices['samples_dir'] = [params['fitted_data_dir']]
                new_lc.load_fit_GP(params['fitted_data_dir'] + obj)

                l1 = [1  if len(new_lc.fitted['GP_fit'][fil]) > 0 else 0
                          for fil in params['user_choices']['filters']]

                if sum(l1) == len(params['user_choices']['filters']):
                    if rfil == 'None':
                        new_lc.normalize()
                    else:
                        new_lc.normalize(ref_filter=rfil)
                    new_lc.mjd_shift()
                    new_lc.check_epoch()

                    if new_lc.epoch_cuts:
                         photo_list.append(rname)

                         # only plot if not already done
                         if params['plot_dir'] is not None and \
                         not os.path.isfile(params['plot_dir'] + 'SN' + \
                                            raw['SNID:'][0] + '.png'):
                             new_lc.plot_fitted(file_out=\
                                                params['plot_dir'] + \
                                                'SN' + raw['SNID:'][0] + \
                                                '.png')
                    else:
                        screen('SN' + raw['SNID:'][0] + ' did not satisfy' + \
                               ' epoch cuts!\n', params['user_choices'])
                        cont = cont + 1
                else:
                    screen('SN' + raw['SNID:'][0] + ' does not exist in ' + \
                           'all filters!\n', params['user_choices'])
                    cont = cont + 1
            else:
                screen('Samples not found for SN' + raw['SNID:'][0], 
                       params['user_choices'])

        else:
            cont = cont + 1

    screen('Missed ' + str(cont) + ' SN.', params['user_choices'])

    # set parameter for file name
    if int(params['user_choices']['epoch_cut'][0]) < 0:
        epoch_min = str(abs(int(params['user_choices']['epoch_cut'][0])))
    else:
        epoch_min = 'p' + \
                    str(abs(int(params['user_choices']['epoch_cut'][0])))

    epoch_max = str(int(params['user_choices']['epoch_cut'][1]) - 1)

    filter_list = params['user_choices']['filters'][0]
    for item in params['user_choices']['filters'][1:]:
        filter_list = filter_list + item

    # save objs list
    if not os.path.isdir(params['list_dir']):
        os.makedirs(params['list_dir'])

    ref_filter = params['user_choices']['ref_filter'][0]
    if ref_filter is None:
        ref_fils = 'global'
    else:
        ref_fils = ref_filter

    op1 = open(params['list_dir'] + params['sample'] + '_' + filter_list + \
               '_' + epoch_min + '_' + epoch_max + '_ref_' + ref_fils + \
               '.list', 'w')
    for item in photo_list:
        op1.write(item + '\n')
    op1.close()

def build_sample(params):
    """
    Build a directory holding all raw data passing selection cuts.

    input: params, dict
           keywords:  'raw_dir' -> new directory to be created
                      'photo_dir' -> photometric LC fitted with GP
                      'spec_dir' -> sectroscopic LC fitted with GP
                      'user_choices' -> output from 
                                        snclass.util.read_user_input
    """
    import shutil
    from snclass.util import read_user_input, read_snana_lc, translate_snid
    from snclass.treat_lc import LC
    from snclass.functions import screen

    # create data directory
    if not os.path.isdir(params['raw_dir']):
        os.makedirs(params['raw_dir'])

    # read fitted light curves
    photo_list = os.listdir(params['photo_dir']) 
    spec_list = os.listdir(params['spec_dir'])

    # build filter list
    fil_list = params['user_choices']['filters'][0]
    for i in xrange(1, len(params['user_choices']['filters'])):
        fil_list = fil_list + params['user_choices']['filters'][i]

    for sn_set in [photo_list, spec_list]:
        for obj in sn_set:
            if 'samples' in obj and '~' not in obj and 'Y' not in obj:

                screen(obj, params['user_choices'])
                
                rname = translate_snid(obj)[0]
                params['user_choices']['path_to_lc'] = [rname]
                params['user_choices']['n_samples'] = ['0']

                # read raw data
                raw = read_snana_lc(params['user_choices'])
                new_lc = LC(raw, params['user_choices'])

                # load GP fit
                if sn_set == photo_list:
                    new_lc.load_fit_GP(photo_dir + 'DES_SN' + 
                                       raw['SNID:'][0] + '_mean.dat')
                else:
                    new_lc.load_fit_GP(spec_dir + 'DES_SN' + raw['SNID:'][0] +
                                       '_mean.dat')

                l1 = [1  if len(new_lc.fitted['GP_fit'][fil]) > 0  else 0 
                      for fil in params['user_choices']['filters']]

                if sum(l1) == len(params['user_choices']['filters']):
                    # treat light curve
                    new_lc.normalize(ref_filter= \
                                     params['user_choices']['ref_filter'][0])
                    new_lc.mjd_shift()
                    new_lc.check_basic()
                    new_lc.check_epoch()

                    # check epoch cuts
                    data_path = params['user_choices']['path_to_obs'][0]
                    if new_lc.epoch_cuts:
                        shutil.copy2(data_path + rname, raw_dir + rname)
                    else:
                        screen('... SN' + raw['SNID:'][0] + \
                               ' fail to pass epoch cuts!',
                               params['user_choices'])

def sample_pop(user_choices, params, type_number):
    """
    Count number of each type in sample.

    input: user_choices, dict
           output from snclass.util.read_user_choices

           params, dict
           keywords: 'list_name', str
                     name of file with list of all objs in this sample

           type_number, dict
           dictionary to translate types between raw data and final
           classification
           keywords -> final classificaton elements
           values -> identifiers in raw data

    output: sample_pop, dict
            keywords -> final types
            values -> number of objects of this type
            'tot' -> total number of objs in the sample
    """
    from snclass.util import read_snana_lc
    
    #get population for each SN type in spec sample
    fsample = read_file(params['list_name'])
    sample_pop = {}
    for name in fsample:
        user_choices['path_to_lc'] = [name[0]]
        raw = read_snana_lc(user_choices)
        for type_name in type_number.keys():
            if raw['SIM_NON1a:'][0] in type_number[type_name]:       
                if  type_name not in sample_pop.keys():
                    sample_pop[type_name] = 1
                else:
                    sample_pop[type_name] = sample_pop[type_name] + 1

    sample_pop['tot'] = sum([val for val in sample_pop.values()])

    return sample_pop


def photo_frac(spec_pop, photo_pop, representation):
    """
    Determine the fraction of each class in photo sample.

    input: spec_pop, dict
           output from sample_pop applied to spectroscopic (training) sample

           photo_pop, dict
           output from sample_pop applied to photometric (test) sample

           representation, str
           if 'original' the number of objects are not changed
           if 'balanced' final sample contain the same number of objs from 
              all types
           if 'representative' final spec sample resambles the proportions in
               photometric sample

    output: photo_frac, dict
            keywords -> final classes
            values -> fraction in photometric (test) sample
    """
    # get fraction of each type in photo sample
    photo_frac = {}
    for name in photo_pop.keys():
        if representation == 'original':
            photo_frac[name] = spec_pop[name]/float(spec_pop['tot'])
        elif representation == 'balanced':
            photo_frac[name] = 1.0/(len(spec_pop.keys()) - 1)
        elif representation == 'representative':
            photo_frac[name] = photo_pop[name]/float(photo_pop['tot'])

    return  photo_frac

def get_names(user_choices, params, type_number):
    """
    Separate object identification according to class.

    input: user_choices, dict
           output from snclass.util.read_user_input

           params, dict
           keywords: 'list_name', str
                      name of file with list of all objs in this sample

           type_number, dict
           dictionary to translate types between raw data and final
           classification
           keywords -> final classificaton elements
           values -> identifiers in raw data

    output: surv_spec_names, dict
            keywords -> final classes identification
            values -> set of objects ids for this class
    """
    from snclass.util import read_snana_lc

    # store name of objs surviving selection cuts
    surv_spec_names = {}
    fsample = read_file(params['list_name'])
    for name in fsample:
        user_choices['path_to_lc'] = [name[0]]

        try_lc = read_snana_lc(user_choices)

        stype = try_lc['SIM_NON1a:'][0]
        for type_name in type_number.keys():
            if stype in type_number[type_name]:   
                if type_name not in surv_spec_names.keys():   
                    surv_spec_names[type_name] = [name[0]]
                else:
                    surv_spec_names[type_name].append(name[0])

    return surv_spec_names

def set_parameters(params):
    """
    Set extra sample parameters and copy raw files to new directory.

    input: params, dict
           keywords: 'sample_size' -> number of objs in synthetic sample
                     'photo_perc' -> output from photo_frac
                     'spec_pop' -> output from sample_pop in spectroscopic
                                   (training) set
                     'synthetic_dir' -> directory to store synthetic sample
                     'list_name' -> file holding list of all objs in 
                                    spectroscopic sample

                     'fitted_data_dir' -> directory with spectroscopic
                                          sample fitted with GP
                     'representation'->
                         if 'original' the number of objects are not changed
                         if 'balanced' final sample contain the same number
                            of objs from all types
                         if 'representative' final spec sample resambles the
                            proportions in photometric sample

    output: update params dict
            new keywords: 'spec_num', dict ->  number of objs expected in the
                                               synthetic sample per class
                           'draw_spec_samples', dict - > number of draws
                                                         in each class
    """
    import numpy as np
    import os
    import shutil

    from snclass.util import translate_snid

    # calculate sample size
    if params['representation'] == 'original' or params['sample_size'] is None:
        params['sample_size'] = params['spec_pop']['tot']

    #construct number of SN expected in spec sample
    params['spec_num'] = {}
    for item in params['photo_perc'].keys():
        params['spec_num'][item] = np.round(params['sample_size'] * \
                                            params['photo_perc'][item])

    #define number of SN to draw from spec sample
    params['draw_spec_samples'] = {}
    for item in params['spec_pop'].keys():
        if item is not 'tot':
            diff = params['spec_num'][item] - params['spec_pop'][item]
            params['draw_spec_samples'][item] = int(diff) if diff > 0 else 0

    #construct synthetic spec data directory
    if not os.path.isdir(params['synthetic_dir']):
        os.makedirs(params['synthetic_dir'])

    # copy all mean files to synthetic directory
    fsample = read_file(params['list_name'])
    for fname in fsample:
        new_name = 'DES_SN' + translate_snid(fname[0]) + '_mean.dat'
        shutil.copy2(params['fitted_data_dir'] + new_name, 
                     params['synthetic_dir'] + new_name)

    return params


def select_GP(params, user_choices):
    """
    Select original objs to build a synthetic spectroscopic sample.

    input: params, dict
           output from set_paramameters

           user_choices, dict
           output from snclass.util.read_user_input
    """
    from snclass.util import translate_snid, read_snana_lc
    from snclass.functions import screen
    from snclass.treat_lc import LC
    from snclass.fit_lc_gptools import save_result

    import os
    import numpy as np
    import sys

    # set reference filter
    if user_choices['ref_filter'][0] == 'None':
        fil_choice = None
    else:
        fil_choice = user_choices['ref_filter'][0] 

    # select extra GP realizations in order to construct 
    # a representative spec sample
    for key in params['draw_spec_samples'].keys():
        cont = 0
        fail = 0

        # check if there are existing objs in this sample
        screen('... Check existing objs', user_choices)
        ready = []
        for obj in params['surv_spec_names'][key]:
            obj_id = translate_snid(obj)

            for j in xrange(params['draw_spec_samples'][key]):
                mean_file = params['synthetic_dir'] + '/' + \
                            user_choices['file_root'][0] + str(j) + \
                            'X' + obj_id + '_mean.dat'

                if os.path.isfile(mean_file) and mean_file not in ready:
                    cont = cont + 1
                    ready.append(mean_file)
                    screen('Found ready SN ' + str(cont) + 'X' + \
                           obj_id, user_choices)

        while cont < params['draw_spec_samples'][key]:
        
            # draw one of the objs in the spec sample
            indx = np.random.randint(0, params['spec_pop'][key])
            name = params['surv_spec_names'][key][indx]

            user_choices['path_to_lc'] = [name]

            # read light curve raw data
            raw = read_snana_lc(user_choices)

            if os.path.isfile(params['fitted_data_dir'] + 'DES_SN' + \
                              raw['SNID:'][0] + '_samples.dat'):

                # initiate light curve object
                my_lc = LC(raw, user_choices)

                screen('Loading SN' + raw['SNID:'][0], user_choices)

                # load GP fit
                my_lc.user_choices['n_samples'] = ['100']
                my_lc.user_choices['samples_dir'] = [params['fitted_data_dir']]
                my_lc.load_fit_GP(params['fitted_data_dir'] + 'DES_SN' + \
                                  raw['SNID:'][0] + '_mean.dat')
  

                l1 = [1  if len(my_lc.fitted['GP_fit'][fil]) > 0  else 0 
                      for fil in user_choices['filters']]
                if sum(l1) == len(user_choices['filters']):

                    # normalize
                    my_lc.normalize(samples=True, ref_filter=fil_choice)

                    # shift to peak mjd
                    my_lc.mjd_shift()

                    # check epoch requirements
                    my_lc.check_epoch()

                    if my_lc.epoch_cuts:

                        screen('... Passed epoch cuts', user_choices)
                        screen('... ... This is SN type ' +  raw['SIM_NON1a:'][0] + \
                               ' number ' + str(cont + 1) + ' of ' + 
                               str(params['draw_spec_samples'][key]), user_choices)

                        # draw one realization
                        size = len(my_lc.fitted['realizations'][user_choices['filters'][0]])
                        indx2 = np.random.randint(0, size)

                        for fil in user_choices['filters']:
                            print '... ... ... filter ' + fil

                            raw['GP_fit'][fil] = my_lc.fitted['realizations'][fil][indx2]
                            raw['GP_std'][fil] = my_lc.fitted['GP_std'][fil]
                            raw['xarr'][fil] = my_lc.fitted['xarr'][fil]
 
                        # set new file root
                        raw['file_root'] = [user_choices['file_root'][0] + \
                                             str(cont) + 'X']
                        raw['samples_dir'] = [params['synthetic_dir'] + '/']
                        save_result(raw)

                        # check epoch for this realization
                        new_lc = LC(raw, user_choices)
                        new_lc.load_fit_GP(params['synthetic_dir'] + '/' + \
                                       user_choices['file_root'][0] + str(cont) + \
                                       'X' + raw['SNID:'][0] + '_mean.dat')
                        new_lc.normalize(ref_filter=fil_choice)
                        new_lc.mjd_shift()
                        new_lc.check_epoch()

                        if new_lc.epoch_cuts:
                            cont = cont + 1
                        else:
                            screen('Samples failed to pass epoch cuts!\n', user_choices)
                            os.remove(params['synthetic_dir'] + '/' +
                                      user_choices['file_root'][0] + str(cont) + \
                                  'X' + raw['SNID:'][0] + '_mean.dat')
                        print '\n'

                    else:
                        screen('Failed to pass epoch cuts!\n', user_choices)
                        fail = fail + 1

                    if fail > 10 * params['spec_pop'][key]:
                        cont = 100000
                        sys.exit()

def build_spec_matrix(params, type_number):
    """
    Build spectroscopic data matrix.

    input: params, dict
           keywords: 'input_file' -> user input file
                     'synthetic_dir' -> directory with GP fitted spec sample
                     'out_dir' -> directory to store classification results
                     'mat_dir' -> directory to store data matrices
                     'range_pcs', list -> min and max number of PCs to be
                                          probed during cross-validation
                     'representation' -> type of modificationin 

           type_number, dict
           dictionary to translate types between raw data and final
           classification
           keywords -> final classificaton elements
           values -> identifiers in raw data
    """
    import os
    import numpy as np

    from snclass.matrix import DataMatrix
    from snclass.functions import screen, set_types

    # check matrices directory
    if not os.path.isdir(params['mat_dir']):
        os.makedirs(params['mat_dir'])

    # initiate data matrix obj 
    d = DataMatrix(params['input_file'])

    # build filter set
    fils = d.user_choices['filters'][0]
    for item in d.user_choices['filters'][1:]:
        fils = fils + item

    # build mjd boundaries
    if d.user_choices['epoch_cut'][0][0] == '-' or \
    d.user_choices['epoch_cut'][0] == '0':
        mjd_min = d.user_choices['epoch_cut'][0][1:]
    else:
        mjd_min = 'p' + d.user_choices['epoch_cut'][0]      

    mjd_max = str(int(d.user_choices['epoch_cut'][1]) - 1)

    # set reference filter flag
    if d.user_choices['ref_filter'][0] == None:
        fil_ref = 'global'
        d.user_choices['ref_filter'] = [None]
        fil_choice = None
    else:
        fil_ref = d.user_choices['ref_filter'][0]        
        fil_choice = fil_ref

    # store user defined number of samples
    nsamples =  d.user_choices['n_samples']

    # set keywords for spec sample
    d.user_choices['samples_dir'] = [params['synthetic_dir']]
    d.user_choices['n_samples'] = ['0']

    # build matrix with user defined parameters
    if fil_choice == 'None':
        fil_choice = None

    d.build(file_out=params['mat_dir'] + params['representation'] + '_' + \
            fils + '_' + mjd_min + '_' + mjd_max  + '_ref_' + fil_ref + \
            '_data_matrix.dat', ref_filter=fil_choice)

    screen('\n Spec sample contain ' + str(d.datam.shape[0]) + ' SNe.\n',
           d.user_choices)

    orig_types = []
    for item in d.sntype:
        for names in type_number.keys():
            if item in type_number[names]:
                orig_types.append(names)

    orig_types = np.array(orig_types)


    for npcs in xrange(params['range_pcs'][0], params['range_pcs'][1]):

        screen('npcs = ' + str(npcs), d.user_choices)

        # create directories if necessary
        hyperpar_file = params['out_dir'] + str(npcs) + 'PC/hyperpar_values.dat'

        if not os.path.isdir(params['out_dir'] + str(npcs) + 'PC/'):
            os.makedirs(params['out_dir']  + str(npcs) + 'PC/')

        # optimize hyperparameters
        d.user_choices['ncomp_lim'] = [str(npcs), str(npcs + 1)]
        d.sntype = set_types(d.sntype)
        d.cross_val()

        screen('Hyperparameter for ' + str(npcs) + ' PCs:', d.user_choices)
        screen('     gamma: ' + str(d.final['gamma']), d.user_choices)

        # update hyperparameter values
        d.final_configuration()

        # keep kpcs
        kpcs = d.transf_test.alphas_

        # save hyperparameter values
        pars = d.transf_test.get_params()

        # open file for hyperparameter value storage
        op1 = open(hyperpar_file, 'w')
        for par in pars.keys():
            op1.write(str(par) + '    ' + str(pars[par]) + '\n')
        op1.write('ncomp    ' + str(npcs))
        op1.write('\n\n\n')
        for k in xrange(kpcs.shape[1]):
            op1.write('alphas_' + str(k + 1) + '    ')
        op1.write('\n')

        for line in kpcs:
            for item in line:
                op1.write(str(item) + '    ')
            op1.write('\n')
        op1.close()

def plot_proj(spec_matrix, data_test, labels, new_lc, plot_dir, pcs,
              true_type):
    """
    Plot kPCA projections of training sample and test objects.

    input: spec_matrix, array
           low dimension spectroscopic matrix

           data_test, array
           low dimensional vector of test object

           labels, vec
           classes of objects in spectroscopic sample

           new_lc, LC object
           processed light curve of test object

           plot_dir, str
           directory to store plots

           pcs, list
           PCs to be plotted

           true_type, str 
           type label of test object
    """
    # create type flags
    snIa = labels == 'Ia'
    snIbc = labels == 'Ibc'
    snII = labels == 'II'

    # plot
    plt.figure()
    plt.title('prob_Ia = ' + str(round(new_lc.prob_Ia, 2)))
    plt.scatter(spec_matrix[snII, pcs[0]], spec_matrix[snII, pcs[1]], color='purple', marker='s', label='II')
    plt.scatter(spec_matrix[snIbc, pcs[0]], spec_matrix[snIbc, pcs[1]], color='green', marker='^', label='Ibc')
    plt.scatter(spec_matrix[snIa, pcs[0]], spec_matrix[snIa, pcs[1]], color='blue', marker='o', label='Ia')
    plt.scatter(data_test[:,0], data_test[:,1], color='red', marker='*', label='test - ' + true_type)
    plt.xlabel('kPC ' + str(pcs[0] + 1), fontsize=15)
    plt.ylabel('kPC ' + str(pcs[1] + 1), fontsize=15)
    plt.legend()
    plt.savefig(plot_dir +  '/proj_SN' + new_lc.raw['SNID:'][0] + '.png')
    plt.close()

def read_matrix(data_matrix, convert_types=True):
    """
    Read spectroscopic sample matrix.

    input: data_matrix, str
           file holding spectroscopic matrix

           convert_types, bool
           if True use set_types function to convert labels to binary
           default is True

    output: data, array
            spectroscopic data matrix

            sntype, list
            types of each objs in spec sample (order is important!)

            binary_types, list
            types converted to binary format
            if convert_types is False this is returned as None
    """
    from snclass.functions import set_types

    import numpy as np

    # read data matrix and obj classes
    op1 = open(data_matrix, 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [elem.split() for elem in lin1[1:]]

    datam = []
    sntype = []
    for line in data1:
        obj = []
        sntype.append(line[1])
        for item in line[3:]:
            obj.append(float(item))
        datam.append(obj)

    data = np.array(datam)

    if convert_types:
        binary_types = set_types(sntype)
    else:
        binary_types = None

    return data, sntype, binary_types

def read_hyperpar(parameters):
    """
    Read hyperparameters result from cross-validation.

    input: cv_file, str
           file holding hyperparameter results from cross-validation

    output: pars, dict
            dictionary of hyperparameter values

            alphas, array
    """
    import numpy as np

    # determine types
    parameters['strings'] = ['eigen_solver', 'kernel']
    parameters['floats'] = ['alpha', 'gamma', 'n_components']
    parameters['ints'] = ['coef0', 'degree', 'ncomp', 'tol']
    parameters['bools'] = ['fit_inverse_transform', 'remove_zero_eig']
    parameters['nones'] = ['kernel_params', 'max_iter']

    # read hyperparameters values
    op2 = open(parameters['cv_file'], 'r')
    lin2 = op2.readlines()
    op2.close()

    data2 = [elem.split() for elem in lin2]
    pars = {}
    for line in data2:
        if len(line) > 0:
            if line[0] in parameters['floats']:
                pars[line[0]] = float(line[1])
            elif line[0] in parameters['ints']:
                pars[line[0]] = int(line[1])
            elif line[0] in parameters['bools'] and line[1] == 'False':
                pars[line[0]] = False
            elif line[0] in parameters['nones']:
                pars[line[0]] = None
            else:
                pars[line[0]] = line[1]
        else:
            indx = data2.index(line) + 3
            break

    alphas = []
    for line in data2[indx:]:
        l1 = []
        for item in line:
            l1.append(float(item))
        alphas.append(l1)

    alphas = np.array(alphas)

    return pars, alphas

def set_kpca_obj(pars, data, sntype, type_number):
    """
    Set kpca object based on cross-validation results.

    input: pars, dict
           output from read_hyperpar

           data, array
           output from read_matrix

           sntype, list
           output from read_matrix

           type_number, dict
           dictionary to translate types between raw data and final
           classification
           keywords -> final classificaton elements
           values -> identifiers in raw data

    output: obj_kpca, KernelPCA obj
            tailored with cross-validation results

            spec_matrix, array
            low dimension spectroscopic matrix

            labels, list
            classes as defined in raw data files
    """
    from sklearn.decomposition import KernelPCA
    import numpy as np

    # start kpca object
    obj_kpca = KernelPCA()
    obj_kpca.eigen_solver = pars['eigen_solver']
    obj_kpca.kernel = pars['kernel']
    obj_kpca.alpha = pars['alpha']
    obj_kpca.gamma = pars['gamma']
    obj_kpca.n_components = pars['n_components']
    obj_kpca.coef0 = pars['coef0']
    obj_kpca.degree = pars['degree']
    obj_kpca.tol = pars['tol']
    obj_kpca.fit_inverse_transform = pars['fit_inverse_transform']
    obj_kpca.remove_zero_eig = pars['remove_zero_eig']
    obj_kpca.kernel_params = pars['kernel_params']
    obj_kpca.max_iter = pars['max_iter']
    spec_matrix = obj_kpca.fit_transform(data)

    # construct label vector
    labels = []
    for elem in sntype:
        for classes in type_number.keys():
            if elem in type_number[classes]:
                labels.append(classes)

    labels = np.array(labels)

    return obj_kpca, spec_matrix, labels

def classify_1obj(din):
    """
    Perform classification of 1 supernova.

    input: din, dict - keywords, value type: 
                     user_input, dict -> output from read_user_input
                     name, str -> name of raw light curve file
                     type_number, dict -> translate between str and numerical
                                          classes identification
                     do_plot, bool -> if True produce plots, default is False

                     p1, dict ->  keywords, value type:
                         fname_photo_list, str: list of all photometric 
                                                sample objects
                         photo_dir, str: directory of GP fitted results
                                         for photo sample
                         range_pcs, list: [min_number_PCs, max_number_PCs]
                                          to be tested through cross-validation
                         SNR_dir, str: directory to store all results from 
                                       this SNR cut
                         out_dir, str: directory to store classification 
                                       results
                         plot_proj_dir, str: directory to store 
                                             projection plots
                         data_matrix, str: file holding spec data matrix

    output: class_results:
               list -> [snid, true_type, prob_Ia] 
    """
    from snclass.functions import screen, nneighbor
    from snclass.util import translate_snid, read_snana_lc
    from snclass.treat_lc import LC

    # update supernova name    
    din['user_input']['path_to_lc'] = [translate_snid(din['name'])[0]]

    # read raw data
    raw = read_snana_lc(din['user_input'])

    # set true type
    for names in din['type_number'].keys():
        if raw['SIM_NON1a:'][0] in din['type_number'][names]:
            true_type = names

    # load GP fit and test epoch cuts
    new_lc = LC(raw, din['user_input'])
    new_lc.user_choices['samples_dir'] = [din['p1']['photo_dir']]
    new_lc.load_fit_GP(din['p1']['photo_dir'] + din['name'])

    l1 = [1  if len(new_lc.fitted['GP_fit'][fil]) > 0  else 0 
          for fil in din['user_input']['filters']]

    fil_choice = din['user_input']['ref_filter'][0]
    if fil_choice == 'None':
        fil_choice = None

    if sum(l1) == len(din['user_input']['filters']):
        new_lc.normalize(samples=True, 
                         ref_filter=fil_choice)
        new_lc.mjd_shift()
        new_lc.check_epoch()

        if new_lc.epoch_cuts:

            screen(new_lc.raw['SNID:'][0], din['user_input'])

            # build matrix lines
            new_lc.build_steps(samples=True)

            # transform samples
            small_matrix = new_lc.samples_for_matrix
            data_test = din['p1']['obj_kpca'].transform(small_matrix)

            #classify samples
            new_label = nneighbor(data_test, din['p1']['spec_matrix'],
                                  din['p1']['binary_types'], din['user_input'])

            # calculate final probability
            ntypes = [1 for item in new_label if item == '0']
            new_lc.prob_Ia = sum(ntypes) / \
                             float(din['user_input']['n_samples'][0])
            
            if din['do_plot']:
                plot_proj(din['p1']['spec_matrix'], data_test, din['p1']['labels'],
                          new_lc, din['p1']['plot_dir'], [0,1], true_type)
            

            # print result to screen
            screen('SN' + new_lc.raw['SNID:'][0] + \
                   ',   True type: ' + true_type + ', prob_Ia = ' + \
                    str(new_lc.prob_Ia), din['user_input'])
 
            class_results = [new_lc.raw['SNID:'][0], true_type,
                             new_lc.prob_Ia]
            return class_results


def classify(p1, user_input, type_number, do_plot=False):
    """
    Classify all objects in photometric sample.

    input: p1, dict
           keywords, value type:
               fname_photo_list, str: list of all photometric sample objects
               photo_dir, str: directory of GP fitted results for photo sample
               range_pcs, list: [min_number_PCs, max_number_PCs] to be tested
                                through cross-validation
               SNR_dir, str: directory to store all results from this SNR cut
               out_dir, str: directory to store classification results
               plot_proj_dir, str: directory to store projection plots
               data_matrix, str: file holding spec data matrix

           user_input, dict
           output from snclass.util.read_user_input

           type_number, dict
           dictionary to translate types between raw data and final
           classification
           keywords -> final classificaton elements
           values -> identifiers in raw data

           do_plot, bool - optional
           if True, creature projection plot for all test objects
           default is False
    """
    from snclass.functions import screen
    from snclass.util import translate_snid, read_snana_lc
    from snclass.treat_lc import LC

    import os
    import sys
    import numpy as np
    from multiprocessing import Pool

    # read photometric sample
    photo_fname = read_file(p1['fname_photo_list'])
    photo_list = ['DES_SN' + translate_snid(item[0]) + '_mean.dat'
                  for item in photo_fname 
                  if os.path.isfile(p1['photo_dir'] + 
                                    'DES_SN' + translate_snid(item[0]) + 
                                    '_samples.dat') and '~' not in item[0]]

    for npcs in xrange(p1['range_pcs'][0], p1['range_pcs'][1]): 

        if int(user_input['epoch_cut'][0]) < 0:
            mjd_min = str(abs(int(user_input['epoch_cut'][0])))
        else:
            mjd_min = 'p' + str(abs(int(user_input['epoch_cut'][0])))

        if int(user_input['epoch_cut'][1]) < 0:
            mjd_max = 'm' + str(abs(int(user_input['epoch_cut'][1]) - 1))
        else:
            mjd_max = str(int(user_input['epoch_cut'][1]) - 1)

        fils = user_input['filters'][0]
        for item in user_input['filters'][1:]:
            fils = fils + item

        out_dir = p1['out_dir']
        plot_proj_dir = p1['plot_proj_dir']
        p1['out_dir'] = p1['out_dir'] + str(npcs) + 'PC/'
        p1['cv_file'] = p1['out_dir'] + 'hyperpar_values.dat'

        if plot_proj_dir is not None:
            p1['plot_proj_dir'] = p1['plot_proj_dir'] + str(npcs) + 'PC/'
            if not os.path.isdir(p1['plot_proj_dir']):
                os.makedirs(p1['plot_proj_dir'])

        if not os.path.isdir(p1['out_dir']):
            os.makedirs(p1['out_dir'])

        p1['pars'], p1['alphas'] = read_hyperpar(p1)
        p1['data'], p1['sntype'], p1['binary_types'] = \
            read_matrix(p1['data_matrix'])
        p1['obj_kpca'], p1['spec_matrix'], p1['labels'] = \
            set_kpca_obj(p1['pars'], p1['data'], p1['sntype'], type_number)

        pars = []
        for name in photo_list:
            ptemp = {}
            ptemp['p1'] = p1
            ptemp['name'] = name
            ptemp['type_number'] = type_number
            ptemp['user_input'] = user_input
            ptemp['do_plot'] = do_plot

            pars.append(ptemp)

        pool = Pool(processes=int(user_input['n_proc'][0]))
        my_pool = pool.map_async(classify_1obj, pars)
        try:
            results = my_pool.get(0xFFFF)
        except KeyboardInterrupt:
            print 'Interruputed by the user!'
            sys.exit()

        pool.close()
        pool.join() 

        p1['out_dir'] = out_dir
        p1['plot_proj_dir'] = plot_proj_dir

        op2 = open(p1['out_dir'] +  'class_res_' + str(npcs) + 'PC.dat', 'w')
        op2.write('SNID    true_type    prob_Ia\n')
        for line in results:
            for item in line:
                op2.write(item + '    ')
            op2.write('\n')
        op2.close()

