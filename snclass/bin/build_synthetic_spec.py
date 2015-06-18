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
This file contains tools to build a synthetic training sample based on the
proportions between supernova types as described in Ishida et al., 2015.

Usage:

$ build_synthetic_spec.py -i <user_input_file> -d <synthetic_sample_directory>

Functions:
- check_pop:
    Check the number of objects from each time in a given sample.
- calc_fraction:
    Calculate the fraction of each type in a given sample.
- check_fitted:
    Check object types on previously fitted sample.
- check_mean_GP:
    Check if there is a GP fit file.
- test_epoch:
    Check if previously calculated GP fit survives selection cuts.
- setup_gp:
    Set up necessary GP object.    
"""

import os
import asciitable
import numpy as np
import argparse
import shutil

from snclass.matrix import DataMatrix
from snclass.functions import screen
from snclass.util import translate_snid, read_user_input, choose_sn
from snclass.util import read_snana_lc
from snclass.treat_lc import LC
from snclass.fit_lc_gptools import imp_gptools, save_result


def check_pop(sample_list, user_choices):
    """
    Check the number of objects from each time in a given sample.

    input: sample_list, str
           list of objects raw data files.

           user_choices, dict
           output from read_user_input

    output: pop, dict
            keys are types and values are the number of objects in the sample
    """
    # read raw data files names 
    flist = asciitable.read(sample_list)

    # count types
    pop = {}
    for name in flist:
        user_choices['path_to_lc'] = [name[0]]
        raw = read_snana_lc(user_choices)
        if raw['SIM_NON1a:'][0] not in pop.keys():
            pop[raw['SIM_NON1a:'][0]] = 1
        else:
            pop[raw['SIM_NON1a:'][0]] = pop[raw['SIM_NON1a:'][0]] + 1

    return pop


def calc_fraction(pop):
    """
    Calculate the fraction of each type in a given sample.

    input: pop, dict
           keys are types and values are the number of objects in the sample

    output: pop_frac, dict
            keys are types and values are corresponding fraction in the entire
            sample
    """
    # get fraction of each type in photo sample
    frac = {}
    tot = sum([val for val in pop.values()])
    for name in pop.keys():
        frac[name] = pop[name]/float(tot)

    return frac


def check_fitted(name_dir, user_choices):
    """
    Check object types on previously fitted sample.

    input: name_dir, str
           directory where fitted samples are stored.

           user_choices, dict
           output from snclass.util.read_user_input

    output: surv_names, dict
            keys are types and values are names of raw data files
    """
    # list classes surviving selection cuts
    surv = os.listdir(name_dir)

    surv_names = {}
    for name in surv:
        user_choices['path_to_lc'] = [translate_snid(name)[0]]

        try_lc = read_snana_lc(user_choices)
        stype = try_lc['SIM_NON1a:'][0]

        if try_lc['SIM_NON1a:'][0] not in surv_names.keys():
            surv_names[stype] = user_choices['path_to_lc']
        else:
            surv_names[stype].append(user_choices['path_to_lc'][0])

    return surv_names


def check_mean_GP(key, surv, draw, user_choices, synthetic_dir):
    """
    Check if there is a GP fit file.

    input: key, str
           object type

           surv, dict
           dictionary of objects surviving basic cuts
           keys are types, values are GP fit mean file name

           draw, dict
           dictionary of number of extra draws necessary for each type
           keys are types, values are number of draws

           user_choices, dict
           output from snclass.util.read_user_input

           synthetic_dir, str
           directory where synthetic sample is stored

    output: ready, list
            list of already existing files
    """
    ready = []
    
    # run through all types
    if key in surv.keys():
        for obj in surv[key]:

            # get object id
            obj_id = translate_snid(obj)

            # run through all types which shall be re-sampled,
            # check if they already exists
            for j in xrange(draw[key]):
                mean_file = synthetic_dir + '/' + \
                            user_choices['file_root'][0] + str(j) + \
                            'X' + obj_id + '_mean.dat'
                if os.path.isfile(mean_file) and mean_file not in ready:
                    ready.append(mean_file)
                    screen('Found ready SN ' + str(len(ready)) + 'X' + obj_id, user_choices)
    else:
        screen('type ' + str(key) + ' not present in surviving sample.')
    
    return ready


def test_epoch(key, surv, user_choices):
    """
    Check if previously calculated GP fit survives selection cuts.
    
    input: key, str
           object type

           surv, dict
           dictionary of objects surviving basic cuts
           keys are types, values are GP fit mean file name

           user_choices, dict
           output from snclass.util.read_user_input

    output: my_lc, LC object
            updated light curve object after checked for epoch cuts            
    """
    # sample a random obj in the training sample
    indx = np.random.randint(0, len(surv[key]))
    name = surv[key][indx]

    # determine fitting method
    fit_method = bool(int(user_choices['do_mcmc'][0]))

    # update path to raw data
    user_choices['path_to_lc'] = [name]

    # read light curve raw data
    raw = read_snana_lc(user_choices)

    # update raw data with user choices
    raw.update(user_choices)

    # set number of samples to 0 (we are only interested in the mean for now)
    raw['n_samples'] = ['0']

    # initiate light curve object
    my_lc = LC(raw, user_choices)

    screen('Fitting SN' + raw['SNID:'][0], user_choices)

    # load GP fit
    my_lc.load_fit_GP(user_choices['samples_dir'][0] + '/DES_SN' + raw['SNID:'][0] + '_mean.dat')

    # normalize
    my_lc.normalize()

    # shift to peak mjd
    my_lc.mjd_shift()

    # check epoch requirements
    my_lc.check_epoch()

    return my_lc, raw


def setup_gp(raw, user_choices, gp_objs):
    """
    Set up necessary GP object. 

    input: raw, dict
           second output from test_epoch.

           user_choices, dict
           output from snclass.util.read_user_input

           gp_objs, dict
           already calculated GP objects

    output: raw, dict
            update light curve data

            gp_objs, dict
            update set of GP objects
    """
    # list of necessary keywords
    key_list = ['realizations', 'xarr', 'GP_std', 'GP_obj']
    fit_method = bool(int(user_choices['do_mcmc'][0]))

    for name in key_list:
        if name not in raw.keys():
            raw[name] = {}

    for fil in user_choices['filters']:
        screen('... ... filter ' + fil, user_choices)

        if raw['SNID:'][0] not in gp_objs.keys():
            raw = imp_gptools(raw, fil, mcmc=fit_method)
            new_obj = raw['GP_obj'][fil]
            draws = new_obj.draw_sample(raw['xarr'][fil],
                                        num_samp=1)
            raw['GP_fit'][fil] = [line[0] for line in draws]

        else:
            new_obj = gp_objs[raw['SNID:'][0]][fil]
            draws = new_obj.draw_sample(raw['xarr'][fil],
                                                num_samp=1)
            raw['GP_fit'][fil] = [line[0] for line in draws]

    if raw['SNID:'][0] not in gp_objs.keys():
        gp_objs[raw['SNID:'][0]] = raw['GP_obj']
        
    return raw, gp_objs


def main(args):
    """
    Construct 'fake' training.

    Use photometric simulated sample to guess proportions between
    spectroscopic sample classes.
    """
    # read user input
    user_choices = read_user_input(args.input)

    ##########################################################################
    # Spec

    # build complete spec list
    screen('Building spectroscopic sample.', user_choices)
    user_choices['sample_cut'] = ['1', '3', '21', '22', '23', '32', '33']
    spec_list = choose_sn(user_choices, output_file='spec.list')

    # check population according to type
    spec_pop = check_pop('spec.list', user_choices)

    # count spec classes surviving selection cuts
    surv_spec = check_fitted(user_choices['samples_dir'][0], user_choices)

    ##########################################################################
    # Photo

    #build complete photo list
    screen('Build photometric samples.', user_choices)
    user_choices['sample_cut'] = ['-9']
    photo_list = choose_sn(user_choices, output_file='photo.list')

    # check population according to type
    photo_pop = check_pop('photo.list', user_choices)
    photo_frac = calc_fraction(photo_pop)

    ##########################################################################
    # Building fake training sample

    screen('Checking compatibility.', user_choices)

    # construct number of SN expected in spec sample
    spec_num = {}
    for item in photo_pop.keys():
        spec_num[item] = int(np.round(sum(spec_pop.values()) * photo_frac[item]))

    #construct synthetic spec data directory
    synthetic_dir = args.dir
    if not os.path.isdir(synthetic_dir):
        os.makedirs(synthetic_dir)

    # collect gp objects
    gp_objs = {}

    #run through all types
    for key in spec_num.keys():

        # start cont of failed tries
        fail = 0

        if key in surv_spec.keys():

            # check which objs and samples were already calculated
            ready = check_mean_GP(key, surv_spec, spec_num, user_choices, synthetic_dir)

            cont = len(ready)

            while cont < spec_num[key]:
  
                my_lc = test_epoch(key, surv_spec, user_choices)

                mean_name = synthetic_dir + '/' + \
                            user_choices['file_root'][0] + my_lc[1]['SNID:'][0] + \
                           '_mean.dat'

                screen('... This is SN type ' +  my_lc[1]['SIM_NON1a:'][0] + \
                       ' number ' + str(len(ready) + 1) + ' of ' + 
                       str(spec_num[key]), user_choices) 

                if my_lc[0].epoch_cuts and mean_name not in ready:
                    ready.append(mean_name)
                    cont = len(ready)

                    shutil.copy2(user_choices['samples_dir'][0] + '/' + \
                                 user_choices['file_root'][0] + \
                                 my_lc[1]['SNID:'][0] + '_mean.dat',
                                 mean_name)
                    screen('\n', user_choices)

                else:                

                    # build GP object
                    raw, gp_objs = setup_gp(my_lc[1], user_choices, gp_objs)

                    # set new file root
                    raw['file_root'] = [user_choices['file_root'][0] + \
                                        str(len(ready)) + 'X']
                    raw['samples_dir'] = [synthetic_dir + '/']
                    save_result(raw)

                    # check epoch for this realization
                    new_lc = LC(raw, user_choices)
                    new_lc.load_fit_GP(synthetic_dir + '/' + \
                                       user_choices['file_root'][0] + str(cont) + \
                                       'X' + raw['SNID:'][0] + '_mean.dat')
                    new_lc.normalize()
                    new_lc.mjd_shift()
                    new_lc.check_epoch()

                    if new_lc.epoch_cuts:
                        ready.append(synthetic_dir + '/' +
                                     user_choices['file_root'][0] + str(cont) + \
                                     'X' + raw['SNID:'][0] + '_mean.dat')
                        cont = len(ready)
                        screen('\n', user_choices)

                    else:
                        os.remove(synthetic_dir + '/' + \
                                  user_choices['file_root'][0] + str(cont) + \
                                  'X' + raw['SNID:'][0] + '_mean.dat')
                        fail = fail + 1
                        screen(str(fail) + ' samples failed to pass epoch cuts!\n', user_choices)

                        if fail > 10 * spec_num[key]:
                            cont = 100000


if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Build fake training sample.')
    parser.add_argument('-i','--input', dest='input', help='Input file name', 
                        required = True)
    parser.add_argument('-d', '--dir', dest='dir', 
                        help='Directory with fitted light curves.', required=True) 

    args = parser.parse_args()
   
    main(args)

