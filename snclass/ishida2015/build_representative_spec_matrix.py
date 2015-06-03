from snclass.matrix import DataMatrix
from snclass.functions import screen
import snclass
import os
import asciitable
import numpy as np
import shutil

#build complete spec list
user_choices = snclass.util.read_user_input("fit_lc_input.dat")
snclass.util.choose_sn(user_choices, output_file='spec.list')

#get population for each SN type in spec sample
fspec = asciitable.read('spec.list')
spec_pop = {}
for name in fspec:
    user_choices['path_to_lc'] = [name[0]]
    raw = snclass.util.read_snana_lc(user_choices)
    if raw['SIM_NON1a:'][0] not in spec_pop.keys():
        spec_pop[raw['SIM_NON1a:'][0]] = 1
    else:
       spec_pop[raw['SIM_NON1a:'][0]] = spec_pop[raw['SIM_NON1a:'][0]] + 1

#build complete photo list
user_choices['sample_cut'] = ['-9']
snclass.util.choose_sn(user_choices, output_file='photo.list')

#get population for each SN type in photo sample
fphoto = asciitable.read('photo.list')
photo_pop = {}
for name in fphoto:
    user_choices['path_to_lc'] = [name[0]]
    raw = snclass.util.read_snana_lc(user_choices)
    if raw['SIM_NON1a:'][0] not in photo_pop.keys():
        photo_pop[raw['SIM_NON1a:'][0]] = 1
    else:
       photo_pop[raw['SIM_NON1a:'][0]] = photo_pop[raw['SIM_NON1a:'][0]] + 1

#get fraction of each type in photo sample
photo_frac = {}
photo_tot = sum([val for val in photo_pop.values()])
for name in photo_pop.keys():
    photo_frac[name] = photo_pop[name]/float(photo_tot)

#count spec classes surviving selection cuts
surv_spec = os.listdir('spec_SNR5/')

surv_spec_pop = {}
surv_spec_names = {}
for name in surv_spec:
    user_choices['path_to_lc'] = [translate_snid(name)[0]]

    try_lc = snclass.util.read_snana_lc(user_choices)

    stype = try_lc['SIM_NON1a:'][0]
    if try_lc['SIM_NON1a:'][0] not in surv_spec_pop.keys():
        surv_spec_pop[stype] = 1
        surv_spec_names[stype] = user_choices['path_to_lc']
    else:
        surv_spec_pop[stype] = surv_spec_pop[stype] + 1
        surv_spec_names[stype].append(user_choices['path_to_lc'][0])

#construct number of SN expected in spec sample
spec_num = {}
for item in photo_pop.keys():
    spec_num[item] = np.round(sum(spec_pop.values()) * photo_frac[item])

#define number of SN to draw from surviving spec sample
draw_spec_samples = {}
for item in surv_spec_pop.keys():
    diff = spec_num[item] - surv_spec_pop[item]
    draw_spec_samples[item] = int(diff) if diff > 0 else 0

#construct synthetic spec data directory
synthetic_dir = 'fake_spec_SNR_v3/'
if not os.path.isdir(synthetic_dir):
    os.makedirs(synthetic_dir)

#collect gp objects
gp_objs = {}

for key in draw_spec_samples.keys():
    cont = 0
    fail = 0
    while cont < draw_spec_samples[key]:

        indx = np.random.randint(0, surv_spec_pop[key])
        name = surv_spec_names[key][indx]

        fit_method = bool(int(user_choices['do_mcmc'][0]))

        user_choices['path_to_lc'] = [name]

        # read light curve raw data
        raw = snclass.util.read_snana_lc(user_choices)

        # update raw data with user choices
        raw.update(user_choices)
        raw['n_samples'] = ['0']

        # initiate light curve object
        my_lc = snclass.treat_lc.LC(raw, user_choices)

        screen('Fitting SN' + raw['SNID:'][0], user_choices)

        # load GP fit
        my_lc.load_fit_GP('spec_SNR5/DES_SN' + raw['SNID:'][0] + '_mean.dat')

        # normalize
        my_lc.normalize()

        # shift to peak mjd
        my_lc.mjd_shift()

        # check epoch requirements
        my_lc.check_epoch()

        if my_lc.epoch_cuts:

            shutil.copy2('spec_SNR5/DES_SN' + raw['SNID:'][0] + '_mean.dat',
                         synthetic_dir)

            screen('... Passed epoch cuts', user_choices)
            screen('... ... This is SN type ' +  raw['SIM_NON1a:'][0] + \
                  ' number ' + str(cont + 1) + ' of ' + \
                  str(draw_spec_samples[key], user_choices)

            key_list = ['realizations', 'xarr', 'GP_std', 'GP_obj']

            for name in key_list:
                if name not in raw.keys():
                    raw[name] = {}

            for fil in user_choices['filters']:
                # fit
                print '... ... ... filter ' + fil

                if raw['SNID:'][0] not in gp_objs.keys():
                    raw = snclass.fit_lc_gptools.imp_gptools(raw, fil,
                                                             mcmc=fit_method)
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

            # set new file root
            raw['file_root'] = [user_choices['file_root'][0] + \
                                str(cont) + 'X']
            raw['samples_dir'] = [synthetic_dir + '/']
            snclass.fit_lc_gptools.save_result(raw)

            # check epoch for this realization
            new_lc = snclass.treat_lc.LC(raw, user_choices)
            new_lc.load_fit_GP(synthetic_dir + '/' + \
                               user_choices['file_root'][0] + str(cont) + \
                               'X' + raw['SNID:'][0] + '_mean.dat')
            new_lc.normalize()
            new_lc.mjd_shift()
            new_lc.check_epoch()

            if new_lc.epoch_cuts:
                cont = cont + 1
            else:
                screen('Samples failed to pass epoch cuts!\n', user_choices)
                os.remove(synthetic_dir + '/' +
                          user_choices['file_root'][0] + str(cont) + \
                          'X' + raw['SNID:'][0] + '_mean.dat')
            print '\n'

        else:
            screen('Failed to pass epoch cuts!\n', user_choices)
            fail = fail + 1

        if fail > 10 * surv_spec_pop[key]:
            cont = 100000
