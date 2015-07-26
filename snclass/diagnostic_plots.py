"""
Created by Emille Ishida in July 17th, 2015.

Collection of functions to generate diagnostic plots for 
supernova classification.

- calc_ROC:
    Calculate and plot ROC for classification with different number of PCs.
- det_threshold:
    Write to file and on screen number of pcs, gamma and threshold results.
- calc_binned_diag:
    Calculate efficiency, purity and figure of merit per redshift bin.
- calc_global_diag:
    Calculate global efficiency, purity and figure of merit.
- count_pop:
    Count original population from spec and photo samples.
- read_literature_results:
    Read results from the literature build on post-SNPCC data.
- plot_diagnostics_aftercuts:
    Build diagnostic plots considering only objects surviving selection cuts.
- plot_diagnostics_before_cuts:
    Build diagnostic plots considering all objects in original data set.
"""

def calc_ROC(params, plot=False):
    """
    Calculate and plot ROC for classification with different number of PCs.

    input: params, dict
           keywords:
               range_pcs, list: [min_number_pcs, max_number_pcs]
               class_res_dir, str: directory to store classification results
               representation, str: 'original' or 'balanced' or 
                                    'representative'
               xmin, str: min mjd for plot title purpouses
               xmax, str: max mjd for plot title purpouses
               fils, str: concatenation of filters for plot title purpouses
               ref_filter, str: reference filter for plot title purpouses
               
            plot, bool
            if True, show on screen ROC plot.
            Default is False

    return: params, dict
            additional keywords:
                tpr_set, array: set of true positive rate results
                fpr_set, array: set of false positive rate results
                threshold_set, array: set of threshold values
                auc_set, array: set of AUC results
                min_npcs, int: minimum number of PCs leading to higher AUC
    """
    from snclass.algorithm import read_file
    from sklearn.metrics import roc_curve, auc

    import pylab as plt

    fpr_set = []
    tpr_set = []
    threshold_set = []
    auc_set = []

    for npcs in xrange(params['range_pcs'][0], params['range_pcs'][1]):

        # read classification results
        data1 = read_file(params['class_res_dir'] + str(npcs) + \
                          'PC/class_res_' + str(npcs) + 'PC.dat')

        types = []
        prob = []
        for line in data1[1:]:
            types.append(1 if line[1] == 'Ia' else 0)
            prob.append(float(line[2]))

        # compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(types, prob)
        roc_auc = auc(fpr, tpr)

        fpr_set.append(fpr)
        tpr_set.append(tpr)
        threshold_set.append(thresholds)
        auc_set.append(roc_auc)

    params['tpr_set'] = tpr_set
    params['fpr_set'] = fpr_set
    params['threshold_set'] = threshold_set
    params['auc_set'] = auc_set

    if plot:
        if params['ref_filter'] == None:
            ref_filter = 'global'
        else:
            ref_filter = params['ref_filter']

        # plot
        plt.clf()
        for i in xrange(len(fpr_set)):
            plt.plot(fpr_set[i], tpr_set[i], label=str(i + 2) + \
                     'PC (auc = %0.3f)' % auc_set[i])
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.0])
        plt.xlabel('False Positive Rate', fontsize=14)
        plt.ylabel('True Positive Rate', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(loc="lower right", title=params['representation'] + \
                   ' spec, -' + params['epoch_min'] + ' to +' + \
                   params['epoch_max'] + ' days, ' + params['filters'] + \
                   ',\n SNR' + params['SNR_cut'] + ', ' + \
                   ref_filter + ' max', fontsize=14)
        plt.show()

    # determine max auc
    max_auc = max(auc_set)

    # determine number of pcs
    n_pcs_max_auc = [i for i in xrange(len(auc_set)) if auc_set[i] == max_auc]
    min_npcs = min(n_pcs_max_auc) + 2

    params['final_npcs'] = min_npcs 

    return params

def det_threshold(params):
    """
    Write to file and on screen number of pcs, gamma and threshold results.

    input: params, dict
           output from calc_ROC
           + screen, bool: if True show results on screen
    """
    from snclass.algorithm import read_file
    from snclass.functions import screen

    # calculate distance for best auc curve
    dist = [params['tpr_set'][3][i] - params['fpr_set'][3][i] 
            for i in xrange(len(params['tpr_set'][3]))]

    max_dist = max(dist)
    threshold_indx = dist.index(max_dist)

    threshold_result = params['threshold_set'][3][threshold_indx]

    # read gamma
    data2 = read_file(params['class_res_dir'] + str(params['final_npcs']) + \
                      'PC/hyperpar_values.dat')

    for line in data2:
        if line[0] == 'gamma':
            gamma = float(line[1])
            break

    screen('Number of PCs: ' + str(params['final_npcs']), params)
    screen('gamma: ' + str(gamma), params)
    screen('threshold: ' + str(threshold_result), params)

    op3 = open(params['class_res_dir'] + 'ROC_results.dat', 'w')
    op3.write('ncomp: ' + str(params['final_npcs']) + '\n')
    op3.write('gamma: ' + str(gamma) + '\n')
    op3.write('threshold: ' + str(threshold_result))
    op3.close()

def calc_binned_diag(params):
    """
    Calculate efficiency, purity and figure of merit per redshift bin.

    input: params, dict
           keywords, value type:
               class_res, dict: results from classification
               config, dict: optimal parameter configurations
               orig_pop, dict: snid: [sample, type, redshift]
               nbins, int: number of redshift bins
               orig_Ia_bins, list: number of SN Ia per redshift bin in
                                   original sample

    output: params, dict
            additional keywords, value type:
                eff, list: global efficiency after cuts per redshift bin
                pur, list: global purity per redshift bin
                fom, list: global figure of merity after cuts per redshift bin
                effB, list: global efficiency before cuts per redshift bin
                fomb, list: global figure of merit before cuts per z bin
    """
    prob_limit =  params['config']['threshold:']

    # calculate diagnostic per redshift bin
    cc_Ia_bins = [sum([1 for key in params['class_res'].keys() 
                        if params['class_res'][key][1] >=  prob_limit and \
                        params['class_res'][key][0] == 'Ia' and \
                        (i * params['dz']) <= params['orig_pop'][key][2] and \
                        params['orig_pop'][key][2] < (i + 1) * params['dz']])
                        for i in xrange(params['nbins'])]

    tot_Ia_bins = [sum([1 for key in params['class_res'].keys() 
                        if params['class_res'][key][0] == 'Ia' and \
                        (i * params['dz']) <= params['orig_pop'][key][2] and \
                        params['orig_pop'][key][2]  < (i + 1) * params['dz']])  
                        for i in xrange(params['nbins'])]

    tot_nonIa_bins = [sum([1 for key in params['class_res'].keys() 
                      if params['class_res'][key][0] in ['Ibc','II'] and \
                      (i * params['dz']) <= params['orig_pop'][key][2] and \
                      params['orig_pop'][key][2] < (i + 1) * params['dz']])
                      for i in xrange(params['nbins'])]

    wc_Ia_bins = [sum([1 for key in params['class_res'].keys() 
                  if params['class_res'][key][1] >= prob_limit and \
                  params['class_res'][key][0] in ['Ibc', 'II'] and \
                  (i * params['dz']) <= params['orig_pop'][key][2] and \
                  params['orig_pop'][key][2] < (i + 1) * params['dz']])
                  for i in xrange(params['nbins'])]

    eff_bin = []
    pur_bin = []
    fom_bin = []
    effb_bin = []
    fomb_bin = []

    for i in xrange(params['nbins']):
        if tot_Ia_bins[i] > 0 and (cc_Ia_bins[i] + wc_Ia_bins[i]) > 0 and \
        params['orig_Ia_bins'][i] > 0:

            eff_temp = cc_Ia_bins[i] / float(tot_Ia_bins[i])
            pur_temp = float(cc_Ia_bins[i]) / (cc_Ia_bins[i] + wc_Ia_bins[i])
            effb_temp = cc_Ia_bins[i] / float(params['orig_Ia_bins'][i])

            eff_bin.append(eff_temp)
            pur_bin.append(pur_temp)
            fom_bin.append(eff_temp * pur_temp)
            effb_bin.append(effb_temp)
            fomb_bin.append(effb_temp * pur_temp)

        else:
            eff_bin.append(0)
            pur_bin.append(0)
            fom_bin.append(0)
            effb_bin.append(0)
            fomb_bin.append(0)


    params['eff_bin'] = eff_bin
    params['pur_bin'] = pur_bin
    params['fom_bin'] = fom_bin
    params['effb_bin'] = effb_bin
    params['fomb_bin'] = fomb_bin

    return params

def calc_global_diag(params):
    """
    Calculate global efficiency, purity and figure of merit.

    input: params, dict
           keywords, value type:
               class_res, dict: results from classification
               config, dict: optimal parameter configurations
               orig_photo_Ia, int: number of SN Ia in original sample

    output: params, dict
            additional keywords, value type:
                eff, float: global efficiency after cuts
                pur, float: global purity
                fom, float: global figure of merity after cuts
                effB, float: global efficiency before cuts
                fomb, float: global figure of merit before cuts
               
    """
    prob_limit = params['config']['threshold:']

    # calculate global diagnostic
    cc_Ia = sum([1 for key in params['class_res'].keys() 
                 if params['class_res'][key][1] >= prob_limit and \
                 params['class_res'][key][0] == 'Ia'])

    tot_Ia = sum([1 for key in params['class_res'].keys() 
                  if params['class_res'][key][0] == 'Ia'])

    tot_nonIa = sum([1 for key in params['class_res'].keys() 
                     if params['class_res'][key][0] in ['Ibc','II']])

    wc_Ia = sum([1 for key in params['class_res'].keys() 
                 if params['class_res'][key][1] >= prob_limit and \
                 params['class_res'][key][0] in ['Ibc', 'II']])

    params['eff'] = cc_Ia/float(tot_Ia)
    params['pur'] = float(cc_Ia) / (cc_Ia + wc_Ia)
    params['fom'] =  params['eff'] * params['pur']

    params['effb'] = cc_Ia / float(params['orig_photo_Ia'])
    params['fomb'] = params['effb'] * params['pur']

    return params

def count_pop(params, type_number):
    """
    Count original population from spec and photo samples.

    input: params, dict
           keywords, value type:
               user_input, str: path to user input file
               nbins, int: number of redshift bins
               dz, float: width of redshift bin

           type_number, dict
           dictionary to translate types between raw data and final
           classification
           keywords -> final classificaton elements
           values -> identifiers in raw data

    output: params, dict
            additional keywords, value type:
                orig_photo_Ia, int: number of SN Ia before cuts
                orig_Ia_bins, list: number of SN Ia per redshift bin
                                    before cuts
                orig_pop, dict: snid: [sample, type, redshift]
    """
    from snclass.util import read_user_input, read_snana_lc

    import os
    import numpy as np

    # count original population
    user_input = read_user_input(params['user_input'])

    # check reference filter
    if 'ref_filter' in params.keys():
        user_input['ref_filter'] = params['ref_filter']

    raw_list = os.listdir(user_input['path_to_obs'][0])

    orig_pop = {}
    for sn in raw_list:
        if '.DAT' in sn:
            user_input['path_to_lc'] = [sn]
            lc = read_snana_lc(user_input)
    
            z = float(lc['REDSHIFT_FINAL:'][0])
            if lc['SNTYPE:'][0] == '-9':
                samp = 'photo'
            else:
                samp = 'spec'
            for name in type_number.keys():
                if lc['SIM_NON1a:'][0] in type_number[name]:
                    label = name
            orig_pop[lc['SNID:'][0]] = [samp, label, z]

    params['orig_pop'] = orig_pop

    all_sn = np.array(orig_pop.values())
    orig_photo_Ia = sum([1 for key in orig_pop.keys() 
                         if orig_pop[key][0] == 'photo' and \
                         orig_pop[key][1] == 'Ia'])

    params['orig_photo_Ia'] = orig_photo_Ia

    orig_Ia_bins = []
    for i in xrange(params['nbins']):
        cont = 0
        for key in orig_pop.keys():
            if orig_pop[key][0] == 'photo' and orig_pop[key][1] == 'Ia':
                if  (orig_pop[key][2] >= i * params['dz']) and  \
                (orig_pop[key][2] < (i + 1) * params['dz']):
                    cont = cont + 1

        orig_Ia_bins.append(cont)

    params['orig_Ia_bins'] = orig_Ia_bins

    return params

def read_literature_results(params):
    """
    Read results from the literature build on post-SNPCC data.

    input: params, dict
           keywords, value type:
               dz, float: width of redshift bin
               path_Spl, str: path to Ishida & de Souza, 2013 results
               path_karp, str: path to Karpenka et al, 2012 results
               path_rich, str: path to Richards et al, 2011 results

    output: params, dict
            additional keywords, value type:
                kpca_spl_dic, dict or None: results from Ishida & deSouza,2013
                -> from Richards2011
               rich_eff_data, array or None: efficiency results 
               rich_pur_data, array or None: purity results
               rich_fom_data, array or None: figure of merit
               -> from Karpenka, 2012
               karp_eff_data, array or None: efficiency results
               karp_pur_data, array or None: purity results
               karp_fom_data, array or None: figure of merit                
    """
    from snclass.algorithm import read_file

    import numpy as np

    # initiate keywords
    params['kpca_spl_dic'] = None
    params['karp_eff_data'] = None
    params['karp_pur_data'] = None
    params['karp_fom_data'] = None
    params['rich_eff_data'] = None
    params['rich_pur_data'] = None
    params['rich_fom_data'] = None

    # determine maximum redshift
    zmax = params['dz'] * params['nbins']

    if params['path_Spl'] is not None:
        # Ishida & de Souza, 2013
        xaxis = np.arange(params['dz']/2, zmax, params['dz'])
        kpca_spl = np.array(read_file(params['path_Spl']))
        params['kpca_spl_dic'] = dict([[kpca_spl[0][i] , [float(kpca_spl[j][i]) 
                                        for j in range(1, len(kpca_spl))]] 
                                        for i in range(len( kpca_spl[0]))])

    # Karpenka et al., 2012
    karp_eff_str = read_file(params['path_karp'] + 'eff.dat')
    params['karp_eff_data'] = np.array([[float(item) for item in line]  
                              for line in karp_eff_str[1:] 
                              if float(line[0]) <= (zmax - params['dz']/2)])

    karp_pur_str = read_file(params['path_karp'] + 'pur.dat')
    params['karp_pur_data'] = np.array([[float(item) for item in line]
                              for line in karp_pur_str[1:]
                              if float(line[0]) <= (zmax - params['dz']/2)])

    karp_fom_str = read_file(params['path_karp'] + 'fom.dat')
    params['karp_fom_data'] = np.array([[float(item) for item in line]
                              for line in karp_fom_str[1:]
                              if float(line[0]) <= (zmax - params['dz']/2)])

    # Richards et al., 2011
    richards_eff_str = read_file(params['path_rich'] + 'eff.dat')
    params['rich_eff_data'] = np.array([[float(item) for item in line]
                              for line in richards_eff_str[1:]
                              if float(line[0]) <= (zmax - params['dz'] / 2)])

    richards_pur_str = read_file(params['path_rich'] + 'pur.dat')
    params['rich_pur_data'] = np.array([[float(item) for item in line]
                              for line in richards_pur_str[1:]
                              if float(line[0]) <= (zmax - params['dz']/2)])

    richards_fom_str = read_file(params['path_rich'] + 'fom.dat')
    params['rich_fom_data'] = np.array([[float(item) for item in line]
                              for line in richards_fom_str[1:]
                              if float(line[0]) <= (zmax - params['dz']/2)])

    return params
 
def plot_diagnostics_aftercuts(params):
    """
    Build diagnostic plots considering only objects surviving selection cuts.

    *********
    WARNING!!
    comparison with Ishida & de Souza, 2013 is valid only if comparing with
    post-SNPCC data results!
    *********

    input: params, dict
           keywords, value type:
               dz, float: redshift bin for plotting
               representation, str: synthetic spectroscopic sample 
                                   'original', 'balanced' or 'representative'
               kpca_spl_dic, dict or None: results from Ishida & de Souza,2013
               eff_bin, list: efficiency after cuts in redshift bins
               pur_bin, list: purity in redshift bins
               fom_bin, list: fom after cuts in redshift bins
               eff, float: global efficiency after cuts
               pur, float: global purity
               fom, float: global figure of merity after cuts
               xmin, str: minimum MJD for title
               xmax, str: max MJD for title
               ref_fils, str: reference filter for title
               SNR_cut, str: SNR cut for title              
    """
    import numpy as np
    import pylab as plt
    import sys

    # determine number of bins
    for k in xrange(len(params['eff_bin'])):
        if params['pur_bin'][k] == 0:
            nbins = k
            break

    # determine reference filter for title
    if params['ref_filter'] == None:
        ref_filter = 'global'
    else:
        ref_filter = params['ref_filter']

    if params['path_Spl'] is not None:
        # determine horizontal axis for Ishida & de Souza, 2013 results
        zmax = len(params['kpca_spl_dic']['eff']) * 0.2
        xaxis = np.arange(params['dz']/2, zmax, params['dz'])

    plt.figure()
    plt.subplot(1,3,1)
    plt.plot([(i * params['dz']) + 0.5 * params['dz']
              for i in xrange(nbins)], params['eff_bin'][:nbins],
              color='red', lw=4)
    if params['dz'] == 0.2 and params['kpca_spl_dic'] is not None:
        plt.plot(xaxis, params['kpca_spl_dic']['eff'][:len(xaxis)],
                 color='blue', ls='--', lw=4)
    ax1 = plt.gca()
    ax1.annotate('eff =' + str(round(params['eff'], 2)), xy=(0.0,0.0),
                 xytext=(0.55 ,0.05), fontsize=14, fontweight='bold')
    plt.xlim(0, 0.8)
    plt.ylim(0,1.1)
    plt.xlabel('redshift', fontsize=14)
    plt.ylabel('eff', fontsize=14)

    plt.subplot(1,3,2)
    plt.title(params['representation'] + ' spec sample, [-' + \
              params['epoch_min'] + ', ' + params['epoch_max'] + \
              '] days since max, ref ' + ref_filter + ' SNR' + \
              params['SNR_cut'] + ', ' + params['filters'], fontsize=14)
    if params['dz'] == 0.2 and params['kpca_spl_dic'] is not None:
        plt.plot(xaxis, params['kpca_spl_dic']['pur'][:len(xaxis)], 
                 color='blue', ls='--', lw=4)
    plt.plot([(i * params['dz']) + 0.5 * params['dz'] 
              for i in xrange(nbins)], params['pur_bin'][:nbins],
              color='red', lw=4)
    ax2 = plt.gca()
    ax2.annotate('pur =' + str(round(params['pur'], 2)), xy=(0.0,0.0),
                 xytext=(0.55 ,0.05), fontsize=14, fontweight='bold')
    plt.xlim(0, 0.8)
    plt.ylim(0,1.1)
    plt.xlabel('redshift', fontsize=14)
    plt.ylabel('purity', fontsize=14)

    plt.subplot(1,3,3)
    plt.plot([(i * params['dz']) + 0.5 * params['dz']
              for i in xrange(nbins)], params['fom_bin'][:nbins],
              color='red', lw=4, label='this work')
    if params['dz'] == 0.2 and params['kpca_spl_dic'] is not None:
        plt.plot(xaxis, params['kpca_spl_dic']['fom'][:len(xaxis)],
                 color='blue', ls='--', lw=4, label='IdS2013')
    ax3 = plt.gca()
    ax3.annotate('fom =' + str(round(params['fom'], 2)), xy=(0.0,0.0),
                 xytext=(0.55 ,0.05), fontsize=14, fontweight='bold')
    plt.xlim(0, 0.8)
    plt.ylim(0,1.1)
    plt.xlabel('redshift', fontsize=14)
    plt.ylabel('fom', fontsize=14)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

def plot_diagnostics_before_cuts(params):
    """
    Build diagnostic plots considering all objects in original data set.

    *********
    WARNING!!
    comparison with  Richards et al., 2011
                     Karpenka et al., 2012
                     Ishida & de Souza, 2013 

    are valid only if comparing with post-SNPCC data results!
    *********

    input: params, dict
           keywords, value type:
               dz, float: redshift bin for plotting
               representation, str: synthetic spectroscopic sample 
                                   'original', 'balanced' or 'representative'
               kpca_spl_dic, dict or None: results from Ishida & de Souza,2013
               effB_bin, list: efficiency per redshift bin before cuts
               pur_bin, list: purity per redshift bin
               fomb_bin, list: figure of merity per redshift bin before cuts
               effB, float: global efficiency before cuts
               pur, float: global purity
               -> from Richards2011
               rich_eff_data, array or None: efficiency results 
               rich_pur_data, array or None: purity results
               rich_fom_data, array or None: figure of merit
               -> from Karpenka, 2012
               karp_eff_data, array or None: efficiency results
               karp_pur_data, array or None: purity results
               karp_fom_data, array or None: figure of merit
    """
    import pylab as plt
    import numpy as np

    # determine number of bins
    for k in xrange(len(params['effb_bin'])):
        if params['eff_bin'][k] == 0:
            nbins = k
            break

    # determine reference filter for title
    if params['ref_filter'] == None:
        ref_filter = 'global'
    else:
        ref_filter = params['ref_filter']

    if params['kpca_spl_dic'] is not None:
        # determine horizontal axis for Ishida & de Souza, 2013 results
        zmax = len(params['kpca_spl_dic']['eff']) * 0.2
        xaxis = np.arange(params['dz']/2, zmax, params['dz']) 

    plt.figure()
    plt.subplot(1,3,1)
    plt.plot([(i * params['dz']) + 0.5 * params['dz']
              for i in xrange(nbins)], params['effb_bin'][:nbins],
              color='red', lw=4)
    if params['dz'] == 0.2 and params['kpca_spl_dic'] is not None:
        plt.plot(xaxis, params['kpca_spl_dic']['effB'][:len(xaxis)],
                 color='blue', ls='--', lw=4)
    if params['rich_eff_data'] is not None:
        plt.plot(params['rich_eff_data'][:,0],
                 params['rich_eff_data'][:,1],
                 color='brown', ls=':', lw=4)
    if params['karp_eff_data'] is not None:
        plt.plot(params['karp_eff_data'][:,0], params['karp_eff_data'][:,1],
                 color='green', ls='-.', lw=4)
    ax4 = plt.gca()
    ax4.annotate('effb =' + str(round(params['effb'], 2)), xy=(0.0,0.0),
                 xytext=(0.05 ,1.0), fontsize=14, fontweight='bold')
    plt.xlim(0, 0.8)
    plt.ylim(0,1.1)
    plt.xlabel('redshift', fontsize=14)
    plt.ylabel('effb', fontsize=14)

    plt.subplot(1,3,2)
    plt.title(params['representation'] + ' spec sample, [-' + \
              params['epoch_min'] + ', ' + params['epoch_max'] + \
              '] days since max, ref ' + ref_filter + ' SNR' + \
              params['SNR_cut'] +', ' + params['filters'], fontsize=14)
    plt.plot([(i * params['dz']) + 0.5 * params['dz'] for i in xrange(nbins)],
              params['pur_bin'][:nbins], color='red', lw=4)
    if params['dz'] == 0.2 and params['kpca_spl_dic'] is not None:
        plt.plot(xaxis, params['kpca_spl_dic']['pur'][:len(xaxis)],
                 color='blue', ls='--', lw=4)
    if params['rich_pur_data'] is not None:
        plt.plot(params['rich_pur_data'][:,0],
                 params['rich_pur_data'][:,1], color='brown', ls=':', lw=4)
    if params['karp_pur_data'] is not None:
        plt.plot(params['karp_pur_data'][:,0], params['karp_pur_data'][:,1],
                 color='green', ls='-.', lw=4)
    ax5 = plt.gca()
    ax5.annotate('pur =' + str(round(params['pur'], 2)), xy=(0.0,0.0),
                 xytext=(0.55 ,0.05), fontsize=14, fontweight='bold' )
    plt.xlim(0, 0.8)
    plt.ylim(0,1.1)
    plt.xlabel('redshift', fontsize=14)
    plt.ylabel('purity', fontsize=14)

    plt.subplot(1,3,3)
    plt.plot([(i * params['dz']) + 0.5 * params['dz'] 
              for i in xrange(nbins)], params['fomb_bin'][:nbins],
              color='red', lw=4, label='this work')
    if params['dz'] == 0.2 and params['kpca_spl_dic'] is not None:
        plt.plot(xaxis, params['kpca_spl_dic']['fomB'][:len(xaxis)],
                 color='blue', ls='--', lw=4, label='IdS2013')
    if params['rich_fom_data'] is not None:
        plt.plot(params['rich_fom_data'][:,0], params['rich_fom_data'][:,1],
                 color='brown', ls=':', lw=4, label='Richards (DM)')
    if params['karp_fom_data'] is not None:
        plt.plot(params['karp_fom_data'][:,0], params['karp_fom_data'][:,1],
                 color='green', ls='-.', lw=4,  label='Karpenka (NN)')
    ax6 = plt.gca()
    ax6.annotate('fomb =' + str(round(params['fomb'], 2)), xy=(0.0,0.0),
                 xytext=(0.55 ,1.0), fontsize=14, fontweight='bold')
    plt.xlim(0, 0.8)
    plt.ylim(0,1.1)
    plt.xlabel('redshift', fontsize=14)
    plt.ylabel('fomb', fontsize=14)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

