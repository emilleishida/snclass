import numpy as np
import os 

#########################################

def screen(message, choices):
    """
    Print message on screen according to users choice

    input:   message, str
             message to be printed

             choices, dict
             dictionary of users choices
    """
    
    if bool(int(choices['screen'][0])):
        print message


def read_user_input(filename):
    """
    Read user input from file and construct initial dictionary parameter. 

    input:    filename (string) -> user input file parameter 

    output:   dictionary with formated user choices
    """

    #read user input data
    op1 = open(filename, 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [elem.split() for elem in lin1]     

    #store options in params dictionary
    params = dict([(line[0], line[2:line.index('#')])  
                   for line in data1 if len(line) > 1])

    #check dimensionality reduction function
    if ('dim_reduction_func' in params.keys()):
        if  (params['dim_reduction_func'][0] == 'kpca'):
            from functions import kpca
            params['dim_reduction_func'] = kpca
            for i1 in xrange(len(params['kpca_pars'])):
                try:
                    params[params['kpca_pars'][i1]] = int(params['kpca_val'][i1])
                except ValueError:
                    try:
                        params[params['kpca_pars'][i1]] = float(params['kpca_val'][i1])  
                    except ValueError:
                        params[params['kpca_pars'][i1]] = params['kpca_val'][i1]

        elif params['dim_reduction_func'][0] == 'None':
            params['dim_reduction_func'] = None

    #check classifier function
    if ('classifier_func' in params.keys()):
        if (params['classifier_func'][0] == 'nn'):
            from functions import nn
            params['classifier_func'] = nn     
            for i2 in xrange(len(params['classifier_pars'])):
                try:
                    params[params['classifier_pars'][i2]] = int(params['classifier_val'][i2])
                except ValueError:
                    try:
                        params[params['classifier_pars'][i2]] = float(params['classifier_val'][i2])
                    except ValueError:
                        params[params['classifier_pars'][i2]] = params['classifier_val'][i2]     
            
        elif params['classifier_func'][0] == 'None':
            params['classifier_func'][0] = None

    #check set types function
    if ('transform_types_func' in params.keys()):
        if (params['transform_types_func'][0] == 'set_types'):
            from functions import set_types
            params['transform_types_func'] = set_types
        elif params['transform_types_func'][0] == 'None':
            params['transform_types_func'][0] = None

    #check cross-validation function
    if ('cross_validation_func' in params.keys()):
        if (params['cross_validation_func'][0] == 'cross_val'):
            from functions import core_cross_val
            params['cross_validation_func'] = core_cross_val
            params['gamma_nparticles'] = int(params['gamma_nparticles'][0])
            params['n_cross_val_particles'] = int(params['n_cross_val_particles'][0])
            for i3 in xrange(len(params['cross_val_par'])):
                try:
                    params[params['cross_val_par'][i3] + '_lim'] = [int(params[params['cross_val_par'][i3] + '_lim'][i4]) for i4 in range(2)]
                except ValueError:
                    params[params['cross_val_par'][i3] + '_lim'] = [float(params[params['cross_val_par'][i3] + '_lim'][i4]) for i4 in range(2)]

        elif params['cross_validation_func'][0] == 'None':
            params['cross_validation_func'] = None

    params['GP_fit'] = {}
    params['realizations'] = {}
    params['xarr'] = {}
    params['GP_obj'] = {}
    params['GP_std'] = {}

    #check if ``observer'' data already exists
    if not os.path.isdir(params['path_to_obs'][0]):
        raise TypeError('Variable "path_to_obs" is not a valid directory!')            

    return params

def read_SNANA_lc( params ):
    """
    Reads light curve in SNANA format and returns a dictionary with the 
    variables chosen in the user input file.

    input:     params -> dictionary of input parameters 

    output:    mdata -> data from light curve
    """

    #read light curve data
    op1 = open(params['path_to_obs'][0] + params['path_to_lc'][0], 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [elem.split() for elem in lin1]
   
    raw_data = dict([[line[0], line[1:]] for line in data1 if len(line) > 1])

    #determine MJD index
    mjd_indx = raw_data[params['param_list'][0]].index(params['mjd_flag'][0]) + 1

    #determine filter index
    filter_indx = raw_data[params['param_list'][0]].index(params['filter_flag'][0]) + 1

    #determine photon count index
    photon_indx = raw_data[params['param_list'][0]].index(params['photon_flag'][0]) + 1

    #determine photon count error index
    photonerr_indx = raw_data[params['param_list'][0]].index(params['photonerr_flag'][0]) + 1

    #determine quality criteria index
    quality_indx = raw_data[params['param_list'][0]].index(params['quality_flag'][0]) + 1

    #build measurement list for each filter
    mdata = dict([[item, np.array([[float(line[ mjd_indx]), 
                                    float(line[photon_indx]), 
                                    float(line[photonerr_indx]), 
                                    float(line[quality_indx])] 
 			for line in data1 
			if len(line) > 1 and
                        line[0] == params['epoch_flag'][0] and 
			line[filter_indx] == item and  	
			float(line[photon_indx]) >= 0.0 and 
                        float(line[quality_indx]) >= float(params['quality_cut'][0])
                        ])] 
		for item in params['filters']])

    #add usefull header information to output dictionary
    for item in params['header']:
        if item not in params['param_list']:
            mdata[ item ] = raw_data[item]

    return mdata

def choose_sn( params, output_file='snlist.dat' ):
    """
    Read through all files within a given directory and 
    select those satisfying criteria in user input file.

    input:  params (dict)

    output: txt file with the name of objects surviving selections cuts.
    """  

    #take all file names in data directory
    filename = os.listdir( params['path_to_obs'][0] )

    #store files for light curves surviving selection cuts
    final_list = []

    for name in filename:

        if params['file_root'][0] in name:

            #read light curve data
            op1 = open( params['path_to_obs'][0] + name, 'r')
            lin1 = op1.readlines()
            op1.close()

            data1 = [ elem.split() for elem in lin1 ]

            #take header parameters
            header = {}
            for line in data1:
                if len( line ) > 1 and line[0] in params['header']:
                    header[ line[0] ] = line[1:]

            #check for request SN type
            if  params['type_cut'][0]  != 'None' :
                if header[ params['type_flag'][0]][0] in params['type_cut']:
                    type_surv = True

                else:
                    type_surv = False 
            else:
                type_surv = True

            #check for requested SN sample  
            if params['sample_cut'][0] != 'None': 

                if str(header[params['sample_flag'][0]][0]) in params['sample_cut']:
                    sample_surv = True
                else:
                    sample_surv = False

            else:
                sample_surv = True 
          
            #store only if all requirements are satisfied
            if type_surv and sample_surv:
                final_list.append( name )    

    op2 = open( output_file, 'w' )
    for item in final_list:
        op2.write( item  + '\n')
    op2.close()

    screen('Found ' + str( len( final_list ) ) + ' SN satisfying sample and type cuts.', params)
    screen('Surviving objects are listed in file ' + output_file, params)


def read_fitted(lc_data):
    """  
    Read GP results previously calculated and populate the correct keywords in the parameter dictionary.

    input:  user_input, dic
            output from read_SNANA_lc()

    output: updated dictionary of parameters.
    """ 

    loaded = {}
  
    if bool(int(lc_data['n_samples'][0])):
        op1 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + lc_data['SNID:'][0] + '_samples.dat', 'r')
        lin1 = op1.readlines()
        op1.close()

        d1 = [elem.split() for elem in lin1]


        loaded['realizations'] = {}
        loaded['xarr'] = {}
        for fil in lc_data['filters']:
            loaded['realizations'][fil] = [[float(d1[kk][jj]) for kk in xrange(len(d1)) if d1[kk][0]==fil] 
                                                for jj in xrange(2, int(lc_data['n_samples'][0]) + 2)]
            loaded['xarr'][fil] = []  
            for i1 in xrange(len(d1)):
                if d1[i1][0] == fil:
                    loaded['xarr'][fil].append(float(d1[i1][1]))
                
    op2 = open(lc_data['samples_dir'][0] + lc_data['file_root'][0] + lc_data['SNID:'][0] + '_mean.dat', 'r')
    lin2 = op2.readlines()
    op2.close()

    d2 = [elem.split() for elem in lin2]

    loaded['GP_std'] = {}
    loaded['GP_fit'] = {}
    loaded['xarr'] = {}
    for fil in lc_data['filters']:
        loaded['xarr'][fil] = [float(d2[j][1]) for j in xrange(1,len(d2)) if d2[j][0] == fil]
        loaded['GP_fit'][fil] = [float(d2[j][2]) for j in xrange(1,len(d2)) if d2[j][0] == fil]
        loaded['GP_std'][fil] = [float(d2[j][3]) for j in xrange(1,len(d2)) if d2[j][0] == fil]

    return loaded



def main():
    print(__doc__)

if __name__=='__main__':
    main()    

