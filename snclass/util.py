import numpy as np
import os 
import matplotlib.pylab as plt

#########################################

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
    params = dict([(line[0], line[2:line.index('#')])  for line in data1 if len(line) > 1])

    params['GP_fit'] = {}
    params['realizations'] = {}
    params['xarr'] = {}
    
    #check if ``observer'' data already exists
    if not os.path.isdir(params['path_to_obs'][0]):
        raise TypeError('Variable "path_to_obs" is not a valid directory!')            

    return params

def read_SNANA_lc( params ):
    """
    Reads light curve in SNANA format and returns a dictionary with the variables chosen in the user input file.

    input:     params -> dictionary of input parameters 

    output:    mdata -> data from light curve
    """

    #read light curve data
    op1 = open(params['path_to_lc'][0], 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [elem.split() for elem in lin1]
   
    raw_data = dict( [[line[0], line[1:]] for line in data1 if len(line) > 1 ])

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
    mdata = dict([[item, np.array([ [ float(line[ mjd_indx]), float(line[ photon_indx ]), float(line[ photonerr_indx]), float(line[quality_indx])] 
 			for line in data1 
			if len( line ) > 1 
			and line[0] == params['epoch_flag'][0] 
			and line[ filter_indx ] == item 	
			and float( line[ quality_indx ] ) >= float(params['quality_cut'][0])
                        ])] 
		for item in params['filters'] ])

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
            if type_surv == True and sample_surv == True:
                final_list.append( name )    

    op2 = open( output_file, 'w' )
    for item in final_list:
        op2.write( item  + '\n')
    op2.close()

    print  'Found ' + str( len( final_list ) ) + ' SN satisfying sample and type cuts.'
    print  'Surviving objects are listed in file ' + output_file


def plot(user_input, lc_data, output_file):
    """  
    Plot GP fit to light curve. 

    input:    user_input -> output from function read_user_input
              lc_data -> output from function fit_LC
              output_file -> str, name of output file to store plot
    """

    #initiate figure
    plt.figure()
    
    for fil in user_input['filters']:

        # Plot the samples in data space.
        ax = plt.subplot(len(user_input['filters']), 1, user_input['filters'].index(fil) + 1)
        for s in lc_data['realizations'][fil]:
            plt.plot(lc_data['xarr'][fil], s, color="#4682b4", alpha=0.3)
        plt.errorbar(lc_data[fil][:,0], lc_data[fil][:,1], yerr=lc_data[fil][:,2], fmt=".k", capsize=0, label=fil)
        plt.plot(lc_data['xarr'][fil], lc_data['GP_fit'][fil], 'r:', linewidth=2)
        plt.ylabel("FLUXCAL")
        plt.xlabel("MJD")
        plt.legend()
        plt.xlim(min(lc_data[fil][:,0]) - 1.0, max(lc_data[fil][:,0]) + 1.0)

    plt.suptitle("Results with Gaussian process noise model")
    plt.savefig(output_file, dpi=350)


def plot_shifted(user_input, lc_data, output_file):
    """  
    Plot GP fit to light curve. 

    input:    user_input -> output from function read_user_input
              lc_data -> output from function fit_LC
              output_file -> str, name of output file to store plot
    """

    #initiate figure
    plt.figure()
    
    for fil in user_input['filters']:

        # Plot the samples in data space.
        ax = plt.subplot(len(user_input['filters']), 1, user_input['filters'].index(fil) + 1)
        for s in lc_data.fitted['norm_realizations'][fil]:
            plt.plot(lc_data.fitted['xarr_shifted'][fil], s, color="#4682b4", alpha=0.3)
        plt.errorbar(lc_data.raw[fil][:,0] - lc_data.fitted['peak_mjd'], 
                     lc_data.raw[fil][:,1]/lc_data.fitted['max_flux'], 
                     yerr=lc_data.raw[fil][:,2]/lc_data.fitted['max_flux'], 
                     fmt=".k", capsize=0, label=fil)
        plt.plot(lc_data.fitted['xarr_shifted'][fil], 
                 lc_data.fitted['norm_fit'][fil], 'r:', linewidth=2)
        plt.ylabel("FLUXCAL")
        plt.xlabel("MJD")
        plt.legend()
        plt.xlim(min(lc_data.raw[fil][:,0]) - lc_data.fitted['peak_mjd'] - 1.0, 
                 max(lc_data.raw[fil][:,0]) - lc_data.fitted['peak_mjd'] + 1.0)

    plt.suptitle("Results with Gaussian process noise model")
    plt.savefig(output_file, dpi=350)


def main():
  print(__doc__)

if __name__=='__main__':
  main()    

