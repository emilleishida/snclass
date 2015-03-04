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
        print params['path_to_obs'][0]
        raise TypeError('Variable "path_to_obs" is not a valid directory!')            

    return params



def choose_sn( params, output_file='snlist.dat' ):

    #take all file names in data directory
    filename = os.listdir( params['path1'] )

    #store files for light curves surviving selection cuts
    final_list = []

    for name in filename:

        if params['file_root'] in name:

            #read light curve data
            op1 = open( params['path1'] + name, 'r')
            lin1 = op1.readlines()
            op1.close()

            data1 = [ elem.split() for elem in lin1 ]

            #take header parameters
            header = {}
            for line in data1:
                if len( line ) > 1 and line[0] in params['usefull_header']:
                    header[ line[0] ] = line[1:]

            #check for request SN type
            if  params['type_cut']  != None :

                if header[ params['type_flag'] ][0] in params['type_cut']:
                    type_surv = True

                else:
                    type_surv = False 
            else:
                type_surv = True
                
            #check for requested SN sample 
            if params['sample_cut'] != None: 

                if header[ params['sample_flag'] ][0] in params['sample_cut']:
                    sample_surv = True  

                else:
                    sample_surv = False

            else:
                sample_surv = True 
          
            #store only if both requirements are satisfied
            if type_surv == True and sample_surv == True:
                final_list.append( name )    

    op2 = open( output_file, 'w' )
    for item in final_list:
        op2.write( item  + '\n')
    op2.close()

    print  'Found ' + str( len( final_list ) ) + ' SN satisfying sample and type cuts.'
    print  'Surviving objects are listed in file ' + output_file


def read_SNANA_lc( params ):
    """
    Reads light curve in SNANA format and returns a dictionary with the variables chosen in the user input file.

    input:     params -> dictionary of input parameters 

    output:    mdata -> data from light curve
    """

    #read light curve data
    op1 = open(params['path_to_obs'][0] + params['path_to_lc'][0], 'r')
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


def main():
  print(__doc__)

if __name__=='__main__':
  main()    

