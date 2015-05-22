import argparse
import os

import numpy as np
from astropy.io import ascii

from treat_lc import LC
from util import read_user_input, choose_sn, read_SNANA_lc, plot

##############################################

class DataMatrix(object):
    """
    Data matrix object.
    """

    def __init__(self, input_file):
        """
        Read user input file. 

        input: input_file -> str
               name of user input file
        """  

        self.user_choices = read_user_input(input_file)


    def build(self, file_out, plot=False):
        """
        Build data matrix according to user input file specifications.

        input:   file_out -> file to store data matrix (str)
                 plot (optional) -> rather to make save plot for each LC
        """

        if not os.path.exists(self.user_choices['samples_dir'][0]):
            os.makedirs(self.user_choices['samples_dir'][0])

        #read list of SN in sample
        snlist = ascii.read(self.user_choices['snlist'][0])

        #open output file
        op1 = open(file_out, 'w')
        op1.write('SNID    SIM_NON1a    z    ')
        for fil in self.user_choices['filters']:
            for day in np.arange(float(self.user_choices['epoch_cut'][0]), 
                                 float(self.user_choices['epoch_cut'][1]), 
                                 float(self.user_choices['epoch_bin'][0])):
                op1.write( fil + '_' + str(abs(int(day))) + '    ')
        op1.write('\n')

        cont = 0

        for sn in snlist: 

            #update object
            self.user_choices['path_to_lc'] = [sn[0]]

            #read light curve raw data
            raw = read_SNANA_lc(self.user_choices)

            #initiate light curve object
            lc = LC(raw, self.user_choices)

            print 'Fitting SN' + raw['SNID:'][0]

            #write SN identification and type
            op1.write(raw['SNID:'][0] + '    ' + raw['SIM_NON1a:'][0] + '    ' 
                      + raw['REDSHIFT_FINAL:'][0] + '     ')            

            #perform basic check
            lc.check_basic()

            #check if satisfy minimum cut
            if lc.basic_cuts:

                print '... Passed basic cuts'

                #fit 
                lc.fit_GP(samples=False)

                #normalize
                lc.normalize()

                #shift to peak mjd
                lc.mjd_shift()

                #check epoch requirements
                lc.check_epoch()

                if lc.epoch_cuts:

                    cont = cont + 1
                    print '... ... Passed epoch cuts. This is SN number ' + str(cont)

                    if int(self.user_choices['n_samples'][0]) > 0:
                         lc.fit_GP(samples=True)  
                         lc.normalize(samples=True)                      
             
                    #build data matrix lines
                    lc.build_steps()

                    #write to file
                    for fil1 in self.user_choices['filters']:
                        for elem in lc.flux_for_matrix[fil1]:
                            op1.write(str(elem) + '    ')
                    op1.write('\n')
              
                    if plot == True:
                        if int(self.user_choices['n_samples'][0]) == 0:
                            lc.plot_fitted(file_out=self.user_choices['samples_dir'][0] + raw['SNID:'][0]+'.png')
                        else:
                            lc.plot_fitted(file_out=self.user_choices['samples_dir'][0] + raw['SNID:'][0]+'.png', 
                                           samples=True, nsamples=int(self.user_choices['n_samples'][0]))

                else:
                    print '... ... Failed to pass epoch cuts!\n'

            else:
                print '... Failed to pass basic cuts!\n'
                    

        op1.close()

        

        
