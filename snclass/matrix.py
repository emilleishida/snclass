import argparse
import os

import numpy as np

from treat_lc import LC
from util import read_user_input, choose_sn, read_SNANA_lc

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

    

    def build(self, file_out=None):
        """
        Build data matrix according to user input file specifications.

        input:   file_out -> str, optional
                 file to store data matrix (str). Default is None
        """

        #list all files in sample directory
        file_list = os.listdir(self.user_choices['samples_dir'][0])

        datam = []
        self.snid = []
        self.redshift = []
        self.sntype = [] 

        for obj in file_list:
            if 'mean' in obj:

                #take object identifier
                name = obj[len('DES_SN'):-len('_mean.dat')]

                if len(name) == 5:
                    name = '0' + name
  
                self.user_choices['path_to_lc'] = ['DES_SN' + name + '.DAT']

                #read light curve raw data
                raw = read_SNANA_lc(self.user_choices)

                #initiate light curve object
                lc = LC(raw, self.user_choices)

                #load GP fit
                lc.load_fit_GP()

                #normalize
                lc.normalize()

                #shift to peak mjd
                lc.mjd_shift()

                #check epoch requirements
                lc.check_epoch()
  
                if lc.epoch_cuts:                     
             
                    #build data matrix lines
                    lc.build_steps()

                    #store
                    obj_line = []
                    for fil in self.user_choices['filters']:
                        for item in lc.flux_for_matrix[fil]: 
                            obj_line.append(item)

                    datam.append(obj_line)
                    self.snid.append(raw['SNID:'][0])
                    self.redshift.append(raw[self.user_choices['redshift_flag'][0]][0])
                    self.sntype.append(raw[self.user_choices['type_flag'][0]][0])

        self.datam = np.array(datam)

        if file_out is not None:     
            op1 = open(file_out, 'w')
            op1.write('SNID    type    z   LC...\n')
            for i in xrange(len(datam)):
                op1.write(str(snid[i]) + '    ' + str(sntype[i]) + '    ' + str(redshift[i]) + '    ')
                for j in xrange(len(datam[i])):
                    op1.write(str(datam[i][j]) + '    ')
                op1.write('\n')
            op1.close()  
  
