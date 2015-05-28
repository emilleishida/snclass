#!/usr/bin/env python

import argparse
import os
import sys

import numpy as np
from multiprocessing import Pool

from treat_lc import LC
from util import read_user_input, choose_sn, read_snana_lc
from functions import core_cross_val


##############################################

class DataMatrix(object):
    """
    Data matrix object.
    """

    def __init__(self, input_file=None):
        """
        Read user input file. 

        input: input_file -> str
               name of user input file
        """  

        if input_file is not None:
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
        redshift = []
        sntype = [] 

        for obj in file_list:
            if 'mean' in obj:

                #take object identifier
                name = obj[len('DES_SN'):-len('_mean.dat')]

                if len(name) == 5:
                    name = '0' + name
                elif len(name) == 4:
                    name = '00' + name
  
                self.user_choices['path_to_lc'] = ['DES_SN' + name + '.DAT']

                #read light curve raw data
                raw = read_snana_lc(self.user_choices)

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
                    redshift.append(raw[self.user_choices['redshift_flag'][0]][0])
                    sntype.append(raw[self.user_choices['type_flag'][0]][0])

        self.datam = np.array(datam)
        self.redshift = np.array(redshift)
        self.sntype = np.array(sntype)

        #write to file 
        if file_out is not None:     
            op1 = open(file_out, 'w')
            op1.write('SNID    type    z   LC...\n')
            for i in xrange(len(datam)):
                op1.write(str(self.snid[i]) + '    ' + str(self.sntype[i]) + '    ' + str(self.redshift[i]) + '    ')
                for j in xrange(len(datam[i])):
                    op1.write(str(datam[i][j]) + '    ')
                op1.write('\n')
            op1.close()  

    def reduce_dimension(self):
        """
        Perform dimensionality reduction with user defined funciton. 

        input: pars - dict
               Dictionary of parameters.
               Must include all keywords required by 
               self.user_choices['dim_reduction_func']() function.
        """

        self.low_dim_matrix = self.user_choices['dim_reduction_func'](self.datam, self.user_choices)
        self.transf_test = self.user_choices['dim_reduction_func'](self.datam,  self.user_choices, transform=True)
    
    def cross_val(self):
        """
        Optimize the hyperparameters for RBF kernel and number 
        of components in kPCA analysis.     
        """ 

        #correct type parameters if necessary
        if self.user_choices['transform_types_func'] is not None:
            self.sntype = self.user_choices['transform_types_func'](self.sntype)
   
        #initialize parameters
        data = self.datam
        types = self.sntype
        choices = self.user_choices 
        
        parameters=[[data, types, choices] for attempt in xrange(self.user_choices['n_cross_val_particles'])]

        if int(self.user_choices['n_proc'][0]) > 0:
            pool = Pool(processes=int(self.user_choices['nproc'][0]))
            p = pool.map_async(self.user_choices['cross_validation_func'], parameters)
            try:
                 results = p.get(0xFFFF)
            except KeyboardInterrupt:
                print 'Interruputed by the user!'
                sys.exit()

            pool.close()
            pool.join()

        else:
            results = np.array([core_cross_val(data, types, choices) for attempt in xrange(self.user_choices['n_cross_val_particles'])])
            
        indx_max = list(results[:,-1]).index(max(results[:,-1]))
        
        self.final = {}
        for l1 in xrange(len(self.user_choices['cross_val_par'])):
            self.final[self.user_choices['cross_val_par'][l1]] = results[indx_max][l1]
        
def main():
  print(__doc__)

if __name__=='__main__':
  main()     

