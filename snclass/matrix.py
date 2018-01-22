"""
Created by Emille Ishida in May, 2015.

Class to implement calculations on data matrix.
"""

import os
import sys

import matplotlib.pylab as plt
import numpy as np
from multiprocessing import Pool

from snclass.treat_lc import LC
from snclass.util import read_user_input, read_snana_lc, translate_snid
from snclass.functions import core_cross_val, screen

##############################################


class DataMatrix(object):

    """
    Data matrix class.

    Methods:
        - build: Build data matrix according to user input file specifications.
        - reduce_dimension: Perform dimensionality reduction.
        - cross_val: Perform cross-validation.

    Attributes:
        - user_choices: dict, user input choices
        - snid: vector, list of objects identifiers
        - datam: array, data matrix for training
        - redshift: vector, redshift for training data
        - sntype: vector, classification of training data
        - low_dim_matrix: array, data matrix in KernelPC space
        - transf_test: function, project argument into KernelPC space
        - final: vector, optimize parameter values
    """

    def __init__(self, input_file=None):
        """
        Read user input file.

        input: input_file -> str
               name of user input file
        """
        self.datam = None
        self.snid = []
        self.redshift = None
        self.sntype = None
        self.low_dim_matrix = None
        self.transf_test = None
        self.final = None
        self.test_projection = []

        if input_file is not None:
            self.user_choices = read_user_input(input_file)

    def check_file(self, filename, epoch=True, ref_filter=None):
        """
        Construct one line of the data matrix.

        input:   filename, str
                 file of raw data for 1 supernova

                 epoch, bool - optional
                 If true, check if SN satisfies epoch cuts
                 Default is True

                 ref_filter, str - optional
                 Reference filter for peak MJD calculation
                 Default is None
        """
        screen('Fitting ' + filename, self.user_choices)

        # translate identifier
        self.user_choices['path_to_lc'] = [translate_snid(filename, self.user_choices['photon_flag'][0])[0]]

        # read light curve raw data
        raw = read_snana_lc(self.user_choices)

        # initiate light curve object
        lc_obj = LC(raw, self.user_choices)

        # load GP fit
        lc_obj.load_fit_GP(self.user_choices['samples_dir'][0] + filename)

        # normalize
        lc_obj.normalize(ref_filter=ref_filter)

        # shift to peak mjd
        lc_obj.mjd_shift()

        if epoch:
            # check epoch requirements
            lc_obj.check_epoch()
        else:
            lc_obj.epoch_cuts = True

        if lc_obj.epoch_cuts:
            # build data matrix lines
            lc_obj.build_steps()

            # store
            obj_line = []
            for fil in self.user_choices['filters']:
                for item in lc_obj.flux_for_matrix[fil]:
                    obj_line.append(item)

            rflag = self.user_choices['redshift_flag'][0]
            redshift = raw[rflag][0]

            obj_class = raw[self.user_choices['type_flag'][0]][0]

            self.snid.append(raw['SNID:'][0])

            return obj_line, redshift, obj_class

        else:
            screen('... Failed to pass epoch cuts!', self.user_choices)
            screen('\n', self.user_choices)
            return None

    def store_training(self, file_out):
        """
        Store complete training matrix.

        input: file_out, str
               output file name
        """
        # write to file
        if file_out is not None:
            op1 = open(file_out, 'w')
            op1.write('SNID    type    z   LC...\n')
            for i in xrange(len(self.datam)):
                op1.write(str(self.snid[i]) + '    ' + str(self.sntype[i]) +
                          '    ' + str(self.redshift[i]) + '    ')
                for j in xrange(len(self.datam[i])):
                    op1.write(str(self.datam[i][j]) + '    ')
                op1.write('\n')
            op1.close()

    def build(self, file_out=None, check_epoch=True, ref_filter=None):
        """
        Build data matrix according to user input file specifications.

        input:   file_out -> str, optional
                 file to store data matrix (str). Default is None

                 check_epoch -> bool, optional
                 If True check if SN satisfies epoch cuts
                 Default is True

                 ref_filter -> str, optional
                 Reference filter for MJD calculation
                 Default is None
        """
        # list all files in sample directory
        file_list = os.listdir(self.user_choices['samples_dir'][0])

        datam = []
        redshift = []
        sntype = []

        for obj in file_list:
            if 'mean' in obj:
                sn_char = self.check_file(obj, epoch=check_epoch,
                                          ref_filter=ref_filter)
                if sn_char is not None:
                    datam.append(sn_char[0])
                    redshift.append(sn_char[1])
                    sntype.append(sn_char[2])

        self.datam = np.array(datam)
        self.redshift = np.array(redshift)
        self.sntype = np.array(sntype)

        # store results
        self.store_training(file_out)

    def reduce_dimension(self):
        """Perform dimensionality reduction with user defined function."""
        # define dimensionality reduction function
        func = self.user_choices['dim_reduction_func']

        # reduce dimensionality
        self.low_dim_matrix = func(self.datam, self.user_choices)

        # define transformation function
        self.transf_test = func(self.datam, self.user_choices, transform=True)

    def cross_val(self):
        """Optimize the hyperparameters for RBF kernel and ncomp."""
        # correct type parameters if necessary
        types_func = self.user_choices['transform_types_func']
        if types_func is not None:
            self.sntype = types_func(self.sntype, self.user_choices['Ia_flag'][0])

        # initialize parameters
        data = self.datam
        types = self.sntype
        choices = self.user_choices

        nparticles = self.user_choices['n_cross_val_particles']
        parameters = []
        for i in xrange(nparticles):
            pars = {}
            pars['data'] = data
            pars['types'] = types
            pars['user_choices'] = choices

            parameters.append(pars)

        if int(self.user_choices['n_proc'][0]) > 0:
            cv_func = self.user_choices['cross_validation_func']
            pool = Pool(processes=int(self.user_choices['n_proc'][0]))
            my_pool = pool.map_async(cv_func, parameters)
            try:
                results = my_pool.get(0xFFFF)
            except KeyboardInterrupt:
                print 'Interruputed by the user!'
                sys.exit()

            pool.close()
            pool.join()

            results = np.array(results)

        else:
            number = self.user_choices['n_cross_val_particles']
            results = np.array([core_cross_val(pars) 
                                for pars in parameters])

        flist = list(results[:,len(results[0])-1])
        max_success = max(flist)
        indx_max = flist.index(max_success)

        self.final = {}
        for i in xrange(len(self.user_choices['cross_val_par'])):
            par_list = self.user_choices['cross_val_par']
            self.final[par_list[i]] = results[indx_max][i]

    def final_configuration(self):
        """Determine final configuraton based on cross-validation results."""
        #update optimized hyper-parameters
        for par in self.user_choices['cross_val_par']:
            indx = self.user_choices['cross_val_par'].index(par)
            self.user_choices[par] = self.final[par]

        #update low dimensional matrix
        self.reduce_dimension()       

    def plot(self, pcs, file_out, show=False, test=None):
        """
        Plot 2-dimensional scatter of data matrix in kPCA space.

        input: pcs, vector of int
               kernel PCs to be used as horizontal and vertical axis

               file_out, str
               file name to store final plot

               show, bool, optional
               if True show plot in screen
               Default is False

               test, dict, optional
               keywords: data, type
               if not None plot the projection of 1 photometric object
               Default is None
        """
        #define vectors to plot
        xdata = self.low_dim_matrix[:,pcs[0]]
        ydata = self.low_dim_matrix[:,pcs[1]]

        if '0' in self.sntype:
            snIa = self.sntype == '0'
            nonIa = self.sntype != '0'
        else:
            snIa = self.sntype == 'Ia'
            snIbc = self.sntype == 'Ibc'
            snII = self.sntype == 'II'

        plt.figure(figsize=(10,10))
        if '0' in self.sntype:
            plt.scatter(xdata[nonIa], ydata[nonIa], color='purple', marker='s',
                        label='spec non-Ia')
            plt.scatter(xdata[snIa], ydata[snIa], color='blue', marker='o',
                        label='spec Ia')
        else:
            plt.scatter(xdata[snII], ydata[snII], color='purple', marker='s',
                        label='spec II')
            plt.scatter(xdata[snIbc], ydata[snIbc], color='green', marker='^',
                        s=30, label='spec Ibc')
            plt.scatter(xdata[snIa], ydata[snIa], color='blue', marker='o',
                        label='spec Ia')

        if test is not None:
            if len(test.samples_for_matrix) > 0:
                plt.title('prob_Ia = ' + str(round(test['prob_Ia'], 2)))

            if test.raw['SIM_NON1a:'][0] == '0':
                sntype = 'Ia'
            else:
                sntype = 'nonIa'

            plt.scatter([test.test_proj[0][pcs[0]]], [test.test_proj[0][pcs[1]]],
                        marker='*', color='red', s=75,
                        label='photo ' + sntype)
        plt.xlabel('kPC' + str(pcs[0] + 1), fontsize=14)
        plt.ylabel('kPC' + str(pcs[1] + 1), fontsize=14)
        plt.legend(fontsize=12)
        if show:
            plt.show()
        if file_out is not None:
            plt.savefig(file_out)
        plt.close()




def main():
    """Print documentation."""
    print __doc__

if __name__ == '__main__':
    main()
