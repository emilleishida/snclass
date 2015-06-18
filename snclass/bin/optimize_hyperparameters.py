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
Optimize parameters for the kernel, perform dimensionality reduction
and project one photometric object.

Usage:

$ proj_kpca_phot.py -i <user_input_file> -t <mean_GP_fit_file_for_photo_obj>
                    -d <output_directory> -c X - p Y

if X == 1, calculate samples, otherwise read them from file
if Y == 1, generate output plot file (optional)
"""

#!/usr/bin/env python

import argparse

from snclass.matrix import DataMatrix
from snclass.matrix_classification import classify_test, store_test_matrix
from snclass.matrix_classification import test_samples

def main(args):
    """
    Optimize kernel parameters and classify one photometric light curve.
    """
    # train spectroscopic matrix
    my_matrix = DataMatrix(args.input)
    my_matrix.user_choices['n_samples'] = ['0']
    my_matrix.build()
    my_matrix.cross_val()
    my_matrix.final_configuration()

    print 'final configuration: ' + str(my_matrix.final)

    new_lc = classify_test(args.test, my_matrix, my_matrix.user_choices,
                           test_dir=args.dir, csamples=bool(int(args.calc)))

    if args.plot:
        #create test dictionary information
        test_point = {}
        test_point['data'] = new_lc.test_proj
        test_point['type'] = 'Ia' if new_lc.raw['SIM_NON1a:'][0] == '0' else 'nonIa'

        pcs = [int(args.pc1), int(args.pc2)]
        my_matrix.plot(pcs, new_lc.user_choices['file_root'][0] +
                       new_lc.raw['SNID:'][0] + '.png', test=test_point)

if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Optimization and ' + \
                                     'classification.')
    parser.add_argument('-i','--input', help='Input file name', 
                        required=True)
    parser.add_argument('-t', '--test',
                        help='File with GP mean fit for test obj.', required=True) 
    parser.add_argument('-c', '--calc', 
                        help='If true, calculate samples', required=True)
    parser.add_argument('-d', '--dir', help='Output directory.', required=True)
    parser.add_argument('-p', '--plot',
                             help='If true, generate plot', required=False)
    args = parser.parse_args()
   
    main(args)



