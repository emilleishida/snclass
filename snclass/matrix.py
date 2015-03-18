import argparse

from snclass.prepare_LC import LC
from snclass.util import read_user_input, choose_sn, read_SNANA_lc

##############################################

def main(args):
    """
    Build data matrix for posterior kernel PCA analysis.
    """

    #read_user_input
    user_input = read_user_input(args.input) 



if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description = 'Supernova photometric classification using KPCA.')
    parser.add_argument('-i','--input', help = 'Input file name', required = True)

    args = parser.parse_args()
   
    main(args)

