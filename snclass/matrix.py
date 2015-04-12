import argparse

from snclass.prepare_LC import LC
from snclass.util import read_user_input, choose_sn, read_SNANA_lc
from prepare_lc import LC

##############################################

def main(args):
    """
    Build data matrix of light curves.
    """

    #read_user_input
    user_input = read_user_input(args.input) 

    #read light curve raw data
    raw = read_SNANA_lc(user_input)

    #initiate light curve object
    lc = LC(raw_data, user_input)

    #check if satisfy minimum cut
    if lc.check_basics():

        #if so, fit 
        lc.fit_GP()

        #normalize
        lc.normalize()

        #shift to peak mjd
        lc.mjd_shift()

        #check epoch requirements
        lc.check_epoch()


if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description = 'Supernova photometric classification using KPCA.')
    parser.add_argument('-i','--input', help = 'Input file name', required = True)

    args = parser.parse_args()
   
    main(args)

