import numpy as np
import gptools
import os

def fit_LC(data, samples=False):
    """
    Gaussian Process fit using gptools. 

    input:  data -> dictionary of raw data
                    output from read_SNANA_lc
                    keys: filters

    output: data -> update dictionary with new key
                    realizations: 
    """

    data['realizations'] = {}
    data['xarr'] = {}
    data['GP_std'] = {}

    for fil in data['filters']:
        print '... filter: ' + fil

        t = data[fil][:,0]
        y = data[fil][:,1]
        yerr = data[fil][:,2]

        #setup GP
        k = gptools.SquaredExponentialKernel(param_bounds=[(0, max(y)), (0, np.std(t))])
        gp = gptools.GaussianProcess(k)
        gp.add_data(t, y, err_y=yerr)
        
        data['xarr'][fil] = np.arange(min(t), max(t), 0.2)
        
        out = gp.predict(data['xarr'][fil], use_MCMC=True, 
                         num_proc=int(data['n_proc'][0]), nsamp=200, plot_posterior=False,
                         plot_chains=False, burn=100, thin=10)
       
        data['GP_fit'][fil] = out[0]
        data['GP_std'][fil] = out[1]

        if samples == True and int(data['n_samples'][0]) > 0:

            print '... ... calculate samples'

            v1 = gp.draw_sample(data['xarr'][fil], num_samp=int(data['n_samples'][0]))

            data['realizations'][fil] = v1.T
            

        gp = None
        v1 = None
        out = None 
        k = None
        del gp, out, v1, k 

    if bool(int(data['save_samples'][0])) == True:

        if samples == True:

            if not os.path.exists(data['samples_dir'][0]):
                os.makedirs(data['samples_dir'][0])

            op1 = open(data['samples_dir'][0] + data['file_root'][0] + 
                       data['SNID:'][0] + '_samples.dat', 'w')
            for fil in data['filters']:
                for i1 in xrange(len(data['xarr'][fil])):   
                    op1.write(fil + str(data['xarr'][fil][i1]) + '    ')
            op1.write('\n')
            for j1 in xrange(len(data['realizations'][fil])):
                for fil in data['filters']:
                    for k1 in xrange(len(data['xarr'][fil])):
                        op1.write(str(data['realizations'][fil][j1][k1]) + '    ')
                op1.write('\n') 
            op1.close()

        op2 = open(data['samples_dir'][0] + data['file_root'][0] + 
                   data['SNID:'][0] + '_mean.dat', 'w')
        for fil in data['filters']:
            for i2 in xrange(len(data['xarr'][fil])):   
                op2.write(fil + str(data['xarr'][fil][i2]) + '    ')
        op2.write('\n')
        for fil in data['filters']:
            for j2 in xrange(len(data['xarr'][fil])):
                op2.write(str(data['GP_fit'][fil][j2]) + '    ')

        op2.close()       

    return data
       

def main():
  print(__doc__)

if __name__=='__main__':
  main()    


