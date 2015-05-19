==========================================
snclass - Supernova Photometric Classifier
==========================================

``snclass`` is a package designed to perform supernova photometric classification.

It works in 3 steps:

* Convert flux/magnitude measurements in a continuous function using Gaussian process (GP). 
* Dimentionality reduction through kernel principal component analysis (kPCA). 
* Classification of a purely photometric sample thorough nearst neighbor (NN).  


## Installation

In the ``snclass`` directory, do::

    $ python setup.py install

If you do not have root privileges, do::

    $ python setup.py install --user

## Fitting and plotting a single light curve


At this point the script will read the raw data given in [SNANA](http://das.sdss2.org/ge/sample/sdsssn/SNANA-PUBLIC/) format
(you can find an example in ``~snclass/examples/DES_SN849359.DAT``) and many GP realizations for all available bands, generated using [george tutorial](https://github.com/dfm/george/blob/master/docs/_code/model.py).

Copy the sample input file  ``~snclass/examples/user.input``  and the data file ``~snclass/examples/DES_SN849359.DAT``  to your current directory.


Do not forget to change the variables ``path_to_obs`` and  ``samples_dir`` to match your system.

Look carefully through all the options in the sample input file.

This file is documented and should be pretty straight forward to interpret. 


To generate the fit and plot from the command line do::

    $ fit_plot_lc_george.py -i user.input -c 1

This will generate a file with the GP mean, another file with the GP realizations and the corresponding jpeg plot.

The ``-c`` option denotes if you want to calculate all GP realizations or if you only want to read a previous calculated result.

If you only want to take a look at a result you calculated before, do::

    $ fit_plot_lc_george.py -i user.input -c 0


This should generate a plot like this:

![Example fitted light curve using george] 
(https://github.com/emilleishida/snclass/blob/emille_dev/snclass/examples/gp-results.png)

Notice that the above plot is only the GP posterior for the given SN in all filters. 

In order to decide if a given object satisfy all requirements stated in the user input file, do

```python
import numpy as np
import snclass

#read user input file
user_input=snclass.read_user_input('fit_lc_input.dat')

#read raw data
lc_data = snclass.read_SNANA_lc(user_input)

#update data object
lc_data.update(user_input)

#create LC object
lc = snclass.LC(lc_data, user_input)

#check SNR and number of epochs cuts
lc.check_basic()
if lc.basic_cuts == True:
    
    #fit GP u
    lc.fit_GP()
 
    #normalize according to larger flux (all filters)
    lc.normalize()

    #shift to peak MJD
    lc.mjd_shift()

    #check minimum and maximum epoch
    lc.check_epoch()
    if lc.epoch_cuts == True:
        
        #build initial data matrix elements
        lc.build_steps()
```

This will make sure the object in question pass all requirements to populate the initial data matrix. 

You can see the graphical output using

```python
lc.plot_fitted(samples=True, nsamples=int(user_input['n_samples'][0]))
```


## Identifying samples

To create a list of all SNe in the initial pool satisfying some selection cuts, update the corresponding keywords in the user input file. 
As an example, for the post-SNPCC data, in order to select all the spectroscopically classified SNe, set::

    type_flag      = SIM_NON1a:	        # type identification	
    type_cut	   = None		# type selection cut

    sample_flag	    = SNTYPE:		        # sample identification	
    sample_cut	    = 1 3 21 22 23 32 33 	# sample selection cut    


Analogously, in order to construct a list of photometric-only SNe, your user input file should contain::
	
    type_flag      = SIM_NON1a:	        # type identification	
    type_cut	   = None		# type selection cut

    sample_flag	    = SNTYPE:		# sample identification	
    sample_cut	   = -9			# sample selection cut	


The list is created iteractively with 

```python
import snclass

user_choices = snclass.read_user_input("user.input")
snclass.choose_sn(user_choices, output_file='my_sample.list')
```

The list of all SNe satisfying your selection cuts will be stored in ``my_sample.list``.


## Requirements


* Python 2.7
* numpy >=1.8.2
* matplotlib >= 1.3.1     
* argparse >= 1.1
* george >= 0.2.1


## License


* GNU General Public License (GPL>=3)


