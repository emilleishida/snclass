==========================================
snclass - Supernova Photometric Classifier
==========================================

``snclass`` is a package designed to perform supernova photometric classification.

It works in 3 steps:

* Convert flux/magnitude measurements in a continuous function using Gaussian process (GP). 
* Dimentionality reduction through kernel principal component analysis (kPCA). 
* Classification of a purely photometric sample thorough nearst neighbor (NN).  


## Installation

Download and expand this repository. 
In the ``snclass`` directory, do::

    $ python setup.py install

If you do not have root privileges, do::

    $ python setup.py install --user

## Fitting and plotting a single light curve


At this point the script will read the raw data given in [SNANA](http://das.sdss2.org/ge/sample/sdsssn/SNANA-PUBLIC/) format
(you can find an example in ``~snclass/examples/DES_SN077317.DAT``) and many GP realizations for all available bands, generated using [gptools](http://gptools.readthedocs.org/en/latest/).

Copy the sample input file  ``~snclass/examples/user.input``  and the data file ``~snclass/examples/DES_SN077317.DAT``  to your current directory.


Do not forget to change the variables ``path_to_obs`` and  ``samples_dir`` to match your system.

Look carefully through all the options in the sample input file.

This file is documented and should be pretty straight forward to interpret. 


To generate the fit and plot from the command line do::

    $ fit_plot_lc.py -i user.input -c 1

This will generate a file with the GP mean, another file with the GP realizations and the corresponding plot.

The ``-c`` option denotes if you want to calculate all GP realizations or if you only want to read a previous calculated result.

If you only want to take a look at a result you calculated before, do::

    $ fit_plot_lc.py -i user.input -c 0


This will generate a plot like this:

![Example GP fitted light curve] 
(https://github.com/emilleishida/snclass/blob/emille_dev/snclass/examples/gp-results.png)

Notice that the above plot is only the GP posterior for the given SN in all filters. 

In order to decide if a given object satisfy all requirements stated in the user input file, do

```python
import numpy as np
import snclass

#read user input file
user_input=snclass.read_user_input('user.input')

#read raw data
lc_data = snclass.read_SNANA_lc(user_input)

#update data object
lc_data.update(user_input)

#create LC object
lc = snclass.LC(lc_data, user_input)

#check SNR and number of epochs cuts
lc.check_basic()
if lc.basic_cuts == True:
    
    #fit GP  - this calculates only mean fit
    lc.fit_GP(mean=True, samples=False)
 
    #normalize according to larger flux (all filters)
    lc.normalize()

    #shift to peak MJD
    lc.mjd_shift()

    #check minimum and maximum epoch
    lc.check_epoch()
    
    print lc.epoch_cuts
```

If the  object given in keyword ``path_to_lc`` passes all requirements to populate the initial data matrix this will return ``True``. 
In this case, you might want to calculate a number of realizations from the constrained GP. 

```python
lc.fit_GP(mean=False, samples=True)
```

You can see the graphical output using

```python
lc.plot_fitted()
```

### Fitting a set of SN

You can also fit a set of SN sequentially. 
In this case, build a ``sn.list`` file, which contains the name of the raw files for all objects you want to fit.
 
In the ``user.input`` file, set the keyword ``snlist`` and do
 
```python
snclass.fit_objs(user_input)
```

Make sure that the keyword ``samples_dir`` is also properly set, as the output files with mean and samples results will be stored in this directory. 

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

***
**WARNING**

The samples separated using this method where only selected through header variables (types, samples, etc.).
No calculations were made in the raw data. 
In order to select a smaller subset satisfying selection cuts which require treatment, use the ``matrix.build`` module.
***

## Building a data matrix




## Requirements


* Python 2.7
* numpy >=1.8.2
* matplotlib >= 1.3.1     
* argparse >= 1.1
* george >= 0.2.1


## License


* GNU General Public License (GPL>=3)


