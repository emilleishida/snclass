==========================================
snclass - Supernova Photometric Classifier
==========================================

``snclass`` is a package designed to perform supernova photometric classification.

It works in 3 steps:

* Convert flux/magnitude measurements into a continuous function using Gaussian process (GP).
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

# read user input file
user_input=snclass.read_user_input('user.input')

# read raw data
lc_data = snclass.read_snana_lc(user_input)

# update data object
lc_data.update(user_input)

# create LC object
lc = snclass.LC(lc_data, user_input)

# check SNR and number of epochs cuts
lc.check_basic()
if lc.basic_cuts == True:
    
    # fit GP  - this calculates only mean fit
    lc.fit_GP(mean=True, samples=False, screen=True)
 
    # normalize according to larger flux (all filters)
    lc.normalize()

    # shift to peak MJD
    lc.mjd_shift()

    # check minimum and maximum epoch
    lc.check_epoch()
    
    print lc.epoch_cuts
```

If the  object given in keyword ``path_to_lc`` passes all requirements to populate the initial data matrix this will return ``True``. 
In this case, you might want to calculate a number of realizations from the constrained GP. 

```python
lc.fit_GP(mean=False, samples=True, screen=True)
```

You can see the graphical output using

```python
lc.plot_fitted()
```


## Identifying samples

To create a list of all SNe in the initial pool satisfying some selection cuts, update the corresponding keywords in the user input file. 
As an example, for the post-SNPCC data, in order to select all the spectroscopically classified SNe, set::

    type_flag      = SIM_NON1a:	        # type identification
    type_cut	   = None	        	# type selection cut

    sample_flag	    = SNTYPE:		        # sample identification
    sample_cut	    = 1 3 21 22 23 32 33 	# sample selection cut


Analogously, in order to construct a list of photometric-only SNe, your user input file should contain::

    type_flag      = SIM_NON1a:	        # type identification
    type_cut	   = None		        # type selection cut

    sample_flag	    = SNTYPE:		# sample identification
    sample_cut	   = -9			    # sample selection cut


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

## Fitting a set of SN

You can also fit a set of SN sequentially.
In this case, build a ``sn.list`` file, which contains the name of the raw files for all objects you want to fit.

In the ``user.input`` file, set the keyword ``snlist`` and do

```python
snclass.treat_lc.fit_objs(user_input)
```

Make sure that the keyword ``samples_dir`` is also properly set, as the output files with mean and samples results will be stored in this directory.

## Building a data matrix

After you have all your desired light curves already fitted through a GP, with means saved in the ``samples_dir`` directory, you can easily build the data matrix using

```python
from snclass.matrix import DataMatrix

d = DataMatrix('user.input')
d.build(file_out='matrix.dat')
```

This will store the complete training data matrix (one row for each object, each row a concatenation of light curves in different filters) in ``d.datam``, the corresponding objects classification in ``d.sntypes`` and will print the complete table in ``matrix.dat`` file.

## Dimensionality reduction and classifier

The current version of ``snclass``  uses Kernel Principal Component Analysis ([KernelPCA](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.KernelPCA.html)) for dimensionality reduction and [1 Nearst Neighbor](http://scikit-learn.org/stable/modules/neighbors.html) algorithm as a classifier.

However, those functions can be replaced if the user provides personalized functions and set the necessary keywords in the ``user.input`` file.

In order to use the built-in dimensionality reduction KernelPCA, the necessary keywords are

```python
dim_reduction_func = kpca                 # name of dimensionality reduction function
kpca_pars          = kernel gamma ncomp   # parameters for dimensionality reduction
kpca_val           = rbf  1.0   2         # value for dimensionality reduction parameters
```

Then, we can reduce the dimensionality of the data matrix simply doing

```python
d.reduce_dimension()
```

This will only reduce the dimensionality of the training sample and calculate the corresponding projections in ``ncomp`` KernelPCs.
Suppose we have a set of photometric light curves (test sample), ``test_LCs``, which was previously submitted to a GP fit. This can be a set of diverse objects or a number of realizations from the final GP for only 1 object.
First we need to project them using the low dimension space built from the spectroscopic sample. 

```python
# project test
test_projections = d.transf_test.transform(test_LCs)
```

The classifier is given as a separate function, which in the case implemented so far requires the following keywords

```python
classifier_func  = nneighbor   # classifier function
classifier_pars  = n weights   # classifier parameters
classifier_val   = 1 distance  # values for classifier parameters
```

In order to classify the test sample, based on the KernelPCA space from the training sample, do

```python
from snclass.functions import nneighbor

new_label = nneighbor(test_LCs, d.datam, d.user_choices)
```

## Cross-validation

In most cases, the dimensionality reduction phase will depend on a few hyperparameters which should be optimized before the actual classification takes place. This is the case for the KernelPCA algorithm. The optimization is done through a [repeated random sub-sampling validation](http://en.wikipedia.org/wiki/Cross-validation_%28statistics%29), where 1/3 of the spectroscopic sample is set asside and considered as test sample.

The input file keywords controlling the cross-validation process are::

    cross_validation_func   = cross_val            # name of cross-validation function
    n_cross_val_particles   = 10                   # number of times to separate training/test set 
    cross_val_par           = ncomp  gamma         # cross_validation parameters
    ncomp_lim               = 2 11                 # limits on number of components to be tested
    gamma_lim               = 0.05  20.0           # limits on parameter hyper_par gamma
    gamma_nparticles        = 100                  # number of gamma random values to be tested

The complete cross-validation process is performed through

```python
from snclass.matrix import DataMatrix

# initiate data matrix object
my_matrix = DataMatrix('user.input')

# build matrix using all data in samples_dir directory
my_matrix.build()

# perform cross-validation
my_matrix.cross_val()

# use optimize hyper-parameters and build final low dimension representation
my_matrix.final_configuration()

# plot scatter spec sample in PCs 0 and 1
my_matrix.plot([0,1],file_out=None, show=True)
```

## Requirements

* Python 2.7
* argparse >= 1.1
* gptools >= 0.1
* matplotlib >= 1.3.1
* multiprocessing >= 0.70a1
* numpy >=1.8.2
* scikit-learn >= 0.16.0

## License

* GNU General Public License (GPL>=3)

