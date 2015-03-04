snclass - Supernova Photometric Classifier
==========================================

``snclass`` is a package designed to perform supernova photometric classification.

It works in 3 steps:

#. Convert flux/magnitude measurements in a continuous function using Gaussian process (GP). 
#. Dimentionality reduction through kernel principal component analysis (kPCA). 
#. Classification of a purely photometric sample thorough nearst neighbor (NN).  


This preliminary version contains only the GP fitting part in the script ``~snclass/bin/fit_lc_george.py``.

Installation
************

In the ``snclass`` directory, do::

    $ python setup.py install

If you do not have root privileges, do::

    $ python setup.py install --user


Running
********

At this point the script will read the raw data given in `SNANA <http://das.sdss2.org/ge/sample/sdsssn/SNANA-PUBLIC/>`_ format
(you can find an example in ``~snclass/examples/DES_SN849359.DAT``) and many GP realizations for all available bands, generated using `george tutorial <https://github.com/dfm/george/blob/master/docs/_code/model.py>`_.

Copy the sample input file  ``~snclass/examples/fit_lc_input.dat``  and the data file ``~snclass/examples/DES_SN849359.DAT``  to your current directory.


Do not forget to change the variables ``path_to_obs`` and  ``samples_dir`` to match your system.

Look carefully through all the options in the sample input file.

This file is documented and should be pretty straight forward to interpret. 


To generate the fit and plot from the command line do::

    $ fit_lc_george.py -i example_input.dat -c 1

This will generate a file with the GP mean, another file with the GP realizations and the corresponding jpeg plot.

The ``-c`` option denotes if you want to calculate all GP realizations or if you only want to read a previous calculated result.

If you only want to take a look at a result you calculated before, do::

    $ fit_lc_george.py -i example_input.dat -c 0


This should generate a plot like this:

`Example fitted light curve using george <https://github.com/emilleishida/snclass/blob/emille_dev/examples/gp-results.png>`_




Requirements
************

* Python 2.7
* numpy >=1.8.2
* matplotlib >= 1.3.1     
* argparse >= 1.1
* george >= 0.2.1


License
********

* GNU General Public License (GPL>=3)


