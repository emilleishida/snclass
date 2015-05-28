"""
Supernova Photometric Classifier in Python
"""

__author__ = "E. E. O. Ishida"
__maintainer__ = "E. E. O. Ishida"
__copyright__ = "Copyright 2015"
__version__ = "0.0.1"
__email__ = "emilleishida@gmail.com"
__status__ = "Prototype"
__license__ = "GPL"


from util import read_user_input, choose_sn, read_snana_lc, read_fitted
from fit_lc_gptools import fit_LC
from treat_lc import LC, fit_objs
from matrix import DataMatrix


