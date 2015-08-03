"""
Created by Emille Ishida in July, 2015

This script contains the necessary steps to perform the dimensionality
reduction, cross-validation and photometric supernova classification
using snclass.

It assumes all the data were previously fitted with GP as explained in
snclass README file. 

IMPORTANT: run this sequence in a dedicated directory, since it will
generate a large number of intermediate files. 
"""
from snclass.algorithm import build_sample, set_lclist, set_parameters
from snclass.algorithm import sample_pop, photo_frac, get_names, select_GP
from snclass.algorithm import build_spec_matrix, plot_proj, read_matrix
from snclass.algorithm import read_hyperpar, set_kpca_obj, classify
from snclass.functions import screen
from snclass.util import read_user_input

##############################################################################
### User choices

input_file = 'user.input'       	# user input file

epoch_min = 'p5'                	# minimum epoch for name purpouses 
				                	# if minimum epoch > 0, 
			                		# use a 'p' in front
			                   		# of the number
					
epoch_max = '45'               		# maximum epoch for name purpouses
ref_filter_name = 'None'          	# reference filter, where peakMJD will
			                		# be defined
					
filters = 'iz'                		# sequence of all filters

# determine which type of spectroscopic sample to use in the classification
# in representative and balanced options, the extra objects needed to complete
# the sample are taken from realizations of the final GP fit
#
# ***WARNING***
# *** Always run first an original version and aftwewards the other options!
#     This will create the necessary auxiliary files.  ***
#
# options are:
# original -> use data as it is
# representative -> build a representative spectroscopic sample using the 
# 		    same proportions as those present in the photometric 
#		    sample (this can only be used with simulations)
# balanced -> construct a spectroscopic sample where all types of SNe are
#	      equally represented
representation = 'representative'

# size of spectroscopic sample
# as a first try,
# user the size of the spectroscopic sample before any selection cuts 
sample_size = 1103

# directory to store lc plots plots
# if None no plots are generated
plot_dir = None

# range of number of PCs to be tested
range_pcs = [2,36]


### File names and directories

# root directory where all results from this SNR cut will be stored
SNR_dir = '~/my_dir/'

# directory to store lists for different samples
list_dir = SNR_dir + 'lists/'

# directory with GP fitted data from spec sample
fitted_spec_data_dir = '~/fitted_data/spec/'

# directory with GP fitted data from photo sample
fitted_photo_data_dir = '~/fitted_data/photo/'

# path to training sample list file
spec_list_name = list_dir + 'spec_' + filters + '_' + epoch_min + \
                    '_' + epoch_max + '_ref_' + ref_filter_name + '.list'

# path to test sample list file
photo_list_name = list_dir + 'photo_' + filters + '_' + epoch_min + \
                  '_' + epoch_max + '_ref_' + ref_filter_name + '.list'


# directory to store GP results for the syntheticaly constructed spec sample
synthetic_dir = SNR_dir + 'fitted_data/' + representation +'_spec_' + \
                filters + '_' + epoch_min + '_' + epoch_max + '_ref_' + \
                ref_filter_name + '/'

# directory to store classification results
out_dir = SNR_dir + 'class_results/' + representation + '_' + filters + \
          '_' + epoch_min + '_' + epoch_max + '_ref_' + ref_filter_name +'/'

# directory to store matrices
mat_dir = SNR_dir + 'matrices/'

# spec data matrix file
data_matrix = mat_dir + representation + '_' + filters + '_' + epoch_min + \
              '_' + epoch_max + '_ref_' + ref_filter_name + '_data_matrix.dat'

# directory to store projections plot
plot_proj_dir = SNR_dir + 'plots/proj/' + representation + '_' + filters + \
                '_' + epoch_min + '_' + epoch_max + '_ref_' + \
                ref_filter_name + '/'

# training sample to be scanned
sample1 = 'spec'

# test sample to be scanned
sample2 = 'photo'

# define types (types and corresponding number codes)
# the values bellow hold for post-SNPCC data only
type_number = {}
type_number['Ibc'] = ['1', '5', '6', '7', '8', '9', '10', '11', '13', '14',
                      '16','18', '22', '23', '28', '29', '45']
type_number['II'] = ['2','3','4','12','15','17','19','20','21','24', '25',
                     '26', '27', '30','31','32','33','34','35','36','37',
                     '38','39','40','41','42', '43','44','46']
type_number['Ia'] = ['0'] 


##############################################################################

### Build input dictionary parameters
p = {}
p['input_file'] = input_file
p['representation'] = representation
p['sample_size'] = sample_size
p['fname_photo_list'] = photo_list_name
p['data_matrix'] = data_matrix

p['plot_dir'] = plot_dir
p['list_dir'] = list_dir
p['mat_dir'] = mat_dir
p['out_dir'] = out_dir
p['SNR_dir'] = SNR_dir
p['plot_dir'] = plot_dir
p['plot_proj_dir'] = plot_proj_dir
p['synthetic_dir'] = synthetic_dir

p['range_pcs'] = range_pcs

p['sample'] = sample1
p['fitted_data_dir'] = fitted_spec_data_dir

### Read user input choices
p['user_choices'] = read_user_input(input_file)

### Build lists of SN for the chosen representation

screen('Build lists', p['user_choices'])

# set first object list
if p['representation'] == 'original':
    set_lclist(p)

# set new parameters for second sample scan
p['sample'] = sample2
p['fitted_data_dir'] = fitted_photo_data_dir

# set second object list
if p['representation'] == 'original':
    set_lclist(p)

### Set parameters for the synthetic spectroscopic sample
screen('Set parameters for the synthetic spec sample', p['user_choices'])

# count spec sample population of each type
p['list_name'] = spec_list_name
p['spec_pop'] = sample_pop(p['user_choices'], p, type_number)

# count photo sample population type
p['list_name'] = photo_list_name
p['photo_pop'] = sample_pop(p['user_choices'], p, type_number)

# determine fractions in photo sample
p['photo_perc'] = photo_frac(p['spec_pop'], p['photo_pop'],
                             p['representation'])

# get object identification for spec sample
p['list_name'] = spec_list_name
p['surv_spec_names'] = get_names(p['user_choices'], p, type_number)

# set extra parameters
p['fitted_data_dir'] = fitted_spec_data_dir
p = set_parameters(p)

# build synthetic sample
screen('Build synthetic sample', p['user_choices'])
select_GP(p, p['user_choices'])

### Build spectroscopic data matrix

screen('Build spectroscopic data matrix', p['user_choices'])
build_spec_matrix(p, type_number)

### Classify photometric sample
screen('Classify photometric sample', p['user_choices'])
p['photo_dir'] = fitted_photo_data_dir
classify(p, p['user_choices'], type_number, do_plot=False)

