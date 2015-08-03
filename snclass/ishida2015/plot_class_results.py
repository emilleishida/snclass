"""
Created by Emille Ishida on July, 2015

This script contains the necessary steps to generate diagnostic plots from
classification results obtained with snclass code.

IMPORTANT: run this script in a dedicated directory, since it will generate
a large number of directories and files with intermediate results!
"""
from snclass.diagnostic_plots import calc_ROC, det_threshold, calc_binned_diag
from snclass.diagnostic_plots import calc_global_diag, count_pop
from snclass.diagnostic_plots import read_literature_results
from snclass.diagnostic_plots import plot_diagnostics_aftercuts
from snclass.diagnostic_plots import plot_diagnostics_before_cuts
from snclass.algorithm import read_file

##############################################################################
### User choices

# user input file
user_input = 'user.input'

# epoch cuts
epoch_min = 'p5'		# minimum epoch for name purpouses 
				# if minimum epoch > 0, 
				# use a 'p' in front
				# of the number
				
epoch_max = '45'		# maximum epoch for name purpouses

filters = 'iz'			# complete filter list

SNR_cut = '5'			# quality cut

ref_filter_name = 'None'	# reference filter
ref_filter = ref_filter_name

representation = 'balanced'	# spec sample type

screen = '1'			# control screen output
				# if 0 no  intermediate steps are 
				# shown on screen

# range of number of PCs for which cross-validation results exist
range_pcs = [2, 36]

# width of redshift bin for plotting diagnostic results
dz = 0.2

# number of bins for plotting diagnostic results
# PS: use a large number (>5) this will be adjust up to where there are
# data satisfying selection cuts
nbins = 8 

### Paths and file names

# root directory where all results from this SNR cut will be stored
SNR_dir = '~/mydir/'

# directory holding classification results
class_res_dir = SNR_dir + 'class_results/' + representation + '_' + \
                filters + '_' + epoch_min + '_' + epoch_max + '_ref_' + \
                ref_filter_name +'/'


# define types (the values bellow hold for post-SNPCC data only)
type_number = {}
type_number['Ibc'] = ['1', '5', '6', '7', '8', '9', '10', '11', '13', '14',
                      '16','18', '22', '23', '28', '29', '45']
type_number['II'] = ['2','3','4','12','15','17','19','20','21','24', '25',
                     '26', '27', '30','31','32','33','34','35','36','37',
                     '38','39','40','41','42', '43','44','46']
type_number['Ia'] = ['0']


# If you do not want to plot or do not have results from the literature in
# the correct form, set these to None
# 
# Remember that there is no point in comparing literature results data 
# which is not the post-SNPCC

# path to Ishida & de Souza, 2013 results
path_Spl = 'ishida_2013_cad1_SNR5.dat'

# path to Karpenka et al., 2012 results
path_karp = 'karpenka_'

# path to Richards et al., 2011 results
path_rich = 'richards_'

##############################################################################

### Build input dictionary parameters
p = {}
p['range_pcs'] = range_pcs
p['ref_filter'] = ref_filter
p['screen'] = screen
p['representation'] = representation
p['epoch_min'] = epoch_min
p['epoch_max'] = epoch_max
p['SNR_cut'] = SNR_cut
p['filters'] = filters
p['dz'] = dz
p['nbins'] = nbins
p['user_input'] = user_input

p['class_res_dir'] = class_res_dir
p['path_karp'] = path_karp
p['path_rich'] = path_rich

if len(p['filters']) == 4:
    p['path_Spl'] = path_Spl
else:
    p['path_Spl'] = None

# calculate ROC curve for a range of number of PCs
p = calc_ROC(p, plot=True)

# determine optimized result and store on file
det_threshold(p)

# read configuration file
data1 = read_file(p['class_res_dir'] + 'ROC_results.dat')

p['config'] = {}
for line in data1:
    p['config'][line[0]] = float(line[1])

# read classification results
data2 = read_file(p['class_res_dir'] + \
                  str(int(p['config']['ncomp:'])) + 'PC/class_res_' + \
                  str(int(p['config']['ncomp:'])) + 'PC.dat')

# feed classification results to dictionary of parameters
p['class_res'] = {}
for line in data2[1:]:
    p['class_res'][line[0]] = [line[1], float(line[2])]

# count population of each type
p = count_pop(p, type_number)

# calculate diagnostics for entire redshift range 
p = calc_global_diag(p)

# calculate diagnostics in redshift bins
p = calc_binned_diag(p)

# read results from the literature
p = read_literature_results(p)

# plot diagnostics considering only objects which survived selection cuts
plot_diagnostics_aftercuts(p)

# plot diagnostics considering all objects from original sample
plot_diagnostics_before_cuts(p)
