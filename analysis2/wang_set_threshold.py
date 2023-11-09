"""


"""

## import statements

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os 

## setup 

sample_name = 'batch' ## will be 'batch' for wang dataset

plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img_wang_runA'

os.chdir(plots_store)

## file paths

wang_dir = '/Users/klockec/Documents/data/wang_data/'

wang_T_dpt_filename = 'wang_T_dpt2.h5ad'

## load the data

adata_wang_T = sc.read_h5ad(wang_dir + wang_T_dpt_filename)

adata = adata_wang_T

## process the data 

### remove PBMC samples
sample_id_list = [x for x in set(list(adata.obs[sample_name])) if x != '0'] ## need to adjust this for wang dataset ########

adata_subsets = [] #####

for sample_id in sample_id_list: 
    ## create patient-specific anndata objects -- subsetted by patient id 
    ## anndata object, subsetted by given sample; selecting dpt scores from adata.obs
    ## result is list of scores 
    adata_subsetted = adata[adata.obs[sample_name] == sample_id].obs['dpt_pseudotime']
    adata_subsets.append(adata_subsetted)


# HISTOGRAMS

"""
## create patient-specific anndata objects -- subsetted by patient id 

## anndata object, subsetted by given sample; selecting dpt scores from adata.obs
## result is list of scores 

p0 = adata_wang_t[adata_wang_t.obs['batch'] == '0'].obs['dpt_pseudotime']
p1 = adata_wang_t[adata_wang_t.obs['batch'] == '1'].obs['dpt_pseudotime']
p2 = adata_wang_t[adata_wang_t.obs['batch'] == '2'].obs['dpt_pseudotime']
p3 = adata_wang_t[adata_wang_t.obs['batch'] == '3'].obs['dpt_pseudotime']
p4 = adata_wang_t[adata_wang_t.obs['batch'] == '4'].obs['dpt_pseudotime']
p5 = adata_wang_t[adata_wang_t.obs['batch'] == '5'].obs['dpt_pseudotime']

"""

## PLOT HISTOGRAMS, SELECT THRESHOLD

threshold = 0.42

n_bins = 100 ##########

fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)

all_points = adata.obs['dpt_pseudotime']




"""
# We can set the number of bins with the *bins* keyword argument.
#axs[0].hist(all_points, bins=n_bins)
#axs[1].hist(patient_n_dpt, bins=n_bins)

plt.axvline(x=threshold, color='red') ## plot threshold line 
#axs.hist(all_points) ## plot histogram with cells from all samples (must be filtered to just T cells already)
axs.hist(all_points, bins=n_bins)

plt.show(block=True)


# function to save plot 
plt.savefig('wang_cd8_hist.pdf')



fig, axs = plt.subplots(1, 4, sharey=True, tight_layout=False)
axs[0].hist(adata_subsets[0]) ## plot histograms for individual samples -- a few at a time, otherwise x axis is too crowded 
axs[1].hist(adata_subsets[1])
axs[2].hist(adata_subsets[2])
axs[3].hist(adata_subsets[3])

#plt.show(block=True)
"""





"""

generate joint histogram / density plot 

"""

sns.distplot(all_points, bins=80)

plt.axvline(x=threshold, color='red') ## plot threshold line 

plt.savefig('wang_cd8_dens_v1.pdf')
#plt.show(block=True)

