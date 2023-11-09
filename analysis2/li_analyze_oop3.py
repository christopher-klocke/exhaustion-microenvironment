"""

script runs in ~30 seconds


NOTE: 

set 'threshold' by looking at histograms in, e.g., 'li_set_threshold.py' 

NOTE: 

replace 'patient' if called something different in dataset (e.g. 'batch')

NOTE: 

be sure to remove healthy controls with 'sample_id_list' list comprehension step 

NOTE: 

For Reactome pathway scores: 

lists_create(full_name=True)

gene_set_type =, e.g., 'X_pathwayAUC' 

NOTE:

set alpha to determine significance level for results 





save-as of 'li_analyze_oop2.py' 

"""

## import statements

import scanpy as sc
import os
from utils import lists_create, add_sample_scores, get_significant_results, run_score_correlation, plot_scores, plot_scores_v2, get_gene_set_results
import utils 

## setup 

plots_store = '/Users/klockec/Documents/data/analysis_files/p3/img_li_runA'

os.chdir(plots_store)

## file paths

li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890

li_filename = 'li_dd_viz_ready_regs.h5ad'

li_T_dpt_filename = 'li_T_dpt2.h5ad'

## setup

threshold = 0.05 ## ADJUST ###############

alpha = 0.001 ## LOWER THIS LATER -- FDR? ############

sample_name = 'patient' ## will be 'batch' for wang dataset

gene_set_type = 'X_regulonsAUC' ## could alternately use, e.g., Reactome pathway gene set (run AUCell on .gmt)

## load the anndata object

adata_li = sc.read_h5ad(li_path + li_filename)

adata = adata_li ## to make it easier to modify for other datasets

adata_li_T = sc.read_h5ad(li_path + li_T_dpt_filename)

adata_t = adata_li_T

## ADD CELL TYPES OF INTEREST 

cell_types_input = [('macrophage', ['3', '11']), ('NK', ['1']), ('B cells', ['5']), ('plasma', ['8']), ('ILC???', ['12']), ('CD8 T cells', ['0', '2', '4', '6', '9', '10']), ('???', ['7'])]

## REMOVE PBMC / HEALTHY CONTROL SAMPLES

sample_id_list = [x for x in set(list(adata_t.obs[sample_name])) if x[-4:] != 'PBMC'] ## need to adjust this for wang dataset ########
#sample_id_list = [x for x in set(list(adata_t.obs[sample_name])) if x != '0'] ##

##################################################################################

## GENERATE OBJECTS

cell_types, samples_list, activity_objects_list = lists_create(adata=adata, types=cell_types_input, sample_name=sample_name, gene_set=gene_set_type)

## GENERATE SCORES (sample-level exhaustion trajectory scores)

add_sample_scores(sample_ids=sample_id_list, samples=samples_list, adata=adata_t, sample_name=sample_name, threshold=threshold) ## should this be an object method instead? 




## ADD SECTION -- CELL COUNTS THRESHOLDING 

### FIRST -- PRINT COUNTS 

for x in activity_objects_list: 
    print('printing object...')
    print(str(x.celltype.type))
    print(str(x.sample.id))
    print(str(x.count))

### NEXT -- THRESHOLD ON COUNTS 

counts_threshold = 10 ## MINIMUM NUMBER OF CELLS OF SAMPLE IN CLUSTER IN ORDER TO CONSIDER ########## MODIFY 

passing_objects = [] ## objects that pass threshold 

for x in activity_objects_list: 
    if x.count > counts_threshold: 
        passing_objects.append(x)

## PRINT MORE DETAIL (BREAK DOWN BY CELL TYPE), BUT THIS IS FINE FOR NOW

print('full objects list: ' + str(len(activity_objects_list)))

print('passing objects: ' + str(len(passing_objects)))

### CAN NOW CONTINUE ANALYSIS, SUBSTITUTING 'passing_objects' for 'activity_objects_list'








## RUN CORRELATION ANALYSIS

#results_dict = run_score_correlation(celltypes=cell_types, samples = samples_list, activities = activity_objects_list) ##########
results_dict = run_score_correlation(celltypes=cell_types, samples = samples_list, activities = passing_objects) ## ONLY OBJECTS THAT PASS COUNTS THRESHOLD 

## IDENTIFY SIGNIFICANT RESULTS 

#get_significant_results(scores_dict=results_dict, alpha=alpha) #########

utils.results_to_file(scores_dict=results_dict, outfile='/Users/klockec/Documents/data/analysis_files/p3/mod1/li_results/results_perm1/dist_comparison_results.csv')





## SELECTED GENE SETS -- ACROSS CELL TYPES

#get_gene_set_results(scores_dict=results_dict, gene_set='E2F2')
#get_gene_set_results(scores_dict=results_dict, gene_set='E2F7')
#get_gene_set_results(scores_dict=results_dict, gene_set='E2F8')

#get_gene_set_results(scores_dict=results_dict, gene_set='NR2C1')
#get_gene_set_results(scores_dict=results_dict, gene_set='PATZ1')
#get_gene_set_results(scores_dict=results_dict, gene_set='SP2')



## VISUALIZE SIGNIFICANT RESULTS

"""
create function to visualize significant results

input will be results_dict, cell type, and gene set

will return scatterplot with sample_score values on the x axis and sample average gene set scores on y

"""





## VIOLIN PLOT GENERATE 

"""
may need to create new data structure -- 

try -- pandas df with cells as rows, column for patient id, column for dpt score 

this is just the adata.obs df -- use this 

"""


"""
import seaborn as sns

import matplotlib.pyplot as plt

sns.set_theme(style="whitegrid")

ax = sns.violinplot(x='patient', y='dpt_pseudotime', data=adata_t.obs)

# giving labels to x-axis and y-axis
ax.set(xlabel ='patient ID', ylabel ='CD8 T cell pseudotime')
 
# giving title to the plot
plt.title('Sample-Level Exhaustion Ordering')
 
# function to show plot
#plt.show(block=True) #############################

# function to save plot 
plt.savefig('li_violin1.pdf') ###################
"""





## new fig type -- generate 

#plot_scores_v2(results=results_dict, celltype='macrophage', gene_set='TCF7', samples=samples_list, activities=passing_objects, save='True', savefile='li_mac_TCF7_scatter_v2.pdf')

#plot_scores_v2(results=results_dict, celltype='macrophage', gene_set='LEF1', samples=samples_list, activities=passing_objects, save='True', savefile='li_mac_LEF1_scatter_v2.pdf')
#plot_scores_v2(results=results_dict, celltype='NK', gene_set='HIC1', samples=samples_list, activities=passing_objects, save='True', savefile='li_nk_HIC1_scatter_v2.pdf')



#plot_scores_v2(results=results_dict, celltype='CD8 T cells', gene_set='SP2', samples=samples_list, activities=passing_objects, save=True, savefile='cd8_sp2_scatter.png')

#plot_scores_v2(results=results_dict, celltype='CD8 T cells', gene_set='E2F2', samples=samples_list, activities=passing_objects, save=True, savefile='cd8_e2f2_scatter.png')
#plot_scores_v2(results=results_dict, celltype='CD8 T cells', gene_set='E2F7', samples=samples_list, activities=passing_objects, save=True, savefile='cd8_e2f7_scatter.png')
#plot_scores_v2(results=results_dict, celltype='CD8 T cells', gene_set='E2F8', samples=samples_list, activities=passing_objects, save=True, savefile='cd8_e2f8_scatter.png')

#plot_scores_v2(results=results_dict, celltype='macrophage', gene_set='NR2C1', samples=samples_list, activities=passing_objects, save=True, savefile='mac_nr2c1_scatter.png')
#plot_scores_v2(results=results_dict, celltype='macrophage', gene_set='PATZ1', samples=samples_list, activities=passing_objects, save=True, savefile='mac_patz1_scatter.png')
#plot_scores_v2(results=results_dict, celltype='macrophage', gene_set='E2F7', samples=samples_list, activities=passing_objects, save=True, savefile='mac_e2f7_scatter.png')

#plot_scores_v2(results=results_dict, celltype='B cells', gene_set='E2F2', samples=samples_list, activities=passing_objects, save=True, savefile='b_e2f2_scatter.png')
#plot_scores_v2(results=results_dict, celltype='B cells', gene_set='E2F7', samples=samples_list, activities=passing_objects, save=True, savefile='b_e2f7_scatter.png')









"""
COMMENT OUT FOR NOW -- #####################

plot_scores(results=results_dict, celltype='macrophage', gene_set='FOXP1', samples = samples_list, activities = activity_objects_list, save=True, savefile='mac_foxp1_scatter.png')
plot_scores(results=results_dict, celltype='macrophage', gene_set='KLF17', samples = samples_list, activities = activity_objects_list, save=True, savefile='mac_klf17_scatter.png')
plot_scores(results=results_dict, celltype='macrophage', gene_set='LEF1', samples = samples_list, activities = activity_objects_list, save=True, savefile='mac_lef1_scatter.png')
plot_scores(results=results_dict, celltype='macrophage', gene_set='NFIB', samples = samples_list, activities = activity_objects_list, save=True, savefile='mac_nfib_scatter.png')
plot_scores(results=results_dict, celltype='macrophage', gene_set='POU6F1', samples = samples_list, activities = activity_objects_list, save=True, savefile='mac_pou6f1_scatter.png')
plot_scores(results=results_dict, celltype='macrophage', gene_set='REST', samples = samples_list, activities = activity_objects_list, save=True, savefile='mac_rest_scatter.png')
plot_scores(results=results_dict, celltype='macrophage', gene_set='TCF7', samples = samples_list, activities = activity_objects_list, save=True, savefile='mac_tcf7_scatter.png')
plot_scores(results=results_dict, celltype='macrophage', gene_set='ZNF224', samples = samples_list, activities = activity_objects_list, save=True, savefile='mac_znf224_scatter.png')
plot_scores(results=results_dict, celltype='macrophage', gene_set='ZNF740', samples = samples_list, activities = activity_objects_list, save=True, savefile='mac_znf740_scatter.png')

plot_scores(results=results_dict, celltype='NK', gene_set='HES4', samples = samples_list, activities = activity_objects_list, save=True, savefile='nk_hes4_scatter.png')
plot_scores(results=results_dict, celltype='NK', gene_set='HIC1', samples = samples_list, activities = activity_objects_list, save=True, savefile='nk_hic1_scatter.png')
plot_scores(results=results_dict, celltype='NK', gene_set='HIVEP3', samples = samples_list, activities = activity_objects_list, save=True, savefile='nk_hivep3_scatter.png')
plot_scores(results=results_dict, celltype='NK', gene_set='HLF', samples = samples_list, activities = activity_objects_list, save=True, savefile='nk_hlf_scatter.png')
plot_scores(results=results_dict, celltype='NK', gene_set='NR5A2', samples = samples_list, activities = activity_objects_list, save=True, savefile='nk_nr5a2_scatter.png')
plot_scores(results=results_dict, celltype='NK', gene_set='PKNOX1', samples = samples_list, activities = activity_objects_list, save=True, savefile='nk_pknox1_scatter.png')
plot_scores(results=results_dict, celltype='NK', gene_set='RFX3', samples = samples_list, activities = activity_objects_list, save=True, savefile='nk_rfx3_scatter.png')
plot_scores(results=results_dict, celltype='NK', gene_set='ZNF366', samples = samples_list, activities = activity_objects_list, save=True, savefile='nk_znf366_scatter.png')

plot_scores(results=results_dict, celltype='plasma', gene_set='E2F7', samples = samples_list, activities = activity_objects_list, save=True, savefile='plasma_e2f7_scatter.png')
plot_scores(results=results_dict, celltype='plasma', gene_set='HOXD9', samples = samples_list, activities = activity_objects_list, save=True, savefile='plasma_hoxd9_scatter.png')
plot_scores(results=results_dict, celltype='plasma', gene_set='JDP2', samples = samples_list, activities = activity_objects_list, save=True, savefile='plasma_jdp2_scatter.png')
plot_scores(results=results_dict, celltype='plasma', gene_set='TBX6', samples = samples_list, activities = activity_objects_list, save=True, savefile='plasma_tbx6_scatter.png')
plot_scores(results=results_dict, celltype='plasma', gene_set='TP73', samples = samples_list, activities = activity_objects_list, save=True, savefile='plasma_tp73_scatter.png')
plot_scores(results=results_dict, celltype='plasma', gene_set='ZBTB33', samples = samples_list, activities = activity_objects_list, save=True, savefile='plasma_zbtb33_scatter.png')
plot_scores(results=results_dict, celltype='plasma', gene_set='ZNF76', samples = samples_list, activities = activity_objects_list, save=True, savefile='plasma_znf76_scatter.png')

plot_scores(results=results_dict, celltype='ILC???', gene_set='GLI3', samples = samples_list, activities = activity_objects_list, save=True, savefile='ilc_gli3_scatter.png')
plot_scores(results=results_dict, celltype='ILC???', gene_set='HOXC11', samples = samples_list, activities = activity_objects_list, save=True, savefile='ilc_hoxc11_scatter.png')
plot_scores(results=results_dict, celltype='ILC???', gene_set='ZNF267', samples = samples_list, activities = activity_objects_list, save=True, savefile='ilc_znf267_scatter.png')

"""