""" run with 'sc2' conda env """

import logging 

logging.basicConfig(filename='dist_median_monocle_wang2.log', filemode='w', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logging.info('run start')

import utils
import scanpy as sc

wang_path = '/Users/klockec/Documents/data/wang_data/' ## file path for rws07890
wang_filename = 'wang_dd_viz_ready_regs_cut.h5ad' ## updated -- clusters 9 and 10 removed
wang_pt_input_filename = 'wang_t_monocle1.h5ad' ## from 'p3/analysis2/monocle3_r/run2_collect_wang.py' 
results_dir = '/Users/klockec/Documents/data/analysis_files/p3/mod1/wang_dist2_monocle_median/'

cell_types_input = [('macrophage', ['4', '8']), ('plasma', ['7']), ('B', ['1']), ('CD8_T', ['0', '2', '3', '5', '6'])]
#cell_types_input = {'macrophage': ['4', '8'], 'plasma': ['7'], 'B': ['1'], 'CD8_T': ['0', '2', '3', '5', '6']}

adata_wang_pt = sc.read_h5ad(wang_path + wang_pt_input_filename) 
adata_wang_gene_sets = sc.read_h5ad(wang_path + wang_filename) 

sample_name = 'batch'
gene_set_type = 'X_regulonsAUC'
threshold = 0.05 ## ADJUST ###############

sample_id_list = [x for x in set(list(adata_wang_pt.obs[sample_name]))]

cell_types, samples_list, activity_objects_list = utils.lists_create(adata=adata_wang_gene_sets, types=cell_types_input, sample_name=sample_name, gene_set=gene_set_type, mode='median')

utils.add_sample_scores(sample_ids=sample_id_list, samples=samples_list, adata=adata_wang_pt, sample_name=sample_name, threshold=threshold, pt_name='monocle3_pseudotime') ## should this be an object method instead? 

counts_threshold = 10 ## MINIMUM NUMBER OF CELLS OF SAMPLE IN CLUSTER IN ORDER TO CONSIDER ########## MODIFY 

passing_objects = [] ## objects that pass threshold 

for x in activity_objects_list: 
    if x.count > counts_threshold: 
        passing_objects.append(x)

results_dict = utils.run_score_correlation(celltypes=cell_types, samples=samples_list, activities=passing_objects, permute=True, n_perm=1000)

utils.results_to_file(scores_dict=results_dict, outfile='/Users/klockec/Documents/data/analysis_files/p3/mod1/wang_run2_dist_monocle_median.csv')

#utils.plot_scores_v2(results=results_dict, celltype='B', gene_set='FOXO3', samples=samples_list, activities=activity_objects_list, save=True, savefile='/Users/klockec/Documents/data/analysis_files/p3/mod1/img_wang2/dist_monocle_median_B_FOXO3.png')
#utils.plot_scores_v2(results=results_dict, celltype='B', gene_set='STAT5A', samples=samples_list, activities=activity_objects_list, save=True, savefile='/Users/klockec/Documents/data/analysis_files/p3/mod1/img_wang2/dist_monocle_median_B_STAT5A.png')


#utils.plot_scores_v2(results=results_dict, celltype='macrophage', gene_set='AEBP2', samples=samples_list, activities=activity_objects_list, save=True, savefile='/Users/klockec/Documents/data/analysis_files/p3/mod1/img_wang2/dist_monocle_median_mac_AEBP2.png')
#utils.plot_scores_v2(results=results_dict, celltype='macrophage', gene_set='MEIS1', samples=samples_list, activities=activity_objects_list, save=True, savefile='/Users/klockec/Documents/data/analysis_files/p3/mod1/img_wang2/dist_monocle_median_mac_MEIS1.png')
#utils.plot_scores_v2(results=results_dict, celltype='macrophage', gene_set='POU5F1', samples=samples_list, activities=activity_objects_list, save=True, savefile='/Users/klockec/Documents/data/analysis_files/p3/mod1/img_wang2/dist_monocle_median_mac_POU5F1.png')
#utils.plot_scores_v2(results=results_dict, celltype='macrophage', gene_set='NFE2L2', samples=samples_list, activities=activity_objects_list, save=True, savefile='/Users/klockec/Documents/data/analysis_files/p3/mod1/img_wang2/dist_monocle_median_mac_NFE2L2.png')





import sys
sys.exit()

ceu.csea_wrapper(
    adata_pt_filename=(wang_path + wang_pt_input_filename), 
    adata_act_filename=(wang_path + wang_filename), 
    results_dir=results_dir, 
    celltype_dict=cell_types_input, 
    sample_id='batch', 
    set_n_perm=1000, 
    alpha=0.2
    )
