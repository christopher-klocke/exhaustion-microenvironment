"""
-turn 'pyscenic_runs_prep.ipynb;' into .py script
-load in 3 .h5ad files, save output of 3 loom files (which will be used as input into pyscenic GRN inference algorithm)
-output goal -- create:
    -expression_converted_li.loom
    -expression_converted_wang.loom
    -Note: add Yost dataset later
-do not logarithmize or normalize, do not restrict to highly variable genes
-will need to go back to older adata objects for pyscenic runs
-no batch integration
-no filtering to just highly variable genes
-with high-dimensional genes (not cut to e.g. protein-coding), runtime will be much longer
"""

# import statements
import scanpy as sc
import loompy as lp
import numpy as np


# file paths
li_dir = '/Users/klockec/Documents/data/li_data/'
adata_li_pyscenic_filename = 'adata_li_pyscenic_input1.h5ad'
#yost_dir = '/Users/klockec/Documents/data/yost_data/' ## file path for rws07890
#yost_bcc_pyscenic_filename = 'yost_bcc_pyscenic_input1.h5ad'
#yost_scc_pyscenic_filename = 'yost_scc_pyscenic_input1.h5ad'
wang_dir = '/Users/klockec/Documents/data/wang_data/'
adata_wang_pyscenic_filename = 'adata_wang_pyscenic_input1.h5ad'
li_out_loom = 'li_pyscenic_input.loom'
#yost_out_loom = 'yost_pyscenic_input.loom'
wang_out_loom = 'wang_pyscenic_input.loom'

# load the data
adata_li = sc.read_h5ad(li_dir + adata_li_pyscenic_filename)
#adata_yost_bcc = sc.read_h5ad(yost_dir + yost_bcc_pyscenic_filename)
#adata_yost_scc = sc.read_h5ad(yost_dir + yost_scc_pyscenic_filename)
adata_wang = sc.read_h5ad(wang_dir + adata_wang_pyscenic_filename)

# re-remove sample 6 from wang data:
adata_wang_no6 = adata_wang[adata_wang.obs['batch'] != '6']

# remove un-needed attributes of .h5ad file ### CHECK THIS STEP ################
adata_li_reduced = adata_li
del adata_li_reduced.uns
del adata_li_reduced.obsm
del adata_li_reduced.varm
del adata_li_reduced.obsp

# YOST SECTION ADD ####
adata_wang_no6_reduced = adata_wang_no6
del adata_wang_no6_reduced.uns
del adata_wang_no6_reduced.obsm
del adata_wang_no6_reduced.varm
del adata_wang_no6_reduced.obsp

# create basic row and column attributes for the loom files:
row_attrs_li = {
    "Gene": np.array(adata_li_reduced.var_names) ,
}
col_attrs_li = {
    "CellID": np.array(adata_li_reduced.obs_names) ,
    "nGene": np.array( np.sum(adata_li_reduced.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata_li_reduced.X.transpose() , axis=0)).flatten() ,
}
# YOST SECTION ADD ####
row_attrs_wang = {
    "Gene": np.array(adata_wang_no6_reduced.var_names) ,
}
col_attrs_wang = {
    "CellID": np.array(adata_wang_no6_reduced.obs_names) ,
    "nGene": np.array( np.sum(adata_wang_no6_reduced.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata_wang_no6_reduced.X.transpose() , axis=0)).flatten() ,
}

# create loom files
lp.create( (li_dir + li_out_loom), adata_li_reduced.X.transpose(), row_attrs_li, col_attrs_li)
# YOST SECTION ADD ####
lp.create( (wang_dir + wang_out_loom), adata_wang_no6_reduced.X.transpose(), row_attrs_wang, col_attrs_wang)
