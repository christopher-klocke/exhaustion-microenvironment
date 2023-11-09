"""
-save-as of 'pyscenic_runs_prep.py'
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
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'
yost_bcc_pyscenic_filename = 'data/yost_data/yost_bcc_pyscenic_input1.h5ad'
yost_scc_pyscenic_filename = 'data/yost_data/yost_scc_pyscenic_input1.h5ad'
yost_bcc_out_loom = 'data/yost_data/yost_bcc_pyscenic_input.loom'
yost_scc_out_loom = 'data/yost_data/yost_scc_pyscenic_input.loom'

# load the data
adata_yost_bcc = sc.read_h5ad(dir + yost_bcc_pyscenic_filename)
adata_yost_scc = sc.read_h5ad(dir + yost_scc_pyscenic_filename)

# remove un-needed attributes of .h5ad file
adata_yost_bcc_reduced = adata_yost_bcc
del adata_yost_bcc_reduced.uns
del adata_yost_bcc_reduced.obsm
del adata_yost_bcc_reduced.varm
del adata_yost_bcc_reduced.obsp

adata_yost_scc_reduced = adata_yost_scc
del adata_yost_scc_reduced.uns
del adata_yost_scc_reduced.obsm
del adata_yost_scc_reduced.varm
del adata_yost_scc_reduced.obsp

# create basic row and column attributes for the loom files:
row_attrs_yost_bcc = {
    "Gene": np.array(adata_yost_bcc_reduced.var_names) ,
}
col_attrs_yost_bcc = {
    "CellID": np.array(adata_yost_bcc_reduced.obs_names) ,
    "nGene": np.array( np.sum(adata_yost_bcc_reduced.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata_yost_bcc_reduced.X.transpose() , axis=0)).flatten() ,
}

row_attrs_yost_scc = {
    "Gene": np.array(adata_yost_scc_reduced.var_names) ,
}
col_attrs_yost_scc = {
    "CellID": np.array(adata_yost_scc_reduced.obs_names) ,
    "nGene": np.array( np.sum(adata_yost_scc_reduced.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata_yost_scc_reduced.X.transpose() , axis=0)).flatten() ,
}

# create loom files
lp.create( (dir + yost_bcc_out_loom), adata_yost_bcc_reduced.X.transpose(), row_attrs_yost_bcc, col_attrs_yost_bcc)
lp.create( (dir + yost_scc_out_loom), adata_yost_scc_reduced.X.transpose(), row_attrs_yost_scc, col_attrs_yost_scc)
