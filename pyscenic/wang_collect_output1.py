"""
-collect results from pyscenic runs
-this script collects the pandas data frame containing regulon activity scores from pyscenic output loom and adding as 'adata.obsm['X_regulonsAUC']' field to adata object
 that was used as input into pyscenic run, and saving as .h5ad NOTE: this adata object might have slightly different pre-processing / cell numbers than main Li immune .h5ad
  file
"""
# import statements
import loompy as lp
import pandas as pd
import scanpy as sc
#import regulon_viz as rv
#import importlib.util


# file paths
run_dir = '/Users/klockec/Documents/data/analysis_files/p3/pyscenic_runs/p3/wang_run/'
wang_pyscenic_aucell_filename = 'aucell_output1.loom'
wang_path = '/Users/klockec/Documents/data/wang_data/' ## file path for rws07890
wang_main_filename = 'wang_dd_viz_ready.h5ad'
wang_out_filename = 'wang_dd_viz_ready_regs.h5ad'

# load the anndata object
adata_wang = sc.read_h5ad(wang_path + wang_main_filename)

# load the loom file
lf = lp.connect( (run_dir + wang_pyscenic_aucell_filename), mode='r+', validate=False)
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID) 
lf.close()

# just shared cells
a = list(auc_mtx.index)
b = list(adata_wang.obs.index)
adata_li = adata_wang[[x in a for x in b]]
auc_mtx = auc_mtx.loc[[x in b for x in a]]

# fix order so they match
auc_mtx = auc_mtx.reindex(adata_wang.obs.index)

# add regulon activity scores to anndata object
adata_wang.obsm['X_regulonsAUC'] = auc_mtx

# write to file
adata_wang.write_h5ad(wang_path + wang_out_filename)
