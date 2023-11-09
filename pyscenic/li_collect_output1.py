"""-collect results from pyscenic runs
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
run_dir = '/Users/klockec/Documents/data/analysis_files/p3/pyscenic_runs/p3/li_run/'
li_pyscenic_aucell_filename = 'aucell_output1.loom'
li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
li_main_filename = 'li_dd_viz_ready.h5ad'
li_out_filename = 'li_dd_viz_ready_regs.h5ad'

# load the anndata object
adata_li = sc.read_h5ad(li_path + li_main_filename)

# load the loom file
lf = lp.connect( (run_dir + li_pyscenic_aucell_filename), mode='r+', validate=False)
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID) 
lf.close()

# just shared cells
a = list(auc_mtx.index)
b = list(adata_li.obs.index)
adata_li = adata_li[[x in a for x in b]]
auc_mtx = auc_mtx.loc[[x in b for x in a]]

# fix order so they match
auc_mtx = auc_mtx.reindex(adata_li.obs.index)

# add regulon activity scores to anndata object
adata_li.obsm['X_regulonsAUC'] = auc_mtx

# write to file
adata_li.write_h5ad(li_path + li_out_filename)
