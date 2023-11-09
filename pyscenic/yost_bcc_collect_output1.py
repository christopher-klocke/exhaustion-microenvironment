"""
-collect results from pyscenic runs
-this script collects the pandas data frame containing regulon activity scores from pyscenic output loom and adding as 'adata.obsm['X_regulonsAUC']' field to adata object
 that was used as input into pyscenic run, and saving as .h5ad NOTE: this adata object might have slightly different pre-processing / cell numbers than main Li immune .h5ad
  file
-save-as of 'li_collect_output1.py'
"""
# import statements
import loompy as lp
import pandas as pd
import scanpy as sc


# file paths
run_dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/yost_bcc_run1/'
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'
yost_pyscenic_aucell_filename = 'aucell_output1.loom'
yost_bcc_main_filename = 'data/yost_data/yost_bcc_full_processed_dd.h5ad'  # from 'analysis2/yost_bcc_subset.py'
yost_bcc_out_filename = 'data/yost_data/yost_bcc_full_processed_dd_regs.h5ad'

# load the anndata object
adata_yost = sc.read_h5ad(dir + yost_bcc_main_filename)

# load the loom file
lf = lp.connect( (run_dir + yost_pyscenic_aucell_filename), mode='r+', validate=False)
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID) 
lf.close()

# just shared cells
a = list(auc_mtx.index)
b = list(adata_yost.obs.index)
adata_yost = adata_yost[[x in a for x in b]]
auc_mtx = auc_mtx.loc[[x in b for x in a]]

# fix order so they match
auc_mtx = auc_mtx.reindex(adata_yost.obs.index)

# add regulon activity scores to anndata object
adata_yost.obsm['X_regulonsAUC'] = auc_mtx

# write to file
adata_yost.write_h5ad(dir + yost_bcc_out_filename)
