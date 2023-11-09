""" run with 'sc2' conda env

"""
# import statements
import scanpy as sc

# function definitions
def get_nums(adata):
    agg = adata.obs.groupby('celltype').size()
    print(agg)

# setup
li_adata_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_dd_viz_ready_regs_updated.h5ad"
yost_adata_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/yost_bcc_full_processed_dd_regs_mod2.h5ad"
wang_adata_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/wang_dd_viz_ready_regs_cut_updated.h5ad"
li_cd8_adata_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_t_monocle_run4.h5ad'
yost_cd8_adata_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/yost_t_monocle_run4.h5ad'
wang_cd8_adata_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/wang_t_monocle_run4.h5ad'

# load the data
adata_li = sc.read_h5ad(li_adata_filename)
adata_yost = sc.read_h5ad(yost_adata_filename)
adata_wang = sc.read_h5ad(wang_adata_filename)
adata_li_t = sc.read_h5ad(li_cd8_adata_filename)
adata_yost_t = sc.read_h5ad(yost_cd8_adata_filename)
adata_wang_t = sc.read_h5ad(wang_cd8_adata_filename)

# get cell counts
#print(str(adata_li.shape))
#print(str(adata_yost.shape))
#print(str(adata_wang.shape))
#print(str(adata_li_t.shape))
#print(str(adata_yost_t.shape))
#print(str(adata_wang_t.shape))

# supplemental cell counts table -- cell counts by cell type, for each dataset
get_nums(adata_li)
get_nums(adata_yost)
get_nums(adata_wang)
