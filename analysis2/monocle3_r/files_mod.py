# import statements
import scanpy as sc

# setup
li_new_filename = '/Users/klockec/Documents/data/li_data/li_t_modified.h5ad'
wang_new_filename = '/Users/klockec/Documents/data/wang_data/wang_t_modified.h5ad'

# load the data
adata_li = sc.read_h5ad('/Users/klockec/Documents/data/li_data/li_T_refined.h5ad')
adata_wang = sc.read_h5ad('/Users/klockec/Documents/data/wang_data/wang_dd_T_viz_ready.h5ad')

# process the data
del adata_li.uns
del adata_li.obsm
del adata_li.varm
del adata_li.obsp
del adata_wang.uns
del adata_wang.obsm
del adata_wang.varm
del adata_wang.obsp

# write to output files
adata_li.write_h5ad(li_new_filename)
adata_wang.write_h5ad(wang_new_filename)
