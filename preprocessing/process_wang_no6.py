"""remove sample 6"""
# import statements
import scanpy as sc

# file paths
wang_dir = '/Users/klockec/Documents/data/wang_data/'
wang_processed_filename = 'adata_wang_processed1.h5ad'
wang_no6_filename = 'adata_wang_no6.h5ad'

# load the data
adata_wang = sc.read_h5ad(wang_dir + wang_processed_filename)

# subset the data -- remove sample 6 cells
adata_wang_no6 = adata_wang[adata_wang.obs['batch'] != '6']

# write anndata object to .h5ad file
adata_wang_no6.write_h5ad(wang_dir + wang_no6_filename)
