# import statements
import scanpy as sc


# file paths
wang_dir = '/Users/klockec/Documents/data/wang_data/'
wang_pc_filename = 'wang_adata_pc.h5ad'
adata_wang_pyscenic_filename = 'adata_wang_pyscenic_input1.h5ad'
wang_processed_filename = 'adata_wang_processed1.h5ad'

# load the data
adata_wang = sc.read_h5ad(wang_dir + wang_pc_filename)

# process the data
adata_wang.var_names_make_unique()
sc.pp.filter_cells(adata_wang, min_genes=200)
sc.pp.filter_genes(adata_wang, min_cells=3)
adata_wang.var['mt'] = adata_wang.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_wang, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_wang = adata_wang[adata_wang.obs.n_genes_by_counts < 2500, :]
adata_wang = adata_wang[adata_wang.obs.pct_counts_mt < 10, :]

# save for input into pyscenic
adata_wang.write_h5ad(wang_dir + adata_wang_pyscenic_filename)

# continue pre-processing
sc.pp.normalize_total(adata_wang, target_sum=1e4)
sc.pp.log1p(adata_wang)
sc.pp.highly_variable_genes(adata_wang)
adata_wang.raw = adata_wang

# write anndata object to .h5ad file
adata_wang.write_h5ad(wang_dir + wang_processed_filename)
