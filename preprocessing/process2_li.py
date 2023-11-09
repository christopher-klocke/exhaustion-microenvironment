# import statements
import scanpy as sc


# file paths
li_out_dir_name = '/Users/klockec/Documents/data/li_data/'
adata_li_subsetted_pc_filename = 'adata_li_subsetted_pc.h5ad'
adata_li_pyscenic_filename = 'adata_li_pyscenic_input1.h5ad'
li_processed_filename = 'li_processed.h5ad'

# load the data
adata_li_subsetted_pc = sc.read_h5ad(li_out_dir_name + adata_li_subsetted_pc_filename)

# continue pre-processing of li dataset
sc.pp.filter_cells(adata_li_subsetted_pc, min_genes=200)
sc.pp.filter_genes(adata_li_subsetted_pc, min_cells=3)
adata_li_subsetted_pc.var['mt'] = adata_li_subsetted_pc.var_names.str.startswith('MT-') | adata_li_subsetted_pc.var_names.str.startswith('mt-') ## NEED TO MOD FOR DIFFERENT MT GENE NOTATION
sc.pp.calculate_qc_metrics(adata_li_subsetted_pc, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_li_subsetted_pc = adata_li_subsetted_pc[adata_li_subsetted_pc.obs.n_genes_by_counts < 5000, :]

# save for input into pyscenic
adata_li_subsetted_pc.write_h5ad(li_out_dir_name + adata_li_pyscenic_filename)

# continue pre-processing
sc.pp.normalize_total(adata_li_subsetted_pc, target_sum=1e4)
sc.pp.log1p(adata_li_subsetted_pc)
sc.pp.highly_variable_genes(adata_li_subsetted_pc)
adata_li_subsetted_pc.raw = adata_li_subsetted_pc

# write anndata object to .h5ad file
adata_li_subsetted_pc.write_h5ad(li_out_dir_name + li_processed_filename)
