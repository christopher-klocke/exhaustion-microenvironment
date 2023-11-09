"""run with 'sc2_magic' conda env

perform imputation with the MAGIC algorithm
"""

# import statements
import scanpy as sc
import magic

# setup
input_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_t_monocle_run4.h5ad'  # from 'analysis2/monocle3_r/run4_collect_li.py'
output_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_t_monocle_run4_imputed.h5ad'
gene_list = ['TCF7', 'TOX', 'CD8A', 'CD8B', 'PDCD1', 'LAG3', 'TIGIT']

# load the data
adata = sc.read_h5ad(input_filename)

# impute with MAGIC algorithm
magic_operator = magic.MAGIC()
adata_magic = magic_operator.fit_transform(adata, genes=gene_list)
adata_magic.uns['umap'] = adata.uns['umap']
adata_magic.obsm['X_umap'] = adata.obsm['X_umap']

# write to output file
adata_magic.write_h5ad(output_filename)
