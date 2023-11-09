"""

"""
#  import statements
import scanpy as sc
import magic

#  setup
dir_name = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'
li_filename = 'data/li_data/li_t_monocle1.h5ad'
yost_filename = 'data/yost_data/yost_t_monocle1_v2.h5ad'
wang_filename = 'data/wang_data/wang_t_monocle1.h5ad'
li_imputed_filename = 'data/li_data/li_t_monocle1_imputed.h5ad'
yost_imputed_filename = 'data/yost_data/yost_t_monocle1_v2_imputed.h5ad'
wang_imputed_filename = 'data/wang_data/wang_t_monocle1_imputed.h5ad'
genes1 = ['TCF7', 'TOX', 'LAG3', 'CD8B', 'PDCD1', 'TIGIT']

#  load the data
adata_li = sc.read_h5ad(dir_name + li_filename)
adata_yost = sc.read_h5ad(dir_name + yost_filename)
adata_wang = sc.read_h5ad(dir_name + wang_filename)

#  function definition
def impute_adata(
        adata,
        output_filename: str,
        impute_genes: list
):
    """
    impute the anndata file with the given gene list
    """
    X_adata = adata
    magic_operator = magic.MAGIC()
    X_magic = magic_operator.fit_transform(X_adata, genes=impute_genes)
    X_magic.uns['umap'] = adata.uns['umap']
    X_magic.obsm['X_umap'] = adata.obsm['X_umap']
    X_magic.write_h5ad(output_filename)

#  impute the data
impute_adata(adata=adata_li, output_filename=(dir_name + li_imputed_filename), impute_genes=genes1)
impute_adata(adata=adata_yost, output_filename=(dir_name + yost_imputed_filename), impute_genes=genes1)
impute_adata(adata=adata_wang, output_filename=(dir_name + wang_imputed_filename), impute_genes=genes1)