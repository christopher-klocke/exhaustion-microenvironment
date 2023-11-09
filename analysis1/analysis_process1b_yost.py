# import statements
import scanpy as sc
import scanpy.external as sce
import pandas as pd


# setup
random1 = 1

# file paths
#yost_path = '/Users/klockec/Documents/' ## file path for rws07890
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # file path for local
yost_bcc_out_filename = 'data/yost_data/yost_bcc_processed2.h5ad'
yost_scc_out_filename = 'data/yost_data/yost_scc_processed2.h5ad'
yost_bcc_full_filename = 'data/yost_data/yost_bcc_full_processed.h5ad'
yost_scc_full_filename = 'data/yost_data/yost_scc_full_processed.h5ad'

# load the data
yost_bcc = sc.read_h5ad(dir + yost_bcc_out_filename)
yost_scc = sc.read_h5ad(dir + yost_scc_out_filename)

# select just immune cells
#yost_bcc_immune_types = ['B_cells_1', 'B_cells_2', 'CD4_T_cells', 'CD8_act_T_cells', 'CD8_ex_T_cells', 'CD8_mem_T_cells', 'DCs', 'Macrophages', 'NK_cells', 'Plasma_cells', 'Tcell_prolif', 'Tregs', 'pDCs']
#yost_bcc_immune = yost_bcc[[x in yost_bcc_immune_types for x in list(yost_bcc.obs['cluster_original'])]]
#yost_scc_immune = yost_scc

# merge the two datasets
#yost_immune_concat = yost_bcc_immune.concatenate(yost_scc_immune, join='outer', fill_value=0)
#del yost_immune_concat.obsm
#yost_immune_concat.var['gene_name'] = yost_immune_concat.var.index
#hold1 = pd.DataFrame(yost_immune_concat.var['gene_name'])
#del yost_immune_concat.var
#yost_immune_concat.var = hold1
#sc.pp.highly_variable_genes(yost_immune_concat)

# run PCA
print('running PCA...')
sc.tl.pca(yost_bcc, random_state=random1)
sc.tl.pca(yost_scc, random_state=random1)

# perform batch correction by sample with Harmony algorithm
print('running Harmony...')
sce.pp.harmony_integrate(yost_bcc, 'sample', max_iter_harmony=50)
yost_bcc.obsm['X_pca_original'] = yost_bcc.obsm['X_pca']
yost_bcc.obsm['X_pca'] = yost_bcc.obsm['X_pca_harmony']

sce.pp.harmony_integrate(yost_scc, 'sample', max_iter_harmony=50)
yost_scc.obsm['X_pca_original'] = yost_scc.obsm['X_pca']
yost_scc.obsm['X_pca'] = yost_scc.obsm['X_pca_harmony']

# compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(yost_bcc, random_state=random1)
sc.pp.neighbors(yost_scc, random_state=random1)

# perform Leiden clustering
print('running Leiden algorithm...')
sc.tl.leiden(yost_bcc, random_state=random1)
sc.tl.leiden(yost_scc, random_state=random1)

# run PAGA trajectory inference
print('inferring trajectories with PAGA...')
sc.tl.paga(yost_bcc)
sc.pl.paga(yost_bcc, plot=False)

sc.tl.paga(yost_scc)
sc.pl.paga(yost_scc, plot=False)

# run UMAP dimensionality reduction
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(yost_bcc, init_pos='paga', random_state=random1)
sc.tl.umap(yost_scc, init_pos='paga', random_state=random1)

# generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(yost_immune, 'leiden', method='t-test')

# save to output .h5ad file
print('writing to output file...')
yost_bcc.write_h5ad(dir + yost_bcc_full_filename)
yost_scc.write_h5ad(dir + yost_scc_full_filename)
