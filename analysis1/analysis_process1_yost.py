# import statements
import scanpy as sc
import scanpy.external as sce
import pandas as pd


# setup
random1 = 1

# file paths
yost_path = '/Users/klockec/Documents/data/yost_data/' ## file path for rws07890
yost_scc_out_filename = 'yost_scc_processed2.h5ad'
yost_bcc_out_filename = 'yost_bcc_processed2.h5ad'
yost_viz_ready_filename = 'yost_viz_ready.h5ad'

# load the data
yost_bcc = sc.read_h5ad(yost_path + yost_bcc_out_filename)
yost_scc = sc.read_h5ad(yost_path + yost_scc_out_filename)

# select just immune cells
yost_bcc_immune_types = ['B_cells_1', 'B_cells_2', 'CD4_T_cells', 'CD8_act_T_cells', 'CD8_ex_T_cells', 'CD8_mem_T_cells', 'DCs', 'Macrophages', 'NK_cells', 'Plasma_cells', 'Tcell_prolif', 'Tregs', 'pDCs']
yost_bcc_immune = yost_bcc[[x in yost_bcc_immune_types for x in list(yost_bcc.obs['cluster_original'])]]
yost_scc_immune = yost_scc

# merge the two datasets
yost_immune_concat = yost_bcc_immune.concatenate(yost_scc_immune, join='outer', fill_value=0) 
del yost_immune_concat.obsm
yost_immune_concat.var['gene_name'] = yost_immune_concat.var.index
hold1 = pd.DataFrame(yost_immune_concat.var['gene_name'])
del yost_immune_concat.var
yost_immune_concat.var = hold1
sc.pp.highly_variable_genes(yost_immune_concat)

# run PCA
print('running PCA...')
sc.tl.pca(yost_immune_concat, random_state=random1)

# perform batch correction by sample with Harmony algorithm
print('running Harmony...')
sce.pp.harmony_integrate(yost_immune_concat, 'sample', max_iter_harmony=50)
yost_immune = yost_immune_concat
yost_immune.obsm['X_pca_original'] = yost_immune.obsm['X_pca']
yost_immune.obsm['X_pca'] = yost_immune.obsm['X_pca_harmony']

# compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(yost_immune, random_state=random1)

# perform Leiden clustering
print('running Leiden algorithm...')
sc.tl.leiden(yost_immune, random_state=random1)

# run PAGA trajectory inference
print('inferring trajectories with PAGA...')
sc.tl.paga(yost_immune)
sc.pl.paga(yost_immune, plot=False)

# run UMAP dimensionality reduction
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(yost_immune, init_pos='paga', random_state=random1)

# generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(yost_immune, 'leiden', method='t-test')

# save to output .h5ad file
print('writing to output file...')
yost_immune.write_h5ad(yost_path + yost_viz_ready_filename)
