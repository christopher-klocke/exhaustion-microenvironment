# import statements
import scanpy as sc
import scanpy.external as sce


# setup
random1 = 1

# file paths
li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
li_processed_filename = 'li_processed.h5ad'
li_viz_ready_filename = 'li_viz_ready.h5ad'

# load the data
adata_li = sc.read_h5ad(li_path + li_processed_filename)

# process the data

# run PCA
print('running PCA...')
sc.tl.pca(adata_li, random_state=random1)

# run Harmony batch correction by sample
print('running Harmony...')
sce.pp.harmony_integrate(adata_li, 'patient')

adata_li.obsm['X_pca_original'] = adata_li.obsm['X_pca']
adata_li.obsm['X_pca'] = adata_li.obsm['X_pca_harmony']

# compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(adata_li, random_state=random1)

# perform Leiden clustering
print('running Leiden algorithm...')
sc.tl.leiden(adata_li, random_state=random1)

# run PAGA trajectory inference
print('inferring trajectories with PAGA...')
sc.tl.paga(adata_li)
sc.pl.paga(adata_li, plot=False)

# run UMAP dimensionality reduction
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(adata_li, init_pos='paga', random_state=random1)

# generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(adata_li, 'leiden', method='t-test')

# save to output .h5ad file
print('writing to output file...')
adata_li.write_h5ad(li_path + li_viz_ready_filename)
