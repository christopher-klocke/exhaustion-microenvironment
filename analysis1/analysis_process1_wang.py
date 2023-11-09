# import statements
import scanpy as sc
import scanpy.external as sce


# setup
random1 = 1

# file paths
wang_dir = '/Users/klockec/Documents/data/wang_data/'
wang_no6_filename = 'adata_wang_no6.h5ad'
wang_viz_ready_filename = 'wang_viz_ready.h5ad'

# load the data
adata_wang = sc.read_h5ad(wang_dir + wang_no6_filename)

# process the data

# run PCA
print('running PCA...')
sc.tl.pca(adata_wang, random_state=random1)

# run Harmony batch correction by sample
print('running Harmony...')
sce.pp.harmony_integrate(adata_wang, 'batch')

adata_wang.obsm['X_pca_original'] = adata_wang.obsm['X_pca']
adata_wang.obsm['X_pca'] = adata_wang.obsm['X_pca_harmony']

# compute neighbor graph
print('computing neighbor graph...')
sc.pp.neighbors(adata_wang, random_state=random1)

# perform Leiden clustering
print('running Leiden algorithm...')
sc.tl.leiden(adata_wang, random_state=random1)

# run PAGA trajectory inference
print('inferring trajectories with PAGA...')
sc.tl.paga(adata_wang)
sc.pl.paga(adata_wang, plot=False)

# run UMAP dimensionality reduction
print('performing dimensionality reduction with UMAP...')
sc.tl.umap(adata_wang, init_pos='paga', random_state=random1)

# generate marker genes per cluster
#print('calculating marker genes...')
#sc.tl.rank_genes_groups(adata_li, 'leiden', method='t-test')

# save to output .h5ad file
print('writing to output file...')
adata_wang.write_h5ad(wang_dir + wang_viz_ready_filename)
