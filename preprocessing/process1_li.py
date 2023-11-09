"""-pre-processing from li dataset
-continuation from data download steps
-input file is .h5ad created via 'data_download_li_add_metadata.py' (called by 'data_download.sh')
"""
# import statements
import scanpy as sc


# file paths
li_out_dir_name = '/Users/klockec/Documents/data/li_data/'
li_out_metadata_filename = 'li_raw_concat2_metadata.h5ad'
adata_li_subsetted_filename = 'adata_li_subsetted.h5ad'

# load the data
adata_li = sc.read_h5ad(li_out_dir_name + li_out_metadata_filename)

# process the data
untreated = ['p1', 'p3', 'p11', 'p12', 'p12', 'p13', 'p13_PBMC', 'p15', 'p16', 'p17', 'p17_PBMC', 'p18', 'p19', 'p21', 'p23', 'p24', 'p25', 'p26', 'p27', 'p27_PBMC']
adata_li_subsetted = adata_li[[x in untreated for x in adata_li.obs['patient']]]

# write the data to .h5ad file
adata_li_subsetted.write_h5ad(li_out_dir_name + adata_li_subsetted_filename)
