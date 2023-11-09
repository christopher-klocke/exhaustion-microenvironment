"""from 'data_download_viral.ipynb' """
# import statements
from download_unzip import unzip
import scanpy as sc


# file paths
wang_dir = '/Users/klockec/Documents/data/wang_data/download/'
wang_out_dir = '/Users/klockec/Documents/data/wang_data/'
zipped_files_list2 = 'wang_zipped_filenames.txt'
wang_save_filename = 'wang_adata_concat1.h5ad'

# process the data
unzipped_files2 = []
with open(wang_dir + zipped_files_list2, 'r') as read_file: 
    for line in read_file: 
        line = line.rstrip()
        if line[:3] == 'GSM': 
            path1 = wang_dir + line
            unzip(path1)
            unzipped_files2.append(path1[:-3])

wang_prefixes = ['GSM4775588_C1', 'GSM4775589_Q2', 'GSM4775590_Q3', 'GSM4775591_Q1', 'GSM4775592_Q4', 'GSM4775593_Q5', 'GSM4775594_Q7']
anndata_objects_wang = []
for elem in wang_prefixes: 
    adata_x = sc.read_10x_mtx(path=wang_dir, prefix=elem)
    anndata_objects_wang.append(adata_x)

first_one_wang = anndata_objects_wang.pop(0)
adata_concat_wang = first_one_wang.concatenate(anndata_objects_wang, join='outer', fill_value=0)
adata_concat_wang.write_h5ad(wang_out_dir + wang_save_filename)
