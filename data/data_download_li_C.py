""" continuation of Li dataset download steps """
# import statements
import scanpy as sc


# paths and filenames
dir_name = '/Users/klockec/Documents/data/li_data/download/' ## file path for rws07890
li_out_dir_name = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
all_filenames_path = '/Users/klockec/Documents/data/li_data/download/sample_files_unzipped.txt'
li_out_filename2 = 'li_raw_concatenated2.h5ad'

# process the data
unzipped_filenames_list = []
with open(all_filenames_path, 'r') as input1: 
    for line in input1: 
        line = line.rstrip()
        if line[0:3] == 'GSM' and line [-3:] == 'txt': 
            unzipped_filenames_list.append(line)

anndata_objects = []
for elem in unzipped_filenames_list: 
    adata_x = sc.read_csv((dir_name + elem), delimiter='\t').T
    anndata_objects.append(adata_x)

first_one = anndata_objects.pop(0)
adata_concat2 = first_one.concatenate(anndata_objects, join='outer', fill_value=0)
adata_concat2.write_h5ad(li_out_dir_name + li_out_filename2)
