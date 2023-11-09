"""add metadata to anndata object for Li et al dataset """
# import statements
import scanpy as sc
import pandas as pd
from download_unzip import download_and_unzip


# file paths
li_out_dir_name = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
dir_name2 = li_out_dir_name ## REPLACE AND DELETE LATER
li_out_filename2 = 'li_raw_concatenated2.h5ad'
li_out_metadata_filename = 'li_raw_concat2_metadata.h5ad'
li_metadata1_download_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123139/matrix/GSE123139_series_matrix.txt.gz'
li_metadata1_zipped = dir_name2 + 'GSE123139_series_matrix.txt.gz'
li_metadata1_unzipped = dir_name2 + 'GSE123139_series_matrix.txt'
all_filenames_path = '/Users/klockec/Documents/data/li_data/download/sample_files_unzipped.txt'
# load in anndata object
adata_concat2 = sc.read_h5ad(li_out_dir_name + li_out_filename2)
# download the metadata
download_and_unzip(url=li_metadata1_download_url, unzipped=li_metadata1_unzipped, zipped=li_metadata1_zipped)
# process the metadata
rows_list = []
with open(li_metadata1_unzipped, 'r') as input1: 
    for line in input1: 
        line = line.rstrip()
        if line.startswith('!Sample_title'): 
            x = line.split('\t')[1:]
            y = []
            for elem in x:
                y.append(elem[1:-1])
            rows_list.append(y)            
        if line.startswith('!Sample_geo_accession'): 
            x = line.split('\t')[1:]
            y = []
            for elem in x:
                y.append(elem[1:-1])
            rows_list.append(y)       
        if line.startswith('!Sample_characteristics_ch1'):
            x = line.split('\t')[1:]
            y = []
            for elem in x:
                y.append(elem[1:-1])
            rows_list.append(y)   

df = pd.DataFrame(rows_list)
df.index = ['ab_id', 'geo_id', 'plate_id', 'amplification_batch', 'sample_source', 'patient_id', 'facs_gate']
metadata_lookup_dict = {}
geo_ids = []
for elem in list(df.loc['geo_id']): 
    geo_ids.append(elem)
plate_ids = []
for elem in list(df.loc['plate_id']): 
    plate_ids.append(elem[10:])
amplification_batches = []
for elem in list(df.loc['amplification_batch']): 
    amplification_batches.append(elem[21:])
sample_sources = []
for elem in list(df.loc['sample_source']): 
    sample_sources.append(elem[15:])
patient_ids = []
patients = []
for elem in list(df.loc['patient_id']): 
    patient_ids.append(elem[12:])
    patients.append(elem[12:].split('-')[0])
facs_ids = []
for elem in list(df.loc['facs_gate']):
    facs_ids.append(elem[11:])
counter = 0
while counter < 204: 
    metadata_lookup_dict[counter] = (geo_ids[counter], plate_ids[counter], amplification_batches[counter], sample_sources[counter], patient_ids[counter], patients[counter], facs_ids[counter])
    counter += 1
geo_ids = []
ab_ids = []
unzipped_filenames_list = []
with open(all_filenames_path, 'r') as input1: 
    for line in input1: 
        line = line.rstrip()
        if line[0:3] == 'GSM' and line [-3:] == 'txt': 
            unzipped_filenames_list.append(line)
for elem in unzipped_filenames_list: 
    ids = elem[:-4].split('_')
    geo_ids.append(ids[0])
    ab_ids.append(ids[1])
ids_dict = {}
counter = 0
while counter < 204: 
    ids_dict[counter] = (geo_ids[counter], ab_ids[counter])
    counter += 1
obs_geo_ids = []
obs_ab_ids = []
for elem in list(adata_concat2.obs['batch']): 
    obs_geo_ids.append(ids_dict[int(elem)][0])
    obs_ab_ids.append(ids_dict[int(elem)][1])
adata_concat2.obs['geo_id'] = obs_geo_ids
adata_concat2.obs['ab_id'] = obs_ab_ids
del adata_concat2.obs['ab_id']
batches = list(adata_concat2.obs['batch'])
cell_geo_ids = []
cell_plate_ids = []
cell_amplification_batches = []
cell_sample_source = []
cell_patient_id = []
cell_patient = []
cell_facs_gate = []
for elem in batches: 
    x = metadata_lookup_dict[int(elem)]
    cell_geo_ids.append(x[0])
    cell_plate_ids.append(x[1])
    cell_amplification_batches.append(x[2])
    cell_sample_source.append(x[3])
    cell_patient_id.append(x[4])
    cell_patient.append(x[5])
    cell_facs_gate.append(x[6])
adata_concat2.obs['plate_id'] = cell_plate_ids
adata_concat2.obs['amplification_batch'] = cell_amplification_batches
adata_concat2.obs['sample_source'] = cell_sample_source
adata_concat2.obs['patient_id'] = cell_patient_id
adata_concat2.obs['patient'] = cell_patient
adata_concat2.obs['facs_gate'] = cell_facs_gate

# write updated anndata object to file
adata_concat2.write_h5ad(li_out_dir_name + li_out_metadata_filename)
