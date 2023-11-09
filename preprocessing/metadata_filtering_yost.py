# import statements
import scanpy as sc
import numpy as np


# file paths
yost_path = '/Users/klockec/Documents/data/yost_data/' ## file path for rws07890
yost_metadata_path = '/Users/klockec/Documents/data/yost_data/download/'
yost_bcc_processed_filename = 'yost_bcc_processed.h5ad'
yost_scc_processed_filename = 'yost_scc_processed.h5ad'
yost_bcc_metadata_filename = 'GSE123813_bcc_all_metadata.txt'
yost_scc_metadata_filename = 'GSE123813_scc_metadata.txt'
yost_scc_out_filename = 'yost_scc_processed2.h5ad'
yost_bcc_out_filename = 'yost_bcc_processed2.h5ad'

# load the data
adata_yost_bcc = sc.read_h5ad(yost_path + yost_bcc_processed_filename)
adata_yost_scc = sc.read_h5ad(yost_path + yost_scc_processed_filename)

# process the data
yost_bcc_sample_ids = adata_yost_bcc.obs.index
yost_scc_sample_ids = adata_yost_scc.obs.index
yost_bcc_just_pre = adata_yost_bcc[[x.split('.')[2][0:3] == 'pre' for x in yost_bcc_sample_ids], :]
yost_scc_just_pre = adata_yost_scc[[x.split('.')[2][0:3] == 'pre' for x in yost_scc_sample_ids], :]
yost_bcc_samples = [x.split('.')[1] for x in yost_bcc_just_pre.obs.index]
yost_bcc_just_pre.obs['sample'] = yost_bcc_samples
yost_bcc_tumor_type = [x.split('.')[0] for x in yost_bcc_just_pre.obs.index]
yost_bcc_just_pre.obs['tumor_type'] = yost_bcc_tumor_type
yost_scc_samples = [x.split('.')[1] for x in yost_scc_just_pre.obs.index]
yost_scc_just_pre.obs['sample'] = yost_scc_samples
yost_scc_tumor_type = [x.split('.')[0] for x in yost_scc_just_pre.obs.index]
yost_scc_just_pre.obs['tumor_type'] = yost_scc_tumor_type

# add yost bcc metadata
yost_bc_meta_dict = {}
with open((yost_metadata_path + yost_bcc_metadata_filename), 'r') as input1: 
    for line in input1: 
        line = line.rstrip()
        line = line.split('\t')
        # line[0] is cell.id, 1 is patient, 2 is treatment, 3 is sort, 4 is cluster, 5 is umap1, 6 is umap2
        if line[0] != 'cell.id': ## ignore first line
            yost_bc_meta_dict[line[0]] = (line[3], line[4], line[5], line[6])

sort_list = []
cluster_list = []
UMAP1_list = []
UMAP2_list = []

for x in yost_bcc_just_pre.obs.index: 
    sort_list.append(yost_bc_meta_dict[x][0])
    cluster_list.append(yost_bc_meta_dict[x][1])
    UMAP1_list.append(yost_bc_meta_dict[x][2])
    UMAP2_list.append(yost_bc_meta_dict[x][3])

yost_bcc_just_pre.obs['sort'] = sort_list
yost_bcc_just_pre.obs['cluster_original'] = cluster_list
yost_bcc_just_pre.obsm['X_umap'] = np.array(list(zip(UMAP1_list, UMAP2_list))).astype(np.float32)
# add yost scc metadata

yost_sc_meta_dict = {}
with open((yost_metadata_path + yost_scc_metadata_filename), 'r') as input1: 
    for line in input1: 
        line = line.rstrip()
        line = line.split('\t')
        # line[0] is cell.id, 1 is patient, 2 is treatment, 3 is cluster, 4 is umap1, 5 is umap2
        if line[0] != 'cell.id': ## ignore first line
            yost_sc_meta_dict[line[0]] = (line[3], line[4], line[5])

cluster_list1 = []
UMAP1_list1 = []
UMAP2_list1 = []
for x in yost_scc_just_pre.obs.index: 
    cluster_list1.append(yost_sc_meta_dict[x][0])
    UMAP1_list1.append(yost_sc_meta_dict[x][1])
    UMAP2_list1.append(yost_sc_meta_dict[x][2])

yost_scc_just_pre.obs['cluster_original'] = cluster_list1
yost_scc_just_pre.obsm['X_umap'] = np.array(list(zip(UMAP1_list1, UMAP2_list1))).astype(np.float32)

# write anndata objects to .h5ad files
yost_bcc_just_pre.write_h5ad(yost_path + yost_bcc_out_filename)
yost_scc_just_pre.write_h5ad(yost_path + yost_scc_out_filename)
