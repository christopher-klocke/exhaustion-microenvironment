"""-run with 'sc2' conda env
-select just CD8 T cells from anndata objects without preprocessing, using index of CD8 T anndata objects to subset
-do not perform preprocessing (e.g. normalize, logarithmize) twice on data during trajectory inference
"""

# import statements
import scanpy as sc

# setup
dir_name = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'
#li_index_filename = 'data/li_data/li_t_modified.h5ad'
yost_index_filename = 'data/yost_data/yost_t_modified_v2.h5ad'
#wang_index_filename = 'data/wang_data/wang_t_modified.h5ad'
#adata_li_pyscenic_filename = 'data/li_data/adata_li_pyscenic_input1.h5ad'
yost_bcc_pyscenic_filename = 'data/yost_data/yost_bcc_pyscenic_input1.h5ad'
yost_scc_pyscenic_filename = 'data/yost_data/yost_scc_pyscenic_input1.h5ad'
#adata_wang_pyscenic_filename = 'data/wang_data/adata_wang_pyscenic_input1.h5ad'
#li_monocle_input = 'data/li_data/li_pt_input.h5ad'
yost_monocle_input = 'data/yost_data/yost_pt_input.h5ad'
#wang_monocle_input = 'data/wang_data/wang_pt_input.h5ad'

# load the data
#adata_li_index = sc.read_h5ad(dir_name + li_index_filename)
adata_yost_index = sc.read_h5ad(dir_name + yost_index_filename)
#adata_wang_index = sc.read_h5ad(dir_name + wang_index_filename)
#adata_li_unprocessed = sc.read_h5ad(dir_name + adata_li_pyscenic_filename)
adata_yost_unprocessed_bcc = sc.read_h5ad(dir_name + yost_bcc_pyscenic_filename)
adata_yost_unprocessed_scc = sc.read_h5ad(dir_name + yost_scc_pyscenic_filename)
#adata_wang_unprocessed = sc.read_h5ad(dir_name + adata_wang_pyscenic_filename)

adata_yost_unprocessed = adata_yost_unprocessed_bcc.concatenate(adata_yost_unprocessed_scc)  # concatenate yost bcc and scc objects

# subset the data
#li_ind = ['-'.join(x.split('-')[:-1]) for x in adata_li_index.obs.index]
#adata_li_cut = adata_li_unprocessed[adata_li_unprocessed.obs.index.isin(li_ind)]
yost_ind = ['-'.join(x.split('-')[:-1]) for x in adata_yost_index.obs.index]
adata_yost_index.obs.index = yost_ind
yost_re_ind = ['-'.join(x.split('-')[:-1]) for x in adata_yost_unprocessed.obs.index]
adata_yost_unprocessed.obs.index = yost_re_ind
adata_yost_cut = adata_yost_unprocessed[adata_yost_unprocessed.obs.index.isin(yost_ind)]

# add 'sample' column
#df = adata_yost_cut.obs
#df['CellID'] = df.index
adata_yost_cut.obs['CellID'] = adata_yost_cut.obs.index
df2 = adata_yost_index.obs
df2['CellID'] = df2.index
df2 = df2[['sample', 'CellID']]
#df3 = df.join(df2, on='CellID', how='left')  # https://stackoverflow.com/questions/26645515/pandas-join-issue-columns-overlap-but-no-suffix-specified
#df3 = df.merge(df2, on='CellID', how='left')
df3 = adata_yost_cut.obs.merge(df2, on='CellID', how='left')
adata_yost_cut.obs = df3
adata_yost_cut.obs.index = df3['CellID']

#wang_ind = ['-'.join(x.split('-')[:-1]) for x in adata_wang_index.obs.index]
#adata_wang_cut = adata_wang_unprocessed[adata_wang_unprocessed.obs.index.isin(wang_ind)]

# write monocle input files
#adata_li_cut.write_h5ad(dir_name + li_monocle_input)
adata_yost_cut.write_h5ad(dir_name + yost_monocle_input)
#adata_wang_cut.write_h5ad(dir_name + wang_monocle_input)


