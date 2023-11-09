# import statements
import scanpy as sc

# file paths
yost_path = '/Users/klockec/Documents/data/yost_data/' ## file path for rws07890
yost_bcc_counts_pc_filename = 'yost_bcc_counts_pc.h5ad'
yost_scc_counts_pc_filename = 'yost_scc_counts_pc.h5ad'
yost_bcc_pyscenic_filename = 'yost_bcc_pyscenic_input1.h5ad'
yost_scc_pyscenic_filename = 'yost_scc_pyscenic_input1.h5ad'
yost_bcc_processed_filename = 'yost_bcc_processed.h5ad'
yost_scc_processed_filename = 'yost_scc_processed.h5ad'

# read in the data
adata_yost_bcc_pc = sc.read_h5ad(yost_path + yost_bcc_counts_pc_filename)
adata_yost_scc_pc = sc.read_h5ad(yost_path + yost_scc_counts_pc_filename)

# perform pre-processing on yost bcc data
sc.pp.filter_cells(adata_yost_bcc_pc, min_genes=200)
sc.pp.filter_genes(adata_yost_bcc_pc, min_cells=3)
adata_yost_bcc_pc.var['mt'] = adata_yost_bcc_pc.var_names.str.startswith('MT-') | adata_yost_bcc_pc.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_yost_bcc_pc, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_yost_bcc_pc = adata_yost_bcc_pc[adata_yost_bcc_pc.obs.n_genes_by_counts < 6500, :]
adata_yost_bcc_pc = adata_yost_bcc_pc[adata_yost_bcc_pc.obs.pct_counts_mt < 10, :]

# perform pre-processing on yost scc data
sc.pp.filter_cells(adata_yost_scc_pc, min_genes=200)
sc.pp.filter_genes(adata_yost_scc_pc, min_cells=3)
adata_yost_scc_pc.var['mt'] = adata_yost_scc_pc.var_names.str.startswith('MT-') | adata_yost_scc_pc.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_yost_scc_pc, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_yost_scc_pc = adata_yost_scc_pc[adata_yost_scc_pc.obs.n_genes_by_counts < 5000, :]
adata_yost_scc_pc = adata_yost_scc_pc[adata_yost_scc_pc.obs.pct_counts_mt < 10, :]

# save yost bcc/scc for input into pyscenic
adata_yost_bcc_pc.write_h5ad(yost_path + yost_bcc_pyscenic_filename)
adata_yost_scc_pc.write_h5ad(yost_path + yost_scc_pyscenic_filename)

# continue yost bcc pre-processing
sc.pp.normalize_total(adata_yost_bcc_pc, target_sum=1e4)
sc.pp.log1p(adata_yost_bcc_pc)
sc.pp.highly_variable_genes(adata_yost_bcc_pc)
adata_yost_bcc_pc.raw = adata_yost_bcc_pc

# continue yost scc pre-processing
sc.pp.normalize_total(adata_yost_scc_pc, target_sum=1e4)
sc.pp.log1p(adata_yost_scc_pc)
sc.pp.highly_variable_genes(adata_yost_scc_pc)
adata_yost_scc_pc.raw = adata_yost_scc_pc

# write to output files
adata_yost_bcc_pc.write_h5ad(yost_path + yost_bcc_processed_filename)
adata_yost_scc_pc.write_h5ad(yost_path + yost_scc_processed_filename)
