"""-pre-processing from yost datasets
-continuation from data download steps
-input file is .h5ad created via 'data_download_yost.py' (called by 'data_download.sh')
"""
# import statements
import scanpy as sc


# file paths
yost_data_path = '/Users/klockec/Documents/data/yost_data/download/' ## file path for rws07890
out_path = '/Users/klockec/Documents/data/yost_data/' ## file path for rws07890
bcc_scRNAseq_counts_filename = 'GSE123813_bcc_scRNA_counts.txt'
scc_scRNAseq_counts_filename = 'GSE123813_scc_scRNA_counts.txt'
yost_bcc_counts_filename = 'yost_bcc_counts.h5ad'
yost_scc_counts_filename = 'yost_scc_counts.h5ad'

# load the data
ad1 = sc.read((yost_data_path + bcc_scRNAseq_counts_filename))
ad2 = sc.read((yost_data_path + scc_scRNAseq_counts_filename))

# write the data to .h5ad files
ad1.T.write_h5ad(out_path + yost_bcc_counts_filename)
ad2.T.write_h5ad(out_path + yost_scc_counts_filename)
