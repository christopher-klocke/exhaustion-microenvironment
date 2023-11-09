"""
perform imputation with the MAGIC algorithm

save-as of 'impute3.py'

input files from 'remove_rerun_yost_1b.py'
"""

## import statements
import scanpy as sc
import magic

## setup
random1 = 1

## file paths
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # file path for local
yost_bcc_full_filename = 'data/yost_data/yost_bcc_full_processed_dd.h5ad'
yost_scc_full_filename = 'data/yost_data/yost_scc_full_processed_dd.h5ad'
yost_bcc_imputed_filename = 'data/yost_data/yost_bcc_dd_imputed.h5ad'
yost_scc_imputed_filename = 'data/yost_data/yost_scc_dd_imputed.h5ad'

## load the data
adata_yost_bcc = sc.read_h5ad(dir + yost_bcc_full_filename)
adata_yost_scc = sc.read_h5ad(dir + yost_scc_full_filename)

## list of genes for imputation
gene_list = ['CD8A', 'CD8B', 'GZMB', 'GZMK', 'TOX', 'TCF7', 'CD4', 'FOXP3', 'LAG3', 'TIGIT', 'PDCD1', 'CTLA4', 'CD3E', 'GNLY', 'CCL5', 'NKG7', 'B2M', 'TYROBP', 'IL32', 'SSR4', 'CD74', 'HLA-DRA', 'IL7R', 'RPS12', 'LYZ', 'FTL', 'GZMH', 'CST7', 'MS4A1', 'CD79A', 'CD79B', 'HLA-DRB1', 'CD37', 'BANK1', 'CST3', 'LTB', 'MZB1', 'IGJ', 'FTH1', 'CD14', 'CSF1R', 'IL2RA', 'BATF', 'CX3CR1', 'CXCR3', 'CXCR4', 'CXCR5', 'IRF4', 'HAVCR2', 'NFATC1', 'CCR7', 'SELL', 'TNFRSF9', 'TBX21', 'EOMES', 'PRDM1', 'SLAMF6', 'CD101', 'CCR5', 'IFNG', 'TOX2', 'CD19', 'PTPRC', 'KLRK1', 'TNFRSF17', 'CD80', 'CD86', 'KIT', 'CD34', 'FLT3', 'SIRPA', 'NCAM1', 'FUT4', 'GATA1', 'GATA2', 'GATA3', 'TNF', 'LY6A', 'KLRB1', 'MERTK', 'CCR3', 'TLR7', 'TLR9', 'XCR1', 'CLEC9A', 'CD27', 'CD38', 'ADGRE1', 'FCGR3', 'ZAP70', 'LCK', 'PLCG1', 'TET2', 'LYZ2', 'PEAR1', 'SPI1', 'CEBPA', 'VSIR', 'ITGAM', 'ITGAX', 'ZBTB46', 'CCR9', 'VCAN', 'SYNE2', 'LCP1', 'PRRC2C', 'GZMA', 'VIM', 'CD63', 'TCF4', 'MIF', 'TRAC', 'CD2', 'IGHG3', 'IGHG1', 'MTRNR2L2', 'CD3G', 'TRBC2', 'TMSB4X', 'IFI30', 'CD22', 'CD68', 'HLA-A', 'TRBV20-1', 'SPOCK2', 'SAT1', 'IGKC', 'LMNA', 'MYADM', 'RGCC', 'LGALS1', 'IGHA1', 'STMN1', 'HMGN2', 'HMGB1', 'FKBP11', 'HLA-DRB5', 'SEC61B', 'PABPC1', 'S100A4', 'CXCL13', 'NR3C1', 'DNAJB1', 'HSPA1B', 'HSPA1A', 'CTSS', 'HLA-DPB1', 'HLA-DPA1', 'TNFAIP3', 'HSP90B1', 'LST1', 'AIF1', 'COTL1', 'FCER1G', 'HBA2', 'HBB', 'ALAS2', 'SDPR', 'HIST1H2AC', 'NAP1L1']

## impute with MAGIC algorithm
X_yost_bcc = adata_yost_bcc
magic_operator_bcc = magic.MAGIC()
X_yost_bcc_magic = magic_operator_bcc.fit_transform(X_yost_bcc, genes=gene_list)
X_yost_bcc_magic.uns['umap'] = adata_yost_bcc.uns['umap']
X_yost_bcc_magic.obsm['X_umap'] = adata_yost_bcc.obsm['X_umap']

X_yost_scc = adata_yost_scc
magic_operator_scc = magic.MAGIC()
X_yost_scc_magic = magic_operator_scc.fit_transform(X_yost_scc, genes=gene_list)
X_yost_scc_magic.uns['umap'] = adata_yost_scc.uns['umap']
X_yost_scc_magic.obsm['X_umap'] = adata_yost_scc.obsm['X_umap']

## write to output files
X_yost_bcc_magic.write_h5ad(dir + yost_bcc_imputed_filename)
X_yost_scc_magic.write_h5ad(dir + yost_scc_imputed_filename)
