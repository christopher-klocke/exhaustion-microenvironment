"""use UniProt gene names to restrict anndata objects to just protein-coding genes """
# import statements
import scanpy as sc

# file paths
genes_dir = '/Users/klockec/Documents/data/gene_names/'
yost_path = '/Users/klockec/Documents/data/yost_data/'
li_out_dir_name = '/Users/klockec/Documents/data/li_data/'
unzipped_filename = genes_dir + 'uniprot_human_reviewed.tsv'
yost_bcc_counts_filename = 'yost_bcc_counts.h5ad'
yost_scc_counts_filename = 'yost_scc_counts.h5ad'
adata_li_subsetted_filename = 'adata_li_subsetted.h5ad'
yost_bcc_counts_pc_filename = 'yost_bcc_counts_pc.h5ad'
yost_scc_counts_pc_filename = 'yost_scc_counts_pc.h5ad'
adata_li_subsetted_pc_filename = 'adata_li_subsetted_pc.h5ad'

# process the data
gene_names1 = []
with open(unzipped_filename, 'r') as input1: 
    for line in input1: 
        line = line.rstrip()
        if line[:5] != 'Entry': 
            line = line.split(sep='\t')
            alts = line[4].split(sep=' ')
            for x in alts: 
                gene_names1.append(x)

gene_names1_set = set(gene_names1)

# load in anndata objects
adata_yost_bcc = sc.read_h5ad(yost_path + yost_bcc_counts_filename)
adata_yost_scc = sc.read_h5ad(yost_path + yost_scc_counts_filename)
adata_li_subsetted = sc.read_h5ad(li_out_dir_name + adata_li_subsetted_filename)

# filter to just protein-coding genes
adata_yost_bcc_pc = adata_yost_bcc[:, [x in gene_names1_set for x in list(adata_yost_bcc.var.index)]]
adata_yost_scc_pc = adata_yost_scc[:, [x in gene_names1_set for x in list(adata_yost_scc.var.index)]]
adata_li_subsetted_pc = adata_li_subsetted[:, [x in gene_names1_set for x in list(adata_li_subsetted.var.index)]]

# write to output files
adata_yost_bcc_pc.write_h5ad(yost_path + yost_bcc_counts_pc_filename)
adata_yost_scc_pc.write_h5ad(yost_path + yost_scc_counts_pc_filename)
adata_li_subsetted_pc.write_h5ad(li_out_dir_name + adata_li_subsetted_pc_filename)
