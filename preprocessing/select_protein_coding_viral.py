"""use UniProt gene names to restrict anndata objects to just protein-coding genes """
# import statements
import scanpy as sc


# file paths
wang_dir = '/Users/klockec/Documents/data/wang_data/'
genes_dir = '/Users/klockec/Documents/data/gene_names/'
wang_save_filename = 'wang_adata_concat1.h5ad'
unzipped_filename = genes_dir + 'uniprot_human_reviewed.tsv'
wang_pc_filename = 'wang_adata_pc.h5ad'

# load the data
adata_wang = sc.read_h5ad(wang_dir + wang_save_filename)

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

# filter to just protein-coding genes
adata_wang_pc = adata_wang[:, [x in gene_names1_set for x in list(adata_wang.var.index)]]

# write to output files
adata_wang_pc.write_h5ad(wang_dir + wang_pc_filename)
