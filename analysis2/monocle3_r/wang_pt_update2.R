# save-as of 'wang_pt_update.R' -- see file for creation of CDS object

library(monocle3)
library(tidyverse)

# LOAD FILE / SETUP
sample_id <- "batch"
genes_subset <- c("TCF7", "TOX", "LAG3", "PDCD1", "TIGIT")

## LOAD FILE
cds <- readRDS(file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/monocle3_runA_cds_wang.RDS")

## OPTIMIZE PSEUDOTIME
cds <- preprocess_cds(cds, method="PCA", num_dim=10) ## PCA is default
cds <- align_cds(cds, alignment_group=sample_id)
cds <- reduce_dimension(cds, umap.metric='euclidean')
cds <- cluster_cells(cds, random_seed=1, cluster_method="leiden", k=20, reduction_method="UMAP", resolution=0.0001)
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds, reduction_method="UMAP")
cds_subset <- cds[genes_subset, ]
plot_genes_in_pseudotime(cds_subset)
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runB/genes_v_pseudotime1.pdf")
pt_save <- pseudotime(cds, reduction_method = "UMAP") ## extract pseudotime values
write.csv(pt_save, file = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runB/pseudotime.csv")
saveRDS(cds, file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/monocle3_runB_cds_wang_processed.RDS")

##
cds <- readRDS(file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/monocle3_runB_cds_wang_processed.RDS")
umap_embeddings <- SingleCellExperiment::reducedDims(cds)[['UMAP']]  # extract UMAP embeddings
partitions <- partitions(cds, reduction_method="UMAP")  # extract partitions
write.csv(umap_embeddings, "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runB/umap_embeddings.csv")
write.csv(partitions, "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runB/partitions.csv")
