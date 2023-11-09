library(reticulate)
library(monocle3)
library(tidyverse)
use_condaenv("sc")
sc <- import("scanpy")

# LOAD FILE / SETUP
sample_id <- "sample"
genes_subset <- c("TCF7", "TOX", "LAG3", "PDCD1", "TIGIT")
#genes_subset <- c("TCF7", "TOX", "LAG3", "PDCD1", "TIGIT", "CD8A", "CD8B", "CD4", "FOXP3")

# prep cds object
##adata <- sc$read_h5ad("/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/yost_t_modified_v2.h5ad")
#adata <- sc$read_h5ad("/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/yost_pt_input.h5ad")
#expression_matrix <- t(adata$X)
#genes_add <- as.data.frame(rownames(adata$var))
#colnames(genes_add) <- c("gene_short_name")
#rownames(genes_add) <- genes_add$gene_short_name
#cells_add <- as.data.frame(rownames(adata$obs))
#colnames(cells_add) <- c("CellID")
#rownames(cells_add) <- cells_add$CellID
#cells_add["sample"] <- adata$obs$sample
#cds <- new_cell_data_set(expression_matrix, gene_metadata = genes_add, cell_metadata = cells_add)
#saveRDS(cds, file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/monocle3_runA_cds_yost.RDS")

## LOAD FILE
cds <- readRDS(file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/monocle3_runA_cds_yost.RDS")

## OPTIMIZE PSEUDOTIME
#cds <- preprocess_cds(cds, method="PCA", num_dim=50) ## PCA is default
cds <- preprocess_cds(cds, method="PCA", num_dim=15) ## PCA is default
#cds <- preprocess_cds(cds, method="PCA", num_dim=7) ## PCA is default
plot_pc_variance_explained(cds)  # update to n PCs used
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/pcs_variance.pdf")
cds <- align_cds(cds, alignment_group=sample_id)
plot_cells(cds, color_cells_by=sample_id, reduction_method="PCA")
cds <- reduce_dimension(cds, umap.metric='euclidean')
plot_cells(cds, color_cells_by=sample_id, reduction_method="UMAP")
plot_cells(cds, genes=genes_subset)
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/genes1.pdf")
cds <- cluster_cells(cds, random_seed=1, cluster_method="leiden", k=20, reduction_method="UMAP", resolution=0.0001)
#cds <- cluster_cells(cds, random_seed=1, cluster_method="leiden", k=30, reduction_method="UMAP", resolution=0.0001)
#plot_genes_by_group(cds, markers=genes_subset)
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/genes_by_cluster1.pdf")
#plot_cells(cds, color_cells_by="partition")
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/partition1.pdf")
#part_df <- as.data.frame(table(partitions(cds)))
#top_part <- part_df %>% filter(Freq == max(Freq)) %>% select(Var1)
#cds <- cds[, partitions(cds) == top_part$Var1]
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/genes2.pdf")
#root_id <- "W687932-31-18"
#print(root_id)  # "WMC623395-171-0"
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by="cluster")
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/clusters1.pdf")
#cds <- order_cells(cds, reduction_method="UMAP", root_cells=root_id)
cds <- order_cells(cds, reduction_method="UMAP")
#genes_subset <- c("TCF7", "TOX", "LAG3", "PDCD1", "TIGIT")
#cds_subset <- cds[genes_subset, ]

cds_part1 <- cds[, partitions(cds) == 1]
cds_part1_subset <- cds_part1[genes_subset, ]
plot_genes_in_pseudotime(cds_part1_subset)

plot_genes_in_pseudotime(cds_subset)  # https://www.rdocumentation.org/packages/monocle3/versions/1.0.0/topics/plot_genes_in_pseudotime
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/genes_v_pseudotime1.pdf")
#plot_cells(cds, color_cells_by="cluster", cell_size=0.6)
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/clusters2.pdf")
pt_save <- pseudotime(cds, reduction_method = "UMAP") ## extract pseudotime values
write.csv(pt_save, file = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/pseudotime.csv")
saveRDS(cds, file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/monocle3_runA_cds_yost_processed.RDS")

##

##
cds <- readRDS(file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/monocle3_runA_cds_yost_processed.RDS")
umap_embeddings <- SingleCellExperiment::reducedDims(cds)[['UMAP']]  # extract UMAP embeddings
partitions <- partitions(cds, reduction_method="UMAP")  # extract partitions
write.csv(umap_embeddings, "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/umap_embeddings.csv")
write.csv(partitions, "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/yost_runA/partitions.csv")

