library(reticulate)
library(monocle3)
library(tidyverse)
use_condaenv("sc")
sc <- import("scanpy")

# LOAD FILE / SETUP
sample_id <- "batch"
genes_subset <- c("TCF7", "TOX", "LAG3", "PDCD1", "TIGIT")
#genes_subset <- c("TCF7", "TOX", "LAG3", "PDCD1", "TIGIT", "CD8A", "CD8B", "CD4", "FOXP3")

# prep cds object
##adata <- sc$read_h5ad("/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/wang_t_modified.h5ad")
#adata <- sc$read_h5ad("/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/wang_pt_input.h5ad")
##expression_matrix <- t(adata$X)
#expression_matrix <- as.matrix(t(adata$X)) ## mod for different data type --
#genes_add <- as.data.frame(rownames(adata$var))
#colnames(genes_add) <- c("gene_short_name")
#rownames(genes_add) <- genes_add$gene_short_name
#cells_add <- as.data.frame(rownames(adata$obs))
#colnames(cells_add) <- c("CellID")
#rownames(cells_add) <- cells_add$CellID
#cells_add["batch"] <- adata$obs$batch
#cds <- new_cell_data_set(expression_matrix, gene_metadata = genes_add, cell_metadata = cells_add)
#saveRDS(cds, file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/monocle3_runA_cds_wang.RDS")

## LOAD FILE
cds <- readRDS(file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/monocle3_runA_cds_wang.RDS")

## OPTIMIZE PSEUDOTIME
#cds <- preprocess_cds(cds, method="PCA", num_dim=50) ## PCA is default; 7 or 10, probably 10
#cds <- preprocess_cds(cds, method="PCA", num_dim=8) ## PCA is default
cds <- preprocess_cds(cds, method="PCA", num_dim=10) ## PCA is default
plot_pc_variance_explained(cds)  # update to 8 PCs used
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/pcs_variance.pdf")
cds <- align_cds(cds, alignment_group=sample_id)
plot_cells(cds, color_cells_by=sample_id, reduction_method="PCA")
cds <- reduce_dimension(cds, umap.metric='euclidean')
plot_cells(cds, genes=genes_subset)
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/genes1a.pdf")
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/genes1b.pdf")
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/genes1.pdf")
#cds <- reduce_dimension(cds, umap.metric='euclidean', preprocess_method='PCA')
#cds <- reduce_dimension(cds, umap.metric='euclidean', max_components=1)
#cds <- reduce_dimension(cds, umap.metric='euclidean', umap.n_neighbors=30)
#cds <- reduce_dimension(cds, umap.metric='euclidean', umap.n_neighbors=10)
#cds <- reduce_dimension(cds, umap.metric='euclidean', umap.n_neighbors=8)
#cds <- cluster_cells(cds, random_seed=1, cluster_method="leiden", k=15, reduction_method="UMAP", resolution=0.001)
#plot_cells(cds, color_cells_by=sample_id, reduction_method="UMAP")
cds <- cluster_cells(cds, random_seed=1, cluster_method="leiden", k=20, reduction_method="UMAP", resolution=0.0001)
#cds <- cluster_cells(cds, random_seed=1, cluster_method="leiden", k=30, reduction_method="UMAP", resolution=0.0001)
#plot_genes_by_group(cds, markers=genes_subset)
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/genes_by_cluster1.pdf")
#plot_cells(cds, color_cells_by="partition")
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/partition1.pdf")
#part_df <- as.data.frame(table(partitions(cds)))
#top_part <- part_df %>% filter(Freq == max(Freq)) %>% select(Var1)
#cds <- cds[, partitions(cds) == top_part$Var1]
#root_id <- ""
#cds <- learn_graph(cds)
cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds, color_cells_by="cluster")
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/clusters1.pdf")
#cds <- order_cells(cds, reduction_method="UMAP", root_cells=root_id)
cds <- order_cells(cds, reduction_method="UMAP")
#cds_part1 <- cds[, partitions(cds) == 2]
#cds_part1_subset <- cds_part1[genes_subset, ]
cds_subset <- cds[genes_subset, ]
#plot_genes_in_pseudotime(cds_part1_subset)
plot_genes_in_pseudotime(cds_subset)
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/genes_v_pseudotime1.pdf")
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/genes_v_pseudotime1b.pdf")
#plot_cells(cds, color_cells_by="cluster")
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/clusters2.pdf")
pt_save <- pseudotime(cds, reduction_method = "UMAP") ## extract pseudotime values
write.csv(pt_save, file = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/wang_runA/pseudotime.csv")
saveRDS(cds, file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/monocle3_runA_cds_wang_processed.RDS")
