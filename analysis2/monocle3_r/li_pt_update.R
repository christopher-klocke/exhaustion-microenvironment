library(reticulate)
library(monocle3)
library(tidyverse)
use_condaenv("sc")
sc <- import("scanpy")

# LOAD FILE / SETUP
#cds <- readRDS(file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/monocle3_run3_cds_li1.RDS")
sample_id <- "patient"
genes_subset <- c("TCF7", "TOX", "LAG3", "PDCD1", "TIGIT")

## JUST LI DATASET, NEW ROOT CELL
#library(reticulate)
#library(monocle3)
#library(tidyverse)
#use_condaenv("sc")
#sc <- import("scanpy")
##adata_li <- sc$read_h5ad("/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_t_modified.h5ad")
#adata_li <- sc$read_h5ad("/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_pt_input.h5ad")
#expression_matrix_li <- t(adata_li$X)
#genes_li <- as.data.frame(rownames(adata_li$var))
#colnames(genes_li) <- c("gene_short_name")
#rownames(genes_li) <- genes_li$gene_short_name
#cells_li <- as.data.frame(rownames(adata_li$obs))
#colnames(cells_li) <- c("CellID")
#rownames(cells_li) <- cells_li$CellID
#cells_li["patient"] <- adata_li$obs$patient #########
#cds_li <- new_cell_data_set(expression_matrix_li, gene_metadata = genes_li, cell_metadata = cells_li)
#saveRDS(cds_li, file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/monocle3_runA_cds_li1.RDS")

## LOAD FILE
cds <- readRDS(file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/monocle3_runA_cds_li1.RDS")

## OPTIMIZE PSEUDOTIME
genes_subset <- c("TCF7", "TOX", "LAG3", "PDCD1", "TIGIT")
#cds <- preprocess_cds(cds, method="PCA", num_dim=50) ## PCA is default
#cds <- preprocess_cds(cds, method="PCA", num_dim=10) ## PCA is default
cds <- preprocess_cds(cds, method="PCA", num_dim=18) ## PCA is default
plot_pc_variance_explained(cds)  # update to 8 PCs used
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/pcs_variance.pdf")
cds <- align_cds(cds, alignment_group=sample_id)
plot_cells(cds, color_cells_by="patient", reduction_method="PCA")
cds <- reduce_dimension(cds, umap.metric='euclidean')
plot_cells(cds, genes=genes_subset)
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/genes1.pdf")
plot_cells(cds, color_cells_by=sample_id, reduction_method="UMAP")
cds <- cluster_cells(cds, random_seed=1, cluster_method="leiden", k=20, reduction_method="UMAP", resolution=0.0001)
plot_genes_by_group(cds, markers=genes_subset)
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/genes_by_cluster1.pdf")
plot_cells(cds, color_cells_by="partition")
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/partition1.pdf")
#part_df <- as.data.frame(table(partitions(cds)))
#top_part <- part_df %>% filter(Freq == max(Freq)) %>% select(Var1)
#cds <- cds[, partitions(cds) == top_part$Var1]
#cds_subset <- cds[genes_subset, ]
#root_id <- "W687932-31-18"
#print(root_id)  # "WMC623395-171-0"
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by="cluster")
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/clusters1.pdf")
#cds <- order_cells(cds, reduction_method="UMAP", root_cells=root_id)
cds <- order_cells(cds, reduction_method="UMAP")
cds_subset <- cds[genes_subset, ]
plot_genes_in_pseudotime(cds_subset)  # https://www.rdocumentation.org/packages/monocle3/versions/1.0.0/topics/plot_genes_in_pseudotime
ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/genes_v_pseudotime1.pdf")
#plot_cells(cds, color_cells_by="cluster")
#ggsave("/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/clusters2.pdf")
pt_save <- pseudotime(cds, reduction_method = "UMAP") ## extract pseudotime values
write.csv(pt_save, file = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/pseudotime.csv")
saveRDS(cds, file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/monocle3_runA_cds_li1_processed.RDS")

##
cds <- readRDS(file="/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/monocle3_runA_cds_li1_processed.RDS")
umap_embeddings <- SingleCellExperiment::reducedDims(cds)[['UMAP']]  # extract UMAP embeddings
partitions <- partitions(cds, reduction_method="UMAP")  # extract partitions
write.csv(umap_embeddings, "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/umap_embeddings.csv")
write.csv(partitions, "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/partitions.csv")

