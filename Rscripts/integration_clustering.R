# Script to integrate SAT snRNA-seq data across batches, run Harmony, cluster integrated data,
# and annotate major cell populations (adipocytes, ASPCs, lymphoid, myeloid, and vascular cells).

library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(future)
library(SingleR)
library(SeuratWrappers)
library(scRNAseq)
library(scater)
library(tidyr)
library(BiocParallel)
library(sctransform)
library(data.table)
library(harmony)
library(ggpubr)
library(clustree)


# Read list of Seurat objects from all batches
seur_list <- readRDS("sat_snrna_list.rds")

# Merge all Seurat objects
seur_merged <- merge(
  seur_list[[1]],
  y = c(seur_list[c(2:length(seur_list))]),
  project = "merged_sat"
)


### PRE-PROCESSING AND DATA INTEGRATION ############################################################
# Read list of ribosomal and mitochondrial genes
rp_gene <- fread("ribosomal_genes.txt", header = TRUE)
mt_gene <- fread("mitochondrial_genes.txt", header = TRUE)

# Keep genes with at least 3 counts in at least 3 cells
gene_cnt <- seur_merged@assays$RNA@counts
keep <- rowSums(gene_cnt>=3)>=3
cnt_expr <- gene_cnt[keep,]
seur_merged <- subset(seur_merged, features = rownames(cnt_expr))

# Pre-process merged object
seur_merged <- NormalizeData(seur_merged, normalization.method = 'LogNormalize', scale.factor = 10000)
seur_merged <- FindVariableFeatures(seur_merged, selection.method = 'vst', nfeatures = 3000)
feat <- VariableFeatures(seur_merged)
feat <- setdiff(feat, rp_gene)
feat <- setdiff(feat, mt_gene)
feat <- head(feat, n = 2000)
VariableFeatures(seur_merged) <- feat
seur_merged <- ScaleData(seur_merged)
seur_merged <- RunPCA(seur_merged, features = VariableFeatures(object = seur_merged))

# Run harmony
seur <- RunHarmony(seur_merged, group.by.vars = "BatchID", reduction = "pca", 
                   assay.use = "RNA", reduction.save = "harmony", max_iter = 50,
                   plot_convergence = TRUE, verbose = TRUE)
p.elbow <- ElbowPlot(seur, ndims = 50, reduction = 'harmony')

# Cluster data
seur <- FindNeighbors(seur, reduction = 'harmony', dims = 1:30)
seur <- FindClusters(seur, resolution = 0.5)
seur <- RunUMAP(seur, reduction = 'harmony', dims = 1:30)

p1 <- DimPlot(seur, reduction = "harmony", pt.size = .1, group.by = "BatchID", shuffle = TRUE)
ggsave(filename = "plot_harmony1_2_batch.pdf", plot = p1, 
       device = 'pdf', width = 7, height = 7, units = 'in')

p2 <- DimPlot(seur, reduction = 'umap', label = TRUE)
ggsave(filename = "plot_umap_cluster_res0.5.pdf", plot = p2,
       device = 'pdf', width = 7, height = 7, units = 'in')


### MAJOR CELL POPULATION ANNOTATION ###############################################################
# Read human SAT snRNA-seq data from Emont et al. Nature 2022
emont <- readRDS('human_all_sat_sce.rds')

# Convert Seurat object to SCE object for SingleR
seur.sce <- as.SingleCellExperiment(seur, assay = 'RNA')

# Annotate each cluster using SingleR with SAT snRNA-seq data from Emont et al. Nature 2022 data as ref
singler_cellty <- SingleR(
  test = seur.sce,
  ref = emont,
  labels = emont$cell_type2,
  method = "cluster",
  genes = 'de',
  cluster = seur.sce$RNA_snn_res.0.5,
  de.method = "wilcox",
  recompute = TRUE,
  fine.tune = TRUE,
  tune.thresh = 0.5,
  prune = TRUE,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  check.missing = TRUE,
  BPPARAM = MulticoreParam(8)
)
annot_cellty <- data.frame(
  'RNA_snn_res.0.5' = rownames(singler_cellty),
  'CellType_SingleR' = singler_cellty$pruned.labels
)

# Assign a major cell population to each cluster based on SingleR annotation results
annot_cellty <- mutate(
  annot_cellty,
  celltype_broad = c(
    'adipocyte', 'ASPC', 'myeloid', 'vascular', 'myeloid',
    'adipocyte', 'lymphoid', 'vascular', 'lymphoid', 'ASPC',
    'adipocyte', 'vascular', 'myeloid', 'myeloid', 'vascular',
    'vascular', 'vascular', 'lymphoid', 'myeloid', 'myeloid',
    'vascular', 'vascular'))

# Add assigned cell population annotation to Seurat object metadata
seur.meta <- seur@meta.data
seur.meta <- left_join(seur.meta, annot_cellty, by = 'RNA_snn_res.0.5')
rownames(seur.meta) <- rownames(seur@meta.data)
seur <- AddMetaData(seur, seur.meta)
saveRDS(seur, "sat_snrna_harmony.rds")

